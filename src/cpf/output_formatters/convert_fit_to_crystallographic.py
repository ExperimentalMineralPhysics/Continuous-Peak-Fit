__all__ = ["fourier_to_crystallographic"]

import glob
import json
import re

import numpy as np
import pandas as pd

from uncertainties import ufloat
import uncertainties.umath as umath

from cpf.output_formatters.crystallographic_operations import plane_indices_4_to_3
from cpf.series_functions import coefficient_expand, series_properties
# from uncertainties import ufloat
from cpf.output_formatters.jcpds import jcpds
from cpf.util.io import make_outfile_name, peak_hkl, replace_value
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.convert_fit_to_crystallographic")


def _sample_geometries(as_list=True):
    """
    Data structure of recognised sample geometries
    
    Parameters
    ----------
    all : bool, default is True
        if true returns dictionary of recognised geometries where the values are their intermal labelling
        if False returns list of recognised geometies.         
    
    Returns
    -------
    geometries : dict, list
    
    """
    geometries =  {"2d": "2d",
           "vertical": "vertical", 
           "horizontal": "horizontal", 
           "n-s": "vertical", 
           "e-w": "horizontal", 
           "nw-se": "135", 
           "ne-sw": "45", 
           "shear": "2d", 
           "compression": "compression",  
           "extension": "extension", 
           "3d": "compression",
           "axial": "2d"}
    
    if as_list:
        return list(geometries.keys())
    else:
        return geometries


def fourier_to_crystallographic(
    coefficients,
    SampleGeometry: str = "3d",
    SampleStrain: str = "CompressionPositive",
    correlation_coeffs=None,
    subpattern=0,
    peak=0,
    debug=False,
    **kwargs,
):
    """
    Convert fourier coefficients into the centroid and differential values expected for crystallogrpahic strains given sample geometry.

    Currently this script uses the first order approximation that d0 = zeroth order Foutier coefficient, differential strain is
    determined from the second order coefficients and the angle is the arctan of the angle from these two coefficients.
    
    This is true if the center of the Debye-Scherer ring is at 0,0. If X0 or y0 are not zero a small systematic error creeps into the
    coefficients. However, accouting for these higher order terms and converting x0 and y0 into real units is non-trivial.

    Parameters
    ----------
    coefficients : dict
        Coeffecient dictionary used in cpf.
    SampleGeometry : str, optional
        Geometry of the sample - 2D or 3D stress/stress field. The default is "3d".
    SampleStrain : str, optional
        Definition of strain for formation of Q. The default is "CompressionPositived".
    correlation_coeffs : dict, optional
        correlation coefficient dictionary as created by lmfit. The default is None.
    subpattern : list, int, optional
        Which subpattern in the coefficients to calulcate parameters for. The default is 0.
    peak : int, optional
        Which peak in the subpattern to calulcate parameters for. The default is 0.
    **kwargs : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    differential_coefficients : dict
        Dictionary of the calculated properties and their errors. The propertues calculated from the
        Fourier d-spacing coefficients are:
            "d_mean" -- mean d-spacing of diffraction ring, given the 'SampleGeometry'
            "differential" -- total differential elastic strain (d_max - d_min)
            "Q" -- lattice strain parameter
            "orientation" -- orientation of the strain field
            "d_max" -- maximum d-sapcing of fitted diffraction peak
            "d_min" -- minimum d-sapcing of fitted diffraction peak
            "x0" -- cartesian coordinate for centre of diffraction peak
            "y0" -- cartesian coordinate for centre of diffraction peak
        the error on each value has *err appended to the name

    SampleGeometry options (indepedent of strain-definition):
        - "2d"
        - "compression" == assumes 3rd (unmeasured) cartesian axis is same as maximum d-spacings (i.e. minimum is unique)
        - "extension" == assumes 3rd (unmeasured) cartesian axis is same as minimum d-spacings (i.e. maximum is unique)
        - "vertical"  == assumes 90def is vertical. Returns angle within +/-45 of 90 and associated Q for this angle (+/- as defined by SampleStrain)
        - "horizontal"  == assumes 0def is horizontal. Returns angle within +/-45 of 90 and associated Q for this angle (+/- as defined by SampleStrain)
        - "3d" == same as "vertical"
        - "shear" == same as "2d"
        - "axial" == same as "2d"
        - "n-s", "e-w", "nw-se", "ne-sw" -- Returns angle within +/-45 of cardinal direction and associated Q for this angle (+/- as defined by SampleStrain). assumes 90 vertical = N / 0 horizontal = E. 
        go get list of Sample Geometries and their maping call:
            cpf.output_formatters.convert_fit_to_crystallographic._sample_geometries(as_list=False)
            
    SampleStrain:
        - either "CompressionPositive" or "ExtensionPositive"
        This does not affect the geometry calculations but does affect the value of Q returned.

    """
    # FIX ME: strictly speaking all the errors here need to include the covarience matrix.
    # This is calculated and stored but not used here. If the values are small then the errors are about correct.
    # However, if the values are large then the errors are not correct.
    # the uncertainties package might be the package to use here.
    # https://uncertainties-python-package.readthedocs.io/en/latest/user_guide.html?highlight=covariance#covariance-matrix

    # %% validate the inputs.
    SampleGeometry = SampleGeometry.lower()
    SampleStrain = SampleStrain.lower()
    
    #geometry
    geometries = _sample_geometries(as_list=False)# override to use old values. Dangerous. needs removing.
    # FIXME: remove below in future release
    if "SampleDeformation" in list(kwargs):
        logger.warning("'SampleDeformation' is old nomenclature and will be removed in future release")
        if "SampleDeformation" == "compression":
            SampleStrain = "CompressionPositive"
        elif "SampleDeformation" == "extension":
            SampleStrain = "CompressionNegative"
    # FIXME: end ==================
    if SampleGeometry not in list(geometries.keys()):
        err_str = f"SampleGeometry '{SampleGeometry}' is not recognised. It must be one of: {list(geometries.keys())}."
        raise ValueError(err_str)
    # collapse to internal names. 
    SampleGeometry = geometries[SampleGeometry].lower()
    
    """ set prefectors for positive or negative compression.
    if compression positive:
        \epsilon = (d_p - d_m)/d_p
        d_m = d_p [1 - Q(1-3cos^2\psi)]
    if compression negative / positive expansion:
        \epsilon = (d_m - d_p)/d_p
        d_m = d_p [1 + Q(1-3cos^2\psi)]
        
    cast with a prefactor:
        d_m = d_p [1 + A* Q(1-3cos^2\psi)]
    where
        A = -1 for CompressionPositive
        A = 1  for ExtensionPositive
    """
    if SampleStrain.lower() == "CompressionPositive".lower():
        strain_factor = -1
    elif SampleStrain.lower() == "ExtensionPositive".lower():
        strain_factor = 1
    else:
        err_str = "SampleStrain is not recognised. It must be 'CompressionPositive' or 'ExtensionPositive'."
        raise ValueError(err_str)
    
    # sort coefficents 
    if isinstance(coefficients, dict):
        coefficients = [coefficients]
    if not isinstance(coefficients, list):
        raise ValueError("The coefficients need to be a list of dictionaries.")
    # catch 'null' terms in fits
    coefficients = replace_value(coefficients, old=None, new=np.nan)
    # catch d-spacing that is too short for 3D geometry to work.
    if len(coefficients[subpattern]["peak"][peak]["d-space"]) <= 3:
        SampleGeometry == "2d"

    # %% differential coefficients, errors and covarience
    # %%% angle
    # a = sin??
    # b = cos??
    # out_ang  = atan(b/a) / 2
    # d (atan(c))/dc = 1/(c^2+1). c = b/a. dc = c.((da/a)^2 + (db/b)^2)^(1/2)
    # out_angerr = dc.
    # FIX ME need to check this.
    if len(coefficients[subpattern]["peak"][peak]["d-space"]) <= 3:
        out_ang = 0
        out_angerr = 0
    elif (
        len(coefficients[subpattern]["peak"][peak]["d-space"]) >= 4
        and coefficients[subpattern]["peak"][peak]["d-space"][4] != 0
        and coefficients[subpattern]["peak"][peak]["d-space"][3] != 0
    ):
        out_ang = (
            np.arctan(
                coefficients[subpattern]["peak"][peak]["d-space"][3]
                / coefficients[subpattern]["peak"][peak]["d-space"][4]
            )
            / 2
        )
        out_angerr = (
            1
            / (
                (
                    coefficients[subpattern]["peak"][peak]["d-space"][3]
                    / coefficients[subpattern]["peak"][peak]["d-space"][4]
                )
                ** 2
                + 1
            )
            * (
                np.abs(out_ang)
                * (
                    (
                        coefficients[subpattern]["peak"][peak]["d-space_err"][3]
                        / coefficients[subpattern]["peak"][peak]["d-space"][3]
                    )
                    ** 2
                    + (
                        coefficients[subpattern]["peak"][peak]["d-space_err"][4]
                        / coefficients[subpattern]["peak"][peak]["d-space"][4]
                    )
                    ** 2
                )
                ** (1 / 2)
            )
            / 2
        )
    elif (
        len(coefficients[subpattern]["peak"][peak]["d-space"]) >= 5
        and coefficients[subpattern]["peak"][peak]["d-space"][3] != 0
    ):
        out_ang = np.pi / 2
        out_angerr = 0
    else:
        out_ang = np.nan
        out_angerr = np.nan
        # FIXME: this is a bodged fix for now. It needs to be calcualted assuming the error is not also zero.
    # correction to make angle correct (otherwise potentially out by pi/2)
    if (
        len(coefficients[subpattern]["peak"][peak]["d-space"]) >= 4
        and coefficients[subpattern]["peak"][peak]["d-space"][4] > 0
    ):
        if coefficients[subpattern]["peak"][peak]["d-space"][3] <= 0:
            out_ang += np.pi / 2
        else:
            out_ang -= np.pi / 2
    # convert into degrees.
    out_ang = np.rad2deg(out_ang)
    out_angerr = np.rad2deg(out_angerr)

    # deal with sin/cos wrapping
    if (
        len(coefficients[subpattern]["peak"][peak]["d-space"]) >= 4
        and coefficients[subpattern]["peak"][peak]["d-space"][3] > 0
    ):
        out_ang += 180

    # %%% differential strain
    # differentail (3d) = (a2^2+b2^2)^(1/2)
    if len(coefficients[subpattern]["peak"][peak]["d-space"]) <= 3:
        out_dd = 0
        out_dderr = 0
    elif len(coefficients[subpattern]["peak"][peak]["d-space"]) >= 5:
        out_dd = np.sqrt(
            coefficients[subpattern]["peak"][peak]["d-space"][3] ** 2
            + coefficients[subpattern]["peak"][peak]["d-space"][4] ** 2
        )
        # out_dderr= [(2.a.da.)^2 + (2.b.db)^2]^(1/2)]^(1/2)
        out_dderr = (
            (
                2
                * coefficients[subpattern]["peak"][peak]["d-space"][3]
                * coefficients[subpattern]["peak"][peak]["d-space_err"][3]
            )
            ** 2
            + (
                2
                * coefficients[subpattern]["peak"][peak]["d-space"][4]
                * coefficients[subpattern]["peak"][peak]["d-space_err"][4]
            )
            ** 2
        ) ** (1 / 4)
    else:
        out_dd = 0
        out_dderr = np.nan

    # %%% d_max and d_min.
    # differential max
    out_dmax = coefficients[subpattern]["peak"][peak]["d-space"][0] + out_dd
    out_dmaxerr = (
        coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2 + out_dderr**2
    ) ** (1 / 2)
    # differential min
    out_dmin = coefficients[subpattern]["peak"][peak]["d-space"][0] - out_dd
    out_dminerr = (
        coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2 + out_dderr**2
    ) ** (1 / 2)

        
    # %%% get centroids (d0), Q etc. based on geometry. 
    
    #all geometries not specified as 2d are assuemd to be 3d. The determine is then just 
    # working out what is the inmeased axis. 
    match SampleGeometry:
        case "2d":
            # d0 is the mean of the d-spacings. In this case it is the middle of the line.
            out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0]
            out_d0err = coefficients[subpattern]["peak"][peak]["d-space_err"][0]

        case "compression":
            # assumes third cartesian axis is larger than mean of diffraction ring
            # 1/3 of the way from the middle to the maximum d-spacing
            out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0] + -strain_factor * out_dd / 3
            out_d0err = (
                coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2
                + (out_dderr / 3) ** 2
            ) ** (1 / 2)
        
        case "extension":
            # assumes third cartesian axis is smaller than mean of diffraction ring
            # 1/3 of the way from the middle to the minimum d-spacing
            out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0] - -strain_factor * out_dd / 3
            out_d0err = (
                coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2
                + (out_dderr / 3) ** 2
            ) ** (1 / 2)
        
            out_dd = -out_dd
        
        case "vertical":
            # assumes third cartesian axis same as horizontal axis
            
            # organise to have the correct signs (+/-)
            out_dd *= (np.round(out_ang/90) % 2 - 0.5)*2
            out_ang = ((np.round(out_ang/90) % 2 -1)*-90 + out_ang)%180
            
            out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0] + out_dd / 3
            out_d0err = (
                coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2
                + (out_dderr / 3) ** 2
            ) ** (1 / 2)
                
        case "horizontal":
            # assumes third cartesian axis same as vertical axis
            
            #need to know which way the strains are 
            d_0 = coefficient_expand(
                np.array([0]),
                param=coefficients[subpattern]["peak"][peak]["d-space"],
                coeff_type="fourier",
                comp_str="d",
            )
            d_90 = coefficient_expand(
                np.array([90]),
                param=coefficients[subpattern]["peak"][peak]["d-space"],
                coeff_type="fourier",
                comp_str="d",
            )
            if d_0 > d_90:
                scl = 1
            else:
                scl = -1
                
            out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0] - scl*out_dd / 3
            out_d0err = (
                coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2
                + (out_dderr / 3) ** 2
            ) ** (1 / 2)

        case "45":
            # assumes third cartesian axis same as 135+/-45 axis
            
            #need to know which way the strains are 
            d_0 = coefficient_expand(
                np.array([45]),
                param=coefficients[subpattern]["peak"][peak]["d-space"],
                coeff_type="fourier",
                comp_str="d",
            )
            d_90 = coefficient_expand(
                np.array([135]),
                param=coefficients[subpattern]["peak"][peak]["d-space"],
                coeff_type="fourier",
                comp_str="d",
            )
            if d_0 > d_90:
                scl = 1
            else:
                scl = -1
                
            out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0] - scl*out_dd / 3
            out_d0err = (
                coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2
                + (out_dderr / 3) ** 2
            ) ** (1 / 2)
            
        case "135":
            # assumes third cartesian axis same as 45+/-45 axis
            
            #need to know which way the strains are 
            d_0 = coefficient_expand(
                np.array([45]),
                param=coefficients[subpattern]["peak"][peak]["d-space"],
                coeff_type="fourier",
                comp_str="d",
            )
            d_90 = coefficient_expand(
                np.array([135]),
                param=coefficients[subpattern]["peak"][peak]["d-space"],
                coeff_type="fourier",
                comp_str="d",
            )
            if d_0 > d_90:
                scl = 1
            else:
                scl = -1
                
            out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0] + scl*out_dd / 3
            out_d0err = (
                coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2
                + (out_dderr / 3) ** 2
            ) ** (1 / 2)

        case _:
            # issue an error
            err_str = f"SampleGeometry is not recognised. It must be one of: {list(geometries.keys())}."
            raise ValueError(err_str)

    
    # %% Q
    if SampleGeometry == "2d":
        scale = 2 / 2
    else:
        scale = 3 / 2
        
    if out_d0 != 0:
        out_Q = out_dd / out_d0 / scale
        out_Qerr = out_dderr / scale
    else:
        # then there is nothing to divide by and will get error
        out_Q = np.nan
        out_Qerr = np.nan
        
    
    # %% Organise values and angles depending on geometry specified. 
    match SampleGeometry:
        case "2d":
            pass
            
        case "compression":
            # orientation aligned with smallest d-spacing 
            pass
        
        case "extension":
            # orientation aligned with largest d-spacing 
            pass
        
        case "vertical":
            # orientation aligned 90 +/-45
            # change Q to match
            
            # correct Q for the strain definition. 
            if strain_factor != -1: # compression negative/extension positive
                out_Q = -out_Q
            else:
                pass
            # if coefficients[subpattern]["peak"][peak]["d-space"][4]<0:
            # out_ang -= (np.round(out_ang/90) ) *90 - 90 # this one drops into range
            # out_Q *=  (np.round(out_ang/90) % 2 - 0.5)*2
        
        case "horizontal":
            # orientation aligned 0 +/-45
            # change Q to match
            if coefficients[subpattern]["peak"][peak]["d-space"][4]>0:
                out_Q = -out_Q
            out_ang -= (np.round(out_ang/90) ) *90  # this one drops into range
            out_Q *=  -(np.round(out_ang/90) % 2 - 0.5)*2

            # correct Q for the strain definition. 
            if strain_factor != -1: # compression negative/extension positive
                out_Q = -out_Q
            else:
                pass


        case "45":
            # orientation aligned 45 +/-45
            # change Q to match
            if coefficients[subpattern]["peak"][peak]["d-space"][3]>0:
                out_Q = -out_Q
            out_ang -= (np.round((out_ang-45)/90) ) *90  # this one drops into range
            out_Q *=  (np.round((out_ang+45)/90) % 2 - 0.5)*2

            # correct Q for the strain definition. 
            if strain_factor != -1: # compression negative/extension positive
                out_Q = -out_Q
            else:
                pass
        case "135":
            # orientation aligned 135 +/-45
            # change Q to match
            if coefficients[subpattern]["peak"][peak]["d-space"][3]<0:
                out_Q = -out_Q
            out_ang -= (np.round((out_ang-45)/90) ) *90 + 90 # this one drops into range
            out_Q *=  (np.round((out_ang+135)/90) % 2 - 0.5)*2
            
            
            # correct Q for the strain definition. 
            if strain_factor != -1: # compression negative/extension positive
                out_Q = -out_Q
            else:
                pass
        
        case _:
            # issue an error
            err_str = f"SampleGeometry is not recognised. It must be one of: {list(geometries.keys())}."
            raise ValueError(err_str)


    # out_dd needs to be differential not deviatoric or something else.
    out_dd *= 2
    out_dderr *= 2


    # %%% x0 and y0
    # values is in d-spacing and needs converting to mm via calibration.
    if not debug:
        # FIX ME: this conversion needs a calibration and as far as I can work out is non-trivial.
        out_x0 = np.nan
        out_x0err = np.nan
        out_y0 = np.nan
        out_y0err = np.nan
    else:
        out_x0 = coefficients[subpattern]["peak"][peak]["d-space"][2]
        out_x0err = coefficients[subpattern]["peak"][peak]["d-space_err"][2]
        out_y0 = coefficients[subpattern]["peak"][peak]["d-space"][1]
        out_y0err = coefficients[subpattern]["peak"][peak]["d-space_err"][1]

    # %% collate values for output

    # catch 'null' as an error
    if out_dderr is None:
        out_dderr = np.nan
    if out_dmaxerr is None:
        out_dmaxerr = np.nan
    if out_dminerr is None:
        out_dminerr = np.nan
    if out_d0err is None:
        out_d0err = np.nan
    if out_x0err is None:
        out_x0err = np.nan
    if out_y0err is None:
        out_y0err = np.nan

    differential_coefficients = {}

    # FIX ME: this should be removed but remains becuase some outputs depend on it.
    differential_coefficients["dp"] = out_d0
    differential_coefficients["dp_err"] = out_d0err

    differential_coefficients["d_mean"] = out_d0
    differential_coefficients["d_mean_err"] = out_d0err

    differential_coefficients["differential"] = out_dd
    differential_coefficients["differential_err"] = out_dderr
    differential_coefficients["Q"] = out_Q
    differential_coefficients["Q_err"] = out_Qerr

    differential_coefficients["orientation"] = out_ang
    differential_coefficients["orientation_err"] = out_angerr

    differential_coefficients["d_max"] = out_dmax
    differential_coefficients["d_max_err"] = out_dmaxerr
    differential_coefficients["d_min"] = out_dmin
    differential_coefficients["d_min_err"] = out_dminerr

    differential_coefficients["x0"] = out_x0
    differential_coefficients["x0_err"] = out_x0err
    differential_coefficients["y0"] = out_y0
    differential_coefficients["y0_err"] = out_y0err

    return differential_coefficients


def fourier_to_unitcellvolume(
    coefficients,
    jcpds_file: None,
    phase: None,
    SampleGeometry: str = "3d",
    SampleDeformation: str = "compression",
    pressure=False,
    reflections_to_use="all",
    correlation_coeffs=None,
    weighted=True,
    temperature=np.nan,
    debug=False,
    **kwargs,
):
    """
    Convert fitted peak centers into unit cell parameters.

    The method builds a jcpds object and then fits the unit cell to the given peaks.

    FIXME: no errors are currently used or returned.

    Parameters
    ----------
    coefficients : dict
        Coeffecient dictionary used in cpf.
    SampleGeometry : str, optional
        Geometry of the sample - 2D or 3D stress/stress field. The default is "3d".
    SampleDeformation : str, optional
        "compression" or "extension" - type of deformation for the experiment. Determines orientation value.
        The default is "compression".
    correlation_coeffs : dict, optional
        correlation coefficient dictionary as created by lmfit. The default is None.
    subpattern : list, int, optional
        Which subpattern in the coefficients to calulcate parameters for. The default is 0.
    peak : int, optional
        Which peak in the subpattern to calulcate parameters for. The default is 0.
    **kwargs : TYPE
        See below.

    :Keyword Arguments:
        * *reflections_to_use* list  or "all"
          Which ranges from the list of fitted ranges to use.
          FIXME: does not work for multiple peaks in the same region

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    unitcell_volumes : dict
        Dictionary of the calculated cell parameters and the errors.


    SampleGeometry options:
        - either 2d or 3d
    SampleDeformation options:
        - either 'compresion' or 'extension'.

        'compression' assumes that shortening strains are positive (correctly in my view)
        and that the differential strain is more compressive strain than in the other two directions.

        'extension' assumes that an extensive strain is positive (well, someone has to) and
        that the differentail strain is more extensive than the strains in the other two directions.


    """

    # %% validate the inputs.
    SampleGeometry = SampleGeometry.lower()
    SampleDeformation = SampleDeformation.lower()
    if not SampleGeometry in _sample_geometries(as_list=True):# == "3d" or SampleGeometry == "2d"):
        err_str = f"SampleGeometry '{SampleGeometry}' is not recognised. It must be one of: {_sample_geometries(as_list=True)}."
        raise ValueError(err_str)

    if not (SampleDeformation == "compression" or SampleDeformation == "extension"):
        err_str = "SampleDeformation is not recognised. It must be 'compression' or 'extension'."
        raise ValueError(err_str)

    if isinstance(coefficients, dict):
        coefficients = [coefficients]
    if not isinstance(coefficients, list):
        raise ValueError("The coefficients need to be a list of dictionaries.")

    # flatten the multipeak parts of the coefficient structure.
    flat_coef = []
    for i in coefficients:
        for j in i["peak"]:
            flat_coef.append({"peak": [j]})

    if reflections_to_use == "all":
        reflections_to_use = list(range(len(flat_coef)))

    # catch 'null' terms in fits
    flat_coef = replace_value(flat_coef, old=None, new=np.nan)

    # get or guess phase
    if phase is None:
        # list all phases in fits
        phases = []
        for i in range(len(flat_coef)):
            for j in range(len(flat_coef[i]["peak"])):
                if "phase" in flat_coef[i]["peak"][j]:
                    phases.append(flat_coef[i]["peak"][j]["phase"])
        phase = max(set(phases), key=phases.count)

    # %% get d0 accounting for difrerential strain and sample geometry
    for i in reflections_to_use:
        crystallographic = fourier_to_crystallographic(
            flat_coef,
            SampleGeometry,
            SampleDeformation,
            correlation_coeffs,
            subpattern=i,
            peak=0,  # in flattened structure always the first peak
            debug=debug,
            **kwargs,
        )
        flat_coef[i]["peak"][0]["cryst_prop"] = crystallographic

    # intial guess (a0, b0, c0 etc) for lattice parameters comes from jcpds file
    # solve for unit cell
    jcpds_obj = jcpds()
    jcpds_obj.load_file(jcpds_file)
    # clear reflections from jcpds
    jcpds_obj.remove_reflection("all")

    # add reflections for unit cell we need to fit.
    hkls = []
    for i in reflections_to_use:
        j = 0  # always the first peak in flattened structure.
        # for j in range(len(flat_coef[1]["peak"])):
        if flat_coef[i]["peak"][j]["phase"] == phase:
            hkl = peak_hkl(flat_coef[i], j, as_string=False)[0]
            if len(hkl) == 4:
                # convert to 3 value Miller indicies
                hkl = plane_indices_4_to_3(hkl)
            jcpds_obj.add_reflection(
                h=hkl[0],
                k=hkl[1],
                l=hkl[2],
                dobs=flat_coef[i]["peak"][j]["cryst_prop"]["dp"],
                dobs_err=flat_coef[i]["peak"][j]["cryst_prop"]["dp_err"],
                delta_dobs=flat_coef[i]["peak"][j]["cryst_prop"]["Q"],
                delta_dobs_err=flat_coef[i]["peak"][j]["cryst_prop"]["Q_err"],
                orientationobs=flat_coef[i]["peak"][j]["cryst_prop"]["orientation"],
                orientationobs_err=flat_coef[i]["peak"][j]["cryst_prop"][
                    "orientation_err"
                ],
            )
            hkls.append(peak_hkl(flat_coef[i], j, as_string=True)[0])

    jcpds_obj.compute_d0()  # compute lattice parameters for unit cell from jcpds, otherwise initiation not complete.

    jcpds_obj.temperature = temperature

    returned = jcpds_obj.fit_lattice_parameters(weighted=weighted)

    uc_parts = jcpds_obj.get_unique_unitcell_params()
    uc_parms = {}

    # force observed value first in order
    if not np.isnan(temperature) and pressure:
        uc_parms["temperature"] = jcpds_obj.temperature

    for ind in range(len(uc_parts)):
        uc_parms[uc_parts[ind]] = getattr(jcpds_obj, uc_parts[ind]).nominal_value
        uc_parms[uc_parts[ind] + "_err"] = getattr(jcpds_obj, uc_parts[ind]).std_dev
    uc_parms["volume"] = jcpds_obj.v.nominal_value
    uc_parms["volume_err"] = jcpds_obj.v.std_dev

    if not np.isnan(temperature) and pressure:
        uc_parms["pressure"] = jcpds_obj.pressure.nominal_value
        uc_parms["pressure_err"] = jcpds_obj.pressure.std_dev
        for m, n, o in zip(
            jcpds_obj.get_reflection_stresses(),
            jcpds_obj.get_reflection_orientations(),
            jcpds_obj.get_reflection_hkls(),
        ):
            o = f"({o[0]}{o[1]}{o[2]})"
            uc_parms[o + " stress"] = m.nominal_value
            uc_parms[o + " stress_err"] = m.std_dev
            uc_parms[o + " angle"] = n.nominal_value
            uc_parms[o + " angle_err"] = n.std_dev

    # uc_parms["residuals"] = np.array2string(returned.residual, separator="; ", max_line_width=np.inf)
    if weighted == True:
        uc_parms["residuals_weighted"] = returned.residual
        rsd = []
        for i in jcpds_obj.reflections:
            rsd.append(i.dobs.nominal_value - i.d)
        uc_parms["residuals_unweighted"] = np.array(rsd)
    else:
        uc_parms["residuals_unweighted"] = returned.residual
    uc_parms["hkl"] = hkls
    return uc_parms
