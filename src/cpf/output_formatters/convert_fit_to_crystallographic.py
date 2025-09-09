__all__ = ["fourier_to_crystallographic"]

import numpy as np
from uncertainties import ufloat

from cpf.output_formatters.jcpds import jcpds
from cpf.IO_functions import peak_hkl
from cpf.output_formatters.crystallographic_operations import indicies4to3
from cpf.IO_functions import replace_null_terms


def fourier_to_crystallographic(
    coefficients,
    SampleGeometry: str = "3d",
    SampleDeformation: str = "compression",
    correlation_coeffs=None,
    subpattern=0,
    peak=0,
    debug=False,
    **kwargs,
):
    """
    Convert fourier coefficients into the centroid and differential values expected for crystallogrpahic strains/stresses.

    Currently this script uses the first order approximation that d0 = z-th order Foutier coefficient, differential strain is
    the length of the second order coefficients as a vector, and the angle is the arctan of the angle from these two coefficients.
    This is true if the center of the Debye-Scherer ring is at 0,0. If X0 or y0 are not zero a small ststematice error creeps into the
    coefficients. However, accoutning for these higher order terms and converting x0 and y0 into real units is non-trivial.

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
            "d_mean" -- mean d-spacing of diffraction ring, assuming 2d or 3d 'SampleGeometry'
            "differential" -- differential 
            "Q" -- differential strain
            "orientation" -- orientation of the strain field, under 'SampleDeformation'
            "d_max" -- maximum d-sapcing of fitted diffraction peak
            "d_min" -- minimum d-sapcing of fitted diffraction peak
            "x0" -- cartesian coordinate for centre of diffraction peak
            "y0" -- cartesian coordinate for centre of diffraction peak
        

    SampleGeometry options:
        - either 2d or 3d
    SampleDeformation options:
        - either 'compresion' or 'extension'.

        'compression' assumes that shortening strains are positive (correctly in my view)
        and that the differential strain is more compressive strain than in the other two directions.

        'extension' assumes that an extensive strain is positive (well, someone has to) and
        that the differentail strain is more extensive than the strains in the other two directions.


    """
    # FIX ME: strictly speaking all the errors here need to include the covarience matrix.
    # This is calculated and stored but not used here. If the values are small then the errors are about correct.
    # However, if the values are large then the errors are not correct.
    # the uncertainties package might be the package to use here.
    # https://uncertainties-python-package.readthedocs.io/en/latest/user_guide.html?highlight=covariance#covariance-matrix

    # %% validate the inputs.
    SampleGeometry = SampleGeometry.lower()
    SampleDeformation = SampleDeformation.lower()
    if not (SampleGeometry == "3d" or SampleGeometry == "2d"):
        err_str = "SampleGeometry is not recognised. It must be '2d' or '3d'."
        raise ValueError(err_str)
    # if SampleDeformation=="extension":
    #     err_str = "SampleDeformation 'extension' is not implemented."
    #     raise ValueError(err_str)
    if not (SampleDeformation == "compression" or SampleDeformation == "extension"):
        err_str = "SampleDeformation is not recognised. It must be 'compression' or 'extension'."
        raise ValueError(err_str)
    # FIX ME: another possitbility exists, where the orientation of the stress/strain is known and
    # the experiment oscillates between compression and extenion. This should perhaps be added to the
    # possibilities.

    if isinstance(coefficients, dict):
        coefficients = [coefficients]

    if not isinstance(coefficients, list):
        raise ValueError("The coefficients need to be a list of dictionaries.")

    # catch 'null' terms in fits
    coefficients = replace_null_terms(coefficients)

    # %% differential coefficients, errors and covarience

    # %%% angle
    # a = sin??
    # b = cos??
    # out_ang  = atan(b/a) / 2
    # d (atan(c))/dc = 1/(c^2+1). c = b/a. dc = c.((da/a)^2 + (db/b)^2)^(1/2)
    # out_angerr = dc.
    # FIX ME need to check this.
    if (
        len(coefficients[subpattern]["peak"][peak]["d-space"]) >= 5 
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
        ) / 2
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
    if len(coefficients[subpattern]["peak"][peak]["d-space"]) >= 5 and coefficients[subpattern]["peak"][peak]["d-space"][4] > 0:
        if coefficients[subpattern]["peak"][peak]["d-space"][3] <= 0:
            out_ang += np.pi / 2
        else:
            out_ang -= np.pi / 2
    # convert into degrees.
    out_ang = np.rad2deg(out_ang)
    out_angerr = np.rad2deg(out_angerr)

    # %%% differential strain
    # differentail (3d) = (a2^2+b2^2)^(1/2)
    if len(coefficients[subpattern]["peak"][peak]["d-space"]) >= 5:
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
        out_dd = np.nan
        out_dderr = np.nan
    
    #%%% d_max and d_min.
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

    # %%% d0  (centroid)
    if SampleGeometry == "2d":
        # d0 is the mean of the d-spacings. In this case it is the middle of the line.
        out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0]
        out_d0err = coefficients[subpattern]["peak"][peak]["d-space_err"][0]

        # organise for compression or extension experiment
        if SampleDeformation == "compression":
            pass
        elif SampleDeformation == "extension":
            pass  # out_ang += 90

    elif SampleGeometry == "3d" and SampleDeformation == "compression":
        # 1/3 of the way from the middle to the maximum d-spacing (miniminm strain)
        out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0] + out_dd / 3
        out_d0err = (
            coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2
            + (out_dderr / 3) ** 2
        ) ** (1 / 2)

    elif SampleGeometry == "3d" and SampleDeformation == "extension":
        # 1/3 of the way from the middle to the maximum d-spacing
        out_d0 = coefficients[subpattern]["peak"][peak]["d-space"][0] - out_dd / 3
        out_d0err = (
            coefficients[subpattern]["peak"][peak]["d-space_err"][0] ** 2
            + (out_dderr / 3) ** 2
        ) ** (1 / 2)

        # move the angle round to account for fact it is extension
        # out_ang += np.pi/2

    else:
        # issue an error
        err_str = "The geometry type is not recognised. Check input file."
        raise ValueError(err_str)

    # %% Q from Singh et al (1998).
    # This is a third of the differential strain.
    # It is defined here directly for potentail convenience.
    # out_Q = out_dd / out_d0
    if SampleGeometry == "2d":
        scale = 1
    elif SampleGeometry == "3d":
        scale = 3 / 2

    # if SampleDeformation=="extension":
    #     scale *= -1

    out_Q = out_dd / out_d0 / scale
    out_Qerr = out_dderr / scale

    # reoridentate extensions back to the right way.
    if SampleDeformation == "extension":
        if out_ang >= 0:
            out_ang -= 90
        else:
            out_ang += 90

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
    correlation_coeffs=None,
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
        DESCRIPTION.

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
    if not (SampleGeometry == "3d" or SampleGeometry == "2d"):
        err_str = "SampleGeometry is not recognised. It must be '2d' or '3d'."
        raise ValueError(err_str)
    # if SampleDeformation=="extension":
    #     err_str = "SampleDeformation 'extension' is not implemented."
    #     raise ValueError(err_str)
    if not (SampleDeformation == "compression" or SampleDeformation == "extension"):
        err_str = "SampleDeformation is not recognised. It must be 'compression' or 'extension'."
        raise ValueError(err_str)
    # FIX ME: another possitbility exists, where the orientation of the stress/strain is known and
    # the experiment oscillates between compression and extenion. This should perhaps be added to the
    # possibilities.

    if isinstance(coefficients, dict):
        coefficients = [coefficients]

    if not isinstance(coefficients, list):
        raise ValueError("The coefficients need to be a list of dictionaries.")

    # catch 'null' terms in fits
    coefficients = replace_null_terms(coefficients)

    # get or guess phase
    if phase is None:
        #list all phases in fits
        phases = []
        for i in range(len(coefficients)):
            for j in range(len(coefficients[i]["peak"])):
                if "phase" in coefficients[i]["peak"][j]:
                    phases.append(coefficients[i]["peak"][j]["phase"])
        phase = max(set(phases), key=phases.count)  
    
    #get or guess jcpds
    


    # %% get d0 accounting for difrerential strain and sample geometry
    for i in range(len(coefficients)):
        for j in range(len(coefficients[i]["peak"])):
        
            crystallographic = fourier_to_crystallographic(coefficients,
                            SampleGeometry,
                            SampleDeformation,
                            correlation_coeffs,
                            subpattern=i,
                            peak=j,
                            debug=debug,
                            **kwargs)

            coefficients[i]["peak"][j]["cryst_prop"] = crystallographic

    # intial guess (a0, b0, c0 etc) for lattice parameters comes from jcpds file
    # solve for unit cell    
    jcpds_obj = jcpds()
    jcpds_obj.load_file(jcpds_file)
    #clear reflections from jcpds
    jcpds_obj.remove_reflection("all")
    
    # add reflections for unit cell we need to fit.
    for i in range(len(coefficients)):
        for j in range(len(coefficients[1]["peak"])):
            if coefficients[i]["peak"][j]["phase"] == phase:
                hkl = peak_hkl(coefficients[i], j, string=False)[0]
                if len(hkl) == 4:
                    # convert to 3 value Miller indicies
                    hkl = indicies4to3(hkl)
                jcpds_obj.add_reflection(h=hkl[0], k=hkl[1], l=hkl[2],
                                 dobs = coefficients[i]["peak"][j]["cryst_prop"]["dp"])
    jcpds_obj.compute_d0() # compute lattice parameters for unit cell from jcpds
    jcpds_obj.fit_lattice_parameters()

    uc_parts = jcpds_obj.get_unique_unitcell_params()
    uc_parms = {}
    for ind in range(len(uc_parts)):
        uc_parms[uc_parts[ind]] = getattr(jcpds_obj, uc_parts[ind])    
    uc_parms["volume"] = jcpds_obj.v        

    return uc_parms
