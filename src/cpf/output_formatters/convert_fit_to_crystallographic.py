__all__ = ["fourier_to_crystallographic"]

import numpy as np

from cpf.IO_functions import replace_null_terms


def fourier_to_crystallographic(
    coefficients,
    SampleGeometry="3d",
    SampleDeformation="compression",
    correlation_coeffs=None,
    subpattern=0,
    peak=0,
    debug=False,
    **kwargs,
):
    """
    Convert fourier coefficients into the centroid and differnetial values expected for crystallogrpahic strains/stresses.

    Currently this script uses the forst order approximation that d0 = z-th order Foutier coefficient, differential strain is
    the length of the second order coefficients as a vector, and the angle is the arctan of the angle from these two coefficients.
    This is true if the center of the Debye-Scherer ring is at 0,0. If X0 or y0 are not zero a small ststematice error creeps into the
    coefficients. However, accoutning for these higher order terms and converting x0 and y0 into real units is non-trivial.

    Parameters
    ----------
    coefficients : TYPE
        DESCRIPTION.
    SampleGeometry : String, optional
        DESCRIPTION. The default is "3d-compression".
    correlation_coeffs : TYPE, optional
        DESCRIPTION. The default is None.
    debug : TYPE, optional
        DESCRIPTION. The default is False.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    crystallographic_properties.

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
        coefficients[subpattern]["peak"][peak]["d-space"][4] != 0
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
    elif coefficients[subpattern]["peak"][peak]["d-space"][3] != 0:
        out_ang = np.pi / 2
        out_angerr = 0
    else:
        out_ang = np.nan
        out_angerr = np.nan
        # FIXME: this is a bodged fix for now. It needs to be calcualted assuming the error is not also zero.
    # correction to make angle correct (otherwise potentially out by pi/2)
    if coefficients[subpattern]["peak"][peak]["d-space"][4] > 0:
        if coefficients[subpattern]["peak"][peak]["d-space"][3] <= 0:
            out_ang += np.pi / 2
        else:
            out_ang -= np.pi / 2
    # convert into degrees.
    out_ang = np.rad2deg(out_ang)
    out_angerr = np.rad2deg(out_angerr)

    # %%% differential strain
    # differentail (3d) = (a2^2+b2^2)^(1/2)
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

    differential_coefficients["dp"] = out_d0
    differential_coefficients["dp_err"] = out_d0err

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
