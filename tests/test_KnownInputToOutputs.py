import unittest

import cpf.lmfit_model as lmm
import cpf.output_formatters.convert_fit_to_crystallographic as cfc
import cpf.series_functions as sf
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Parameters

# class TestCreateNewparams(unittest.TestCase):
#    def setUp(self):
#        pass#0.071*2
#
#    def test_FitFourierToDspacing(self):
"""
This is somewhere between a test case and a proove that the calcualtions are correct. 

The script takes a number of randomly generated diffraction ring centers, fits a Fourier series to them and 
tries to replicate the origininal values again. 

It works for the centroid (d0), Q (Singh et al 1998's differential stress/strain factor) and its orientation.
I cannot recover the x0 and y0 at the time of writing. These are a complex function of the parameters 
which I don't know
"""


debug = False


# set up problem
SampleGeometry = 3  # dimensions
trig_orders = 2

wavelength = 0.1234  # angstroms
distance = 300  # mm

opts = ["extension", "compression"]

# arrays to stash variables
angles = []
second_order_corefs = []
xs = []
ys = []
d0s = []
dmaxs = []
dmins = []
Qs = []
chsqs = []
all_params = []

# set of parameters that defind the d-spacing location with azimuth.
num_peaks = 1
if 0:
    # run the problem once for set parameters
    d0 = 1  # Angstrons
    Q = 0.02  # as per Singh et al., 1998
    orientation = [60.5]  # np.linspace(-90, 89, 90) #degrees
    x0 = -2.0  # 1E-15#0.000 #mm
    y0 = 0.5  # 10 #mm

    d0 = np.array([d0])
    Q = np.array([Q])
    orientation = np.array(orientation)
    x0 = np.array([x0])
    y0 = np.array([y0])

else:
    # run the problem lots of times.
    num_ang = 500
    small = 1
    d0 = np.random.default_rng().uniform(0.5, 7, num_ang)  # Angstrons
    Q = np.random.default_rng().uniform(
        0, 0.4, num_ang
    )  # *d0 # as per Singh et al., 1998
    orientation = np.random.default_rng().uniform(-90, 90, num_ang)  # degrees
    x0 = np.random.default_rng().uniform(-small, small, num_ang)  # mm
    y0 = np.random.default_rng().uniform(-small, small, num_ang)  # mm

    SampleGeometry = (2 + np.round(np.random.uniform(size=num_ang))).astype(
        int
    )  # dimensions
    SampleDeformation = np.random.choice(opts, num_ang)

    d0 = np.array(d0)
    Q = np.array(Q)
    orientation = np.array(orientation)
    x0 = np.array(x0)
    y0 = np.array(y0)


# make a set of aximuths for the Fourier data to work around.
azm = np.linspace(0, 360, 3601) + 0.1


for i in range(len(orientation)):
    ## make a set of centroid positions.

    # extension is this nomenclature is negative. therefore Q is negative
    pos_neg = -1 if np.bool_(SampleDeformation[i] == "extension") == True else 1

    # Calcualte d-spacing from parameters
    # Singh et al 1998.
    d_azm1 = d0[i] * (
        1
        + Q[i]
        * pos_neg
        * (1 - SampleGeometry[i] * (np.cos(np.deg2rad(azm - orientation[i]))) ** 2)
    )
    # d_azm1 = d0[i] * (1+ Q[i]*(1-SampleGeometry[i]*(np.cos(np.deg2rad(azm-orientation[i])))**2))

    # convert to two theta diffraction angle
    # n lambda = d2 sin (theta)
    tth_azm1 = 2 * np.arcsin(wavelength / 2 / d_azm1)  # radians
    r = distance * np.arctan(tth_azm1)  # mm

    # convert to x, y, on detector assuming perfectly flat and normal
    x = r * np.cos(np.radians(azm))  # mm
    y = r * np.sin(np.radians(azm))  # mm

    # move diffraction ring to offset position
    x = x - x0[i]  # mm
    y = y - y0[i]  # mm

    # convert to new parameters that are effectivly detector positions for perfect diffraction ring.
    r_new = np.sqrt(x**2 + y**2)  # mm
    tth_new = np.rad2deg(np.tan(r_new / distance))  # degrees

    # convert back to d-spacing
    d_new = wavelength / 2 / np.sin(np.deg2rad(tth_new / 2))

    # account for wrapping of tan(x) when reconstruct new angles.
    azm_new = np.rad2deg(np.unwrap(2 * (np.arctan(y / x) + np.pi / 2)) / 2 - np.pi / 2)

    # initate a set of lmfit parameters using the lowest level functions in CPF.
    master_params = Parameters()
    master_params = lmm.initiate_params(
        master_params,
        param_str="peak_" + str(0),
        comp="d",
        coeff_type="fourier",
        limits=[np.min(d_new), np.max(d_new)],
        value=d_new,
        trig_orders=trig_orders,
    )

    # fit azimuths with a Fourier series.
    fout = lmm.coefficient_fit(
        azimuth=azm_new,
        ydata=d_new,
        inp_param=master_params,
        param_str="peak_" + str(0) + "_d",
        fit_method="leastsq",
    )

    ## plot the data and fit just to make sure.
    if 0:
        # expand the data.
        fit_expanded = sf.coefficient_expand(azm_new, fout.params, comp_str="d")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        ax.plot(x, y, color="red")
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(str(Q[i]))

        fig, ax = plt.subplots(1, 2)
        # ax.plot(azm, r, 'r')
        ax[0].plot(d_new, azm_new, ".", color="red", label="obs")
        ax[0].plot(fit_expanded, azm_new, color="blue", label="fit")
        ax[0].legend()
        ax[1].plot(d_new - fit_expanded, azm_new, ".", color="red", label="diffs")
        ax[1].legend()

    # Write master_params to new_params dict object
    new_params = lmm.params_to_new_params(fout.params)

    dmax = np.max(d_azm1)
    dmin = np.min(d_azm1)

    if SampleDeformation[i] == "compression":
        dmin_ang = azm[d_azm1.argmin(axis=0)]
    else:
        dmin_ang = azm[d_azm1.argmax(axis=0)]

    # get converted aprameters.
    Singh_coef = cfc.fourier_to_crystallographic(
        new_params,
        SampleGeometry=f"{SampleGeometry[i]:1d}d",
        SampleDeformation=SampleDeformation[i],
        debug=debug,
    )

    if 0:
        print(
            f'Centroid:    Initial value {d0[i]:.3f}; recovered value {Singh_coef["dp"]:.3f}; difference {(d0[i]-Singh_coef["dp"])/d0[i]*100:4.4g}%'
        )
        print(
            f'Q:           Initial value {Q[i]:.3f}; recovered value {Singh_coef["differential"]:.3f}; difference {(Q[0]-Singh_coef["differential"])/Q[i]*100:4.4g}%'
        )
        print(
            f'orientation: Initial value {orientation[i]:.3f}; recovered value {Singh_coef["orientation"]:.3f}; difference {(orientation[i]-Singh_coef["orientation"])/360*100:4.4g}% of circle'
        )
        print(
            f'x:           Initial value {x0[i]:.4g}; recovered value {Singh_coef["x0"]:.3g}; difference {(x0[i]-Singh_coef["x0"])/x0[i]*100:4.4g}%'
        )
        print(
            f'y:           Initial value {y0[i]:.4g}; recovered value {Singh_coef["y0"]:.3g}; difference {(y0[i]-Singh_coef["y0"])/y0[i]*100:4.4g}%'
        )

    # store the parameters in arrays
    all_params.append(new_params["peak"][0]["d-space"])
    angles.append([orientation[i], Singh_coef["orientation"]])
    Qs.append([Q[i], Singh_coef["Q"]])
    d0s.append([d0[i], Singh_coef["dp"]])
    dmaxs.append([dmax, Singh_coef["d_max"]])
    dmins.append([dmin, Singh_coef["d_min"]])
    xs.append([x0[i], Singh_coef["x0"]])
    ys.append([y0[i], Singh_coef["y0"]])
    second_order_corefs.append(
        [fout.params.valuesdict()["peak_0_d3"], fout.params.valuesdict()["peak_0_d4"]]
    )
    chsqs.append(fout.chisqr)


# %% plot the outputs
# this makes alot of plots to illustrate what is going on with the fits.
if 1:
    angles = np.array(angles)
    d0s = np.array(d0s)
    dmaxs = np.array(dmaxs)
    dmins = np.array(dmins)
    Qs = np.array(Qs)
    xs = np.array(xs)
    ys = np.array(ys)
    second_order_corefs = np.array(second_order_corefs)
    chsqs = np.array(chsqs)
    all_params = np.array(all_params)

    # discard rubbish fits.
    keep = (chsqs) < 1
    angles = angles[keep, :]
    d0s = d0s[keep, :]
    Qs = Qs[keep, :]
    dmaxs = dmaxs[keep, :]
    dmins = dmins[keep, :]
    xs = xs[keep, :]
    ys = ys[keep, :]
    chsqs = chsqs[keep]
    all_params = all_params[keep, :]

    tths = 2 * np.rad2deg(np.arcsin(wavelength / 2 / d0s))  # degrees

    # plot the fitted and known angles.
    # colour against Q which is the elipticity. Q=0 no constraint on angles.
    fig, ax = plt.subplots()
    SC = ax.scatter(angles[:, 0], angles[:, 1], s=2, c=Qs[:, 0])
    ax.plot(
        [np.min(angles[:, 0]), np.max(angles[:, 0])],
        [np.min(angles[:, 0]), np.max(angles[:, 0])],
        ":",
        color="black",
    )
    ax.set_title("orientation")
    ax.set_xlabel("Known orientation")
    ax.set_ylabel("Fitted orientation")
    fig.colorbar(SC, ax=ax, label="Q")

    # plot differences between known and fit.
    angles2 = angles
    angles2[angles[:, 1] - angles[:, 0] > 90, 1] = (
        angles2[angles[:, 1] - angles[:, 0] > 90, 1] - 180
    )
    angles2[angles[:, 1] - angles[:, 0] < -90, 1] = (
        angles2[angles[:, 1] - angles[:, 0] < -90, 1] + 180
    )
    fig, ax = plt.subplots()
    ax.scatter(angles[:, 0], (angles2[:, 1] - angles2[:, 0]), s=2, c=Qs[:, 0])
    ax.set_title("Delta orientation angle")
    ax.set_xlabel("Known orientation (degrees)")
    ax.set_ylabel("Difference between fit and known orientations (degrees)")
    fig.colorbar(SC, ax=ax, label="Q")

    # plot log abs difference
    fig, ax = plt.subplots()
    ax.scatter(angles[:, 0], np.abs(angles2[:, 1] - angles2[:, 0]), s=2, c=Qs[:, 0])  #
    ax.set_title("Delta orientation angle")
    ax.set_xlabel("Known orientation (degrees)")
    ax.set_yscale("log")
    ax.set_ylabel("log absolute difference between fit and known orientations")
    fig.colorbar(SC, ax=ax, label="Q")

    # plot log abs difference vs. length of x0y0
    # shows and increase in misorientation with increasesing x0y0. -- the further the data is from the calubratio the worse the recovered values.
    fig, ax = plt.subplots()
    ax.scatter(
        np.sqrt(xs[:, 0] ** 2 + ys[:, 0] ** 2),
        np.abs(angles2[:, 1] - angles2[:, 0]),
        s=2,
        c=Qs[:, 0],
    )
    ax.set_title("Delta orientation angle")
    ax.set_xlabel("sqrt(x0^2 + y0^2)")
    ax.set_yscale("log")
    ax.set_ylabel("Difference between fit and known orientations (degrees)")
    fig.colorbar(SC, ax=ax, label="Q")

    # plot d0 known vs fitted
    fig, ax = plt.subplots()
    # ax.plot(azm, r, 'r')
    SC = ax.scatter(d0s[:, 0], d0s[:, 1], s=2, c=np.sqrt(xs[:, 0] ** 2 + ys[:, 0] ** 2))
    ax.plot(
        [np.min(d0s[:, 0]), np.max(d0s[:, 0])],
        [np.min(d0s[:, 0]), np.max(d0s[:, 0])],
        ":",
        color="black",
    )
    ax.set_title("d0")
    ax.set_xlabel("Known d0")
    ax.set_ylabel("Fitted d0")
    fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)")

    # plot differneces
    fig, ax = plt.subplots()
    ax.scatter(
        d0s[:, 0],
        (d0s[:, 1] - d0s[:, 0]) / d0s[:, 0] * 100,
        s=2,
        c=np.sqrt(xs[:, 0] ** 2 + ys[:, 0] ** 2),
    )
    ax.set_title("Delta d0")
    ax.set_xlabel("Known d0")
    ax.set_ylabel("% Difference between fit and known d0")
    fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)")

    # plot differences vs. two theta angle.
    # the smaller the angle the greater the difference.
    # the larger the x0y0 offset the larger the difference.
    fig, ax = plt.subplots()
    ax.scatter(
        tths[:, 0],
        np.abs(d0s[:, 1] - d0s[:, 0]),
        s=2,
        c=np.sqrt(xs[:, 0] ** 2 + ys[:, 0] ** 2),
    )
    ax.set_title("Delta d0")
    ax.set_xlabel("two theta (degrees)")
    ax.set_yscale("log")
    ax.set_ylabel("absolute difference between fit and known d0")
    fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)")

    # plot fitted Q falues
    fig, ax = plt.subplots()
    # ax.plot(azm, r, 'r')
    if 1:
        SC = ax.scatter(
            Qs[:, 0], Qs[:, 1], s=2, c=np.sqrt(xs[:, 0] ** 2 + ys[:, 0] ** 2)
        )
        fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)")
    else:
        SC = ax.scatter(Qs[:, 0], Qs[:, 1], s=2, c=d0s[:, 0])
        fig.colorbar(SC, ax=ax, label="d0")
    ax.plot([0, np.max(Qs[:, 0])], [0, np.max(Qs[:, 0])], ":", color="black")
    ax.set_title("Q")
    ax.set_xlabel("Known Q")
    ax.set_ylabel("Fitted Q")

    # plot difference in Q values.
    # shows that if you have a large centre offset relative to the strain you get a large diference in values.
    fig, ax = plt.subplots()
    SC = ax.scatter(
        Qs[:, 0],
        (Qs[:, 1] - Qs[:, 0]),
        s=2,
        c=np.sqrt(xs[:, 0] ** 2 + ys[:, 0] ** 2) / tths[:, 0],
    )
    # plot liness ot equal magnitude.
    ax.plot(
        [0, np.min([np.max(Qs[:, 0]), np.max(Qs[:, 1] - Qs[:, 0])])],
        [0, np.min([np.max(Qs[:, 0]), np.max(Qs[:, 1] - Qs[:, 0])])],
        ":",
        color="black",
    )
    ax.plot(
        [0, np.min([np.max(Qs[:, 0]), np.max(Qs[:, 1] - Qs[:, 0])])],
        [0, -np.min([np.max(Qs[:, 0]), np.max(Qs[:, 1] - Qs[:, 0])])],
        ":",
        color="black",
    )
    ax.set_title("Delta Q")
    ax.set_xlabel("Known Q")
    ax.set_ylabel("Difference between fit and known Q")
    fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)/tth")

    # plot max and mins
    fig, ax = plt.subplots()
    ax.plot([0, np.max(dmaxs[:, 0])], [0, np.max(dmaxs[:, 0])], ":", color="black")
    SC = ax.scatter(dmaxs[:, 0], dmaxs[:, 1], s=2, c=Qs[:, 0])
    fig.colorbar(SC, ax=ax, label="Q")
    ax.set_title("d max")
    ax.set_xlabel("Known d_max (approximated)")
    ax.set_ylabel("Fitted d_max")

    fig, ax = plt.subplots()
    ax.plot([0, np.max(dmins[:, 0])], [0, np.max(dmins[:, 0])], ":", color="black")
    SC = ax.scatter(dmins[:, 0], dmins[:, 1], s=2, c=Qs[:, 0])
    fig.colorbar(SC, ax=ax, label="Q")
    ax.set_title("d min")
    ax.set_xlabel("Known d_min (approximated)")
    ax.set_ylabel("Fitted d_min")

    # %% x0 and y0
    # the relationship between the Fourier coefficients and the offset of the centre of the DEbye-Scherer ring is not obvious (at least to me).
    # these plots are the beginings of the plots to show the conversion.
    if debug:
        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(xs[:, 0], xs[:, 1], s=2, c=Qs[:, 0])
        ax.set_ylabel("Fitted x0s")
        fig.colorbar(SC, ax=ax, label="Q")
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("x0")
        ax.set_xlabel("Known x0s")
        # ax.set_yscale('log')

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(ys[:, 0], ys[:, 1], s=2, c=Qs[:, 0])
        ax.set_ylabel("Fitted y0s")
        fig.colorbar(SC, ax=ax, label="Q")
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("y0")
        ax.set_xlabel("Known y0s")
        # ax.set_yscale('log')

        # x conversion
        # works for q=0
        x_conv = (xs[:, 1]) / all_params[:, 0] ** 2 * distance * wavelength
        # finite Q
        semiQ = np.sqrt(all_params[:, 3] ** 2 + all_params[:, 4] ** 2)
        # x_conv = (xs[:,1])/all_params[:,0]**2*distance*wavelength *(1+xs[:,1]/np.tan(all_params[:,4]/all_params[:,3]))
        # x_conv = (xs[:,1])/d0s[:,1]**2*distance*wavelength + (Qs[:,1])*np.absolute(xs[:,0])
        # x_conv = (xs[:,1])/d0s[:,1]**2*distance*wavelength + 3/2*(Qs[:,1])*(xs[:,0])# best so far angles = 0, ys = 0
        # x_conv = (xs[:,1])/d0s[:,1]**2*distance*wavelength works for q=0

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(xs[:, 0], x_conv, s=2, c=Qs[:, 0])
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("x0 -- 1st  plot")
        ax.set_xlabel("Known x0s")
        ax.set_ylabel("Fitted foutier term/d0s[:,1]**2*distance*wavelength")
        fig.colorbar(SC, ax=ax, label="Q")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(xs[:, 0], x_conv - xs[:, 0], s=2, c=Qs[:, 1])
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("x0")
        ax.set_xlabel("Known x0s")
        # ax.set_yscale('log')
        ax.set_ylabel("(Fitted-known)/known x0s")
        fig.colorbar(SC, ax=ax, label="Qs")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(xs[:, 0], x_conv, s=2, c=d0s[:, 0])
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("x0")
        ax.set_xlabel("Known x0s")
        # ax.set_yscale('log')
        ax.set_ylabel("(Fitted-known)/known x0s")
        fig.colorbar(SC, ax=ax, label="d0")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(xs[:, 0], (x_conv) - xs[:, 0], s=2, c=Qs[:, 0])
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("x0")
        ax.set_xlabel("Known x0s")
        # ax.set_yscale('log')
        ax.set_ylabel("(Fitted-known)/known x0s")
        fig.colorbar(SC, ax=ax, label="Q")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(xs[:, 0], (x_conv) - xs[:, 0], s=2, c=angles[:, 0])
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("x0")
        ax.set_xlabel("Known x0s")
        # ax.set_yscale('log')
        ax.set_ylabel("(Fitted-known)/known x0s")
        fig.colorbar(SC, ax=ax, label="angle")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(xs[:, 0], (x_conv) - xs[:, 0], s=2, c=d0s[:, 0])
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("x0")
        ax.set_xlabel("Known x0s")
        # ax.set_yscale('log')
        ax.set_ylabel("(Fitted-known)/known x0s")
        fig.colorbar(SC, ax=ax, label="d0")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(
            ys[:, 0], ys[:, 1] / d0s[:, 1] ** 2 * distance * wavelength, s=2, c=Qs[:, 0]
        )
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("y0")
        ax.set_xlabel("Known y0s")
        ax.set_ylabel("Fitted foutier term/d0s[:,1]**2*distance*wavelength")
        fig.colorbar(SC, ax=ax, label="Q")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(ys[:, 0], (ys[:, 1] / d0s[:, 1] ** 2), s=2, c=Qs[:, 0])
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("y0")
        ax.set_xlabel("Known y0s")
        ax.set_yscale("log")
        ax.set_ylabel("(Fitted-known)/known y0s")
        fig.colorbar(SC, ax=ax, label="Q")

        fig, ax = plt.subplots()
        # ax.plot(azm, r, 'r')
        SC = ax.scatter(
            ys[:, 0],
            (ys[:, 1] / d0s[:, 1] ** 2 * distance * wavelength - ys[:, 0]) / ys[:, 0],
            s=2,
            c=Qs[:, 0],
        )
        # ax.plot(xs[:,0], xs[:,1], '.', color='red')
        # ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
        ax.set_title("y0")
        ax.set_xlabel("Known y0s")
        ax.set_yscale("log")
        ax.set_ylabel("(Fitted-known)/known y0s")
        fig.colorbar(SC, ax=ax, label="Q")

    """
    # plot the length of the x,y vectors    
    fig, ax = plt.subplots()
    #ax.plot(azm, r, 'r')
    #SC = ax.scatter(ys[:,0], ys[:,1], s=2, c=d0s[:,0])
    SC = ax.scatter(np.sqrt(ys[:,0]**2+xs[:,0]**2), np.sqrt(ys[:,1]**2+xs[:,1]**2), s=2, c=d0s[:,0])
    #SC = ax.scatter(np.sqrt(ys[:,0]**2+ys[:,0]**2)*d0s[:,0], np.sqrt(ys[:,1]**2+ys[:,1]**2)*d0s[:,1], s=2, c=Qs[:,0])
    #ax.plot(ys[:,0], ys[:,1], '.', color='red')
    #ax.plot([np.min(ys[:,0]),np.max(ys[:,0])], [np.min(ys[:,0]),np.max(ys[:,0])], ':', color='black')
    ax.set_title("centroid postion")
    ax.set_xlabel("Known sqrt(x0^2 + y0^2)")
    ax.set_ylabel("Fitted sqrt(x0^2 + y0^2)")
    fig.colorbar(SC, ax=ax, label="d0s")
    
    # plot the length of the x,y vectors    
    fig, ax = plt.subplots()
    #ax.plot(azm, r, 'r')
    #SC = ax.scatter(ys[:,0], ys[:,1], s=2, c=d0s[:,0])
    SC = ax.scatter(np.sqrt(ys[:,0]**2+xs[:,0]**2)/d0s[:,0], np.sqrt(ys[:,1]**2+xs[:,1]**2)/d0s[:,1], s=2, c=Qs[:,0])
    #ax.plot(ys[:,0], ys[:,1], '.', color='red')
    #ax.plot([np.min(ys[:,0]),np.max(ys[:,0])], [np.min(ys[:,0]),np.max(ys[:,0])], ':', color='black')
    ax.set_title("centroid postion")
    ax.set_xlabel("Known sqrt(x0^2 + y0^2)")
    ax.set_ylabel("Fitted sqrt(x0^2 + y0^2)")
    fig.colorbar(SC, ax=ax, label="d0s")
    """
    """
    fig, ax = plt.subplots()
    #ax.plot(azm, r, 'r')
    SC = ax.scatter(xs[:,0], xs[:,1]/d0s[:,1]**2, s=2, c=np.sqrt(xs[:,0]**2+ys[:,0]**2)*d0s[:,0])
    #ax.plot(xs[:,0], xs[:,1], '.', color='red')
    #ax.plot([np.min(xs[:,0]),np.max(xs[:,0])], [np.min(xs[:,0]),np.max(xs[:,0])], ':', color='black')
    ax.set_title("xs")
    ax.set_xlabel("Known xs")
    ax.set_ylabel("Fitted xs")
    fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)")
    
    
    fig, ax = plt.subplots()
    #ax.plot(azm, r, 'r')
    SC = ax.scatter(np.sqrt(xs[:,0]**2+ys[:,0]**2)/d0s[:,0], xs[:,1]/d0s[:,1]**2, s=2, c=np.sqrt(xs[:,0]**2+ys[:,0]**2))# d0s[:,0])
    #ax.plot(angles[:,0], d0s[:,0]-d0s[:,1], '.', color='red')
    #ax.set_title("Delta d0 (known-recovered)")
    #ax.set_xlabel("strain orientation (degrees)")
    ax.set_xlabel("Known lenght(x0,y0)/d0")
    ax.set_ylabel("x0/d0")
    fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)")
    
    
    
    
    
    fig, ax = plt.subplots()
    #ax.plot(azm, r, 'r')
    SC = ax.scatter(ys[:,0], ys[:,1], s=2, c=d0s[:,0])
    #SC = ax.scatter(ys[:,0], ys[:,1]/d0s[:,1]**2*distance/6, s=2, c=d0s[:,0])
    #ax.plot(ys[:,0], ys[:,1], '.', color='red')
    #ax.plot([np.min(ys[:,0]),np.max(ys[:,0])], [np.min(ys[:,0]),np.max(ys[:,0])], ':', color='black')
    ax.set_title("ys")
    ax.set_xlabel("Known ys")
    ax.set_ylabel("Fitted ys")
    fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)")    
        
    fig, ax = plt.subplots()
    #ax.plot(azm, r, 'r')
    #SC = ax.scatter(ys[:,0], ys[:,1], s=2, c=tths[:,0])
    SC = ax.scatter(ys[:,0], ys[:,1]/d0s[:,1], s=2, c=tths[:,0])
    #ax.plot(ys[:,0], ys[:,1], '.', color='red')
    #ax.plot([np.min(ys[:,0]),np.max(ys[:,0])], [np.min(ys[:,0]),np.max(ys[:,0])], ':', color='black')
    ax.set_title("ys")
    ax.set_xlabel("Known ys")
    ax.set_ylabel("Fitted ys")
    fig.colorbar(SC, ax=ax, label="tths")
    
    fig, ax = plt.subplots()
    #ax.plot(azm, r, 'r')
    #SC = ax.scatter(ys[:,0], ys[:,1], s=2, c=d0s[:,0])
    SC = ax.scatter(np.sqrt(ys[:,0]**2+ys[:,0]**2)/d0s[:,0], np.sqrt(xs[:,1]**2+s[:,1]**2)/d0s[:,1], s=2, c=d0s[:,0])
    #ax.plot(ys[:,0], ys[:,1], '.', color='red')
    #ax.plot([np.min(ys[:,0]),np.max(ys[:,0])], [np.min(ys[:,0]),np.max(ys[:,0])], ':', color='black')
    ax.set_title("ys")
    ax.set_xlabel("Known ys")
    ax.set_ylabel("Fitted ys")
    fig.colorbar(SC, ax=ax, label="sqrt(x0^2 + y0^2)")
    
    
    """
