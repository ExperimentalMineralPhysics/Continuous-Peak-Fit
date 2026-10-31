# # Proofs and maths for cpf fits to crystallographic strains.

# Models for converting between Fourier series and differentail strain formultions. 
# 
# For continuous peak fit. 


import cpf
import cpf.output_formatters.convert_fit_to_crystallographic as cftc
from cpf.lmfit_model import coefficient_fit, initiate_params, gather_param_errs_to_list
import cpf.series_functions as sf
from lmfit import Parameters
from cpf.series_functions import series_properties
from cpf.output_formatters.convert_fit_to_crystallographic import (
    fourier_to_crystallographic,
)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
# from scipy import optimize
from tqdm import trange
from scipy.optimize import fmin, minimize, leastsq



# ## Define strain and stress equations

# modified after Singh et al., 1998.
# 
# 
# The stress state in a generic sample is:
# $$\sigma_{ij} = \begin{bmatrix} \sigma_{11} & 0  & 0 \\  0 & \sigma_{22} & 0 \\    0  & 0  & \sigma_{33} \end{bmatrix} $$
# this can be split into pressure and differential stress ($t=\sigma_{11}-\sigma_{33}$) terms:
# $$\sigma_{ij} = \begin{bmatrix} \sigma_{P} & 0  & 0 \\  0 & \sigma_{P} & 0 \\    0  & 0  & \sigma_{P} \end{bmatrix} + 
# \begin{bmatrix} 2t/3 & 0  & 0 \\  0 & -t/3 & 0 \\    0  & 0  & -t/3 \end{bmatrix} $$
# 
# 
# The differentail elastic strain produced from the differential stress with respect to the stress state is (Merkel and Hilaret):
# 
# $$\varepsilon_D(hkl, \Psi) = \frac{d_m(hkl, \Psi) - d_P(hkl)}{d_P(hkl)} = Q(hkl)(1-3cos^2\Psi) $$
# 
# where $d_P(hkl)$ is the d-spacing from pressure alone -- *is the volumetric mean d-spacing of the (hkl) at the  conditions (my interpretation)*. And $Q(hkl)$ is the lattice strain parameter which equals (Singh et al., 1998):
# 
# $$ Q(hkl) = \frac{t}{3}G$$
# 
# <!-- Under a differential stress variation in lattice strains ($d_{hkl}$) as a function of azimuth are (Singh et al., 1998):
# 
# $$d_\chi(h k l) = d_P(h k l)[1+(1-3\cos^2\Psi)Q(hkl)]$$
# 
# and 
# 
# $$ Q(hkl) = (t/3)G$$
# 
# where $\Psi$ is the azimuthal angle from the maximum stress direction and $d = n\lambda/2\sin(2\theta)$. $G$ is the shear modulus, withing he Reuss-Voight bounds. -->
# 
# This is defined for a three dimensional system. In a 2 dimensional system the scale factor of 3 should be two and so the equation generalsises to:
# 
# $$\varepsilon_D(hkl, \Psi) = Q(hkl)(1-Rcos^2\Psi) $$
# 
# $$ Q(hkl) = \frac{t}{R}G$$
# 
# where $R$ is the (real) dimensions of the system. For any diffraction system $R = 2$ or 3. 
# 
# 
# The measured values $d_m(hkl, \Psi)$ in in $R$ dimensional system are therefore: 
# $$d_m(hkl, \Psi) = d_P(hkl) + d_P(hkl) Q(hkl)(1-Rcos^2\Psi) $$
# drop $\Psi$ and (hkl)...
# $$d_m = d_P + d_P Q(1-Rcos^2\Psi) $$

# ### apply Sin identiies
# 
# Appling the sin identity:  
# $$\cos(2\theta) = 2cos^2\theta -1 $$ 
# 
# to 
# $$d_m = d_P + d_P Q(1-Rcos^2\Psi) $$
# 
# gives:
# \begin{eqnarray}
# d_m  = & d_P + d_P Q(1-Rcos^2\Psi) \\
#  = & d_P + d_P Q(1-R( \frac{\cos(2\Psi) +1 }{2} ) ) \\
#  = & d_P + d_P Q - \frac{d_P Q R}{2} + \frac{d_P Q R}{2}\cos(2\Psi)  \\
# \end{eqnarray}
# 

# However, the stress field is likely to be misaligned to the detector. Therefore:
# \begin{equation}
# \Psi = \chi + \psi
# \end{equation}
# Expanding the sin identity for $\cos^2(\chi + \psi)$ gives:
# \begin{equation}
# cos^2(\chi + \psi) = \frac{1}{2} + \frac{1}{2}\cos(2\psi)\cos(2\chi) - \frac{1}{2}\sin(2\psi)\sin(2\chi)
# \end{equation}
# Substituting into equation for $d_m$ and rearranging gives:
# \begin{equation}
# d_m  = d_P + d_P Q - \frac{d_P Q R}{2} - \frac{d_P Q R}{2}\cos(2\psi)\cos(2\chi) + \frac{d_P Q R}{2}\sin(2\psi)\sin(2\chi)
# \end{equation}

# ### Fourier series
# Assuming in the equation above the terms $\cos(2\chi)$ and $\sin(2\chi)$ are terms of a Fourier series. 
# And using the sin indentity:
# \begin{equation}
# A\cos(n\omega) + B\sin(n\omega) = \sqrt{(A^2+B^2)} \cos[(n\omega) - \tan^{-1}(B/A)]
# \end{equation}
# 
# ### $A, B$
# combining terms and substitiuting
# \begin{eqnarray}
# \sqrt{(A^2+B^2)}  =& \sqrt{(-\frac{d_P Q R}{2}\cos(2\psi))^2+(\frac{d_P Q R}{2}\sin(2\psi))^2} \\
# =& \sqrt{\frac{d_P^2 Q^2 R^2}{4} (\cos(2\psi)^2+\sin(2\psi)^2)} \\
# =& \frac{d_P Q R}{2} \\
# \end{eqnarray}
# 
# 
# ### $\tan^{-1}(B/A)$
# \begin{eqnarray}
# \tan^{-1}(B/A)  =& \tan^{-1}\frac{\frac{d_P Q R}{2}\sin(2\psi)}{-\frac{d_P Q R}{2}\cos(2\psi)} \\
# =& \tan^{-1}\tan(-2\psi) \\
# =& -2\psi
# \end{eqnarray}
# 
# 
# ### resubstituting
# restubstituting back into the equation for $d_m$:
# \begin{eqnarray}
# d_m  =& d_P + d_P Q - \frac{d_P Q R}{2} - \frac{d_P Q R}{2}\cos(2\psi)\cos(2\chi) + \frac{d_P Q R}{2}\sin(2\psi)\sin(2\chi)\\
# =& d_P + d_P Q - \frac{d_P Q R}{2} +\frac{d_P Q R}{2} \cos[2\chi - -2\psi] \\
# =& d_P + d_P Q - \frac{d_P Q R}{2} +\frac{d_P Q R}{2} \cos[2 (\chi + \psi) ] \\
# =& d_P + d_P Q - \frac{d_P Q R}{2} +\frac{d_P Q R}{2} \cos[2\Psi] 
# \end{eqnarray}
# 
# which gives the equation we thought of earlier. Therefore can use Fourier series with these coefficients. 
# 

# ## Confirm implementations in code work.
# 

# and limitations thereof. 


def Differential_strain(d0, Q, orientation, azimuths, dimensionality=3, compression="CompressionPositive"):
    """
    retrun d-spacings around azimuth given d0, Q and orientation

    after Singh et al (1998) and others.

    base equation that we should be reproducting with the fourier series.
    """
    if azimuths is None:
        azimuths = np.linspace(0, 360, num=361)
    if compression=="CompressionPositive":
        out = d0 * (1 - (1 - dimensionality*(np.cos(np.deg2rad(azimuths-orientation))**2))* Q)
    else:
        out = d0 * (1 +  (1 - dimensionality*(np.cos(np.deg2rad(azimuths-orientation))**2))* Q)
    return out



def off_centre_strain(d0, Q, orientation, x0=0,y0=0, wavelength=0.1, det_dist=1000, azimuths=None, geometry="3d", compression="CompressionPositive"):
    """
    Given a generic set of parameters make diffraction ring positions.

    calculates d-spacing using d0, Q and orientation, 
    Uses this to get two-theta using wavelengths,
    converts to x, y, (pixels) using detector distance
    and offset for x,y
    the transforms back into an apparent two theta from x,y, positions in detector. 

    Parameters
    ----------
    d0 : TYPE
        DESCRIPTION.
    Q : TYPE
        DESCRIPTION.
    orientation : TYPE
        DESCRIPTION.
    x0 : TYPE
        DESCRIPTION.
    y0 : TYPE
        DESCRIPTION.
    azimuths : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    if azimuths is None:
        azimuths = np.linspace(0, 360, num=361)
    if cftc._sample_geometries(as_list=False)[geometry.lower()] == "2d":
        dims = 2
    else:
        dims=3
    dspace = Differential_strain(d0, Q, orientation, azimuths, dimensionality=dims, compression=compression)
    tth = np.rad2deg(2*np.arcsin(wavelength/2/dspace))
    # print(tth)

    # x_centred = dspace * np.cos(np.deg2rad(azimuths))
    # y_centred = dspace * np.sin(np.deg2rad(azimuths))

    x_centred = det_dist * np.tan(np.deg2rad(tth)) * np.cos(np.deg2rad(azimuths))
    y_centred = det_dist * np.tan(np.deg2rad(tth)) * np.sin(np.deg2rad(azimuths))
    x_offcentre = x_centred + x0
    y_offcentre = y_centred + y0
    # d_app = np.sqrt(x_offcentre**2 + y_offcentre**2)
    r = np.sqrt(x_offcentre**2 + y_offcentre**2)
    tth_app = np.rad2deg(np.arctan(r/det_dist))
    # plt.plot(tth, tth_app)
    azimuth_app = np.rad2deg(np.arctan(y_offcentre/x_offcentre))
    # unwrap new azimuths
    for i in range(len(azimuth_app)):
        if i==0:
            pass
        else:
            if azimuth_app[i-1] > azimuth_app[i]:
                azimuth_app[i:] += 180
    azimuth_app = azimuth_app % 360

    return tth_app, azimuth_app, x_offcentre, y_offcentre,



def fourier_series(coeffs, azimuths):
    """
    Calculate fourier series at azimuths given coeffs.
    """
    if azimuths is None:
        azimuths = np.linspace(0, 360, num=361)    
    # coeffs = np.abs(coeffs)
    # print(coeffs)
    out = np.ones(azimuths.shape) * coeffs[0]
    orders = (len(coeffs)-1) /2
    for i in range(int(orders)):
        out += coeffs[2*i+1] * np.sin(np.deg2rad(azimuths*[i+1]))
        out += coeffs[2*i+2] * np.cos(np.deg2rad(azimuths*[i+1]))

    return out


def fourier_series_tth(coeffs, azimuths, wavelength, log=False):

    if azimuths is None:
        azimuths = np.linspace(0, 360, num=361)

    if log == True:
        coeffs = 10**coeffs

    # print(coeffs)
    dspace = fourier_series(coeffs, azimuths)
    # print(dspace.shape)
    # plt.plot(dspace)

    # print(np.min(dspace))
    out_tth = np.rad2deg(2*np.arcsin(wavelength/2/dspace))
    # print("\n\n", wavelength/2/dspace, out_tth)
    # plt.plot(out_tth)
    # stop
    return out_tth


# %% set up major problem parameters. Related to the detector
wavelength = 0.2 # angstroms
det_dist = 500 #mm    


# ### Show Fourier series works for single case
# show Fourier series can fit a single diffraction centroid. 
if 0:
    # set parameters. 
    num_samples = 1
    if num_samples == 1:
        d0 = 2.1
        deltad = .00
        orientation = 0
        x0 = 18
        y0 = 18
    else:
        d0 = np.random.uniform(low=.5, high=3, size= num_samples)
        deltad = np.random.uniform(low=0.00, high=0.05, size= num_samples) 
        orientation = np.random.uniform(low=1, high=180, size= num_samples) 
        x0 = np.random.uniform(low=0.0, high=5, size= num_samples)
        y0 = np.random.uniform(low=0.0, high=5, size= num_samples) 
    
    # error function to minimize
    log = False
    e = lambda coef, d, azim: ((fourier_series_tth(coef,azim,wavelength,log=log)-d)**2).sum()
    
    # convert crystallographic values in to apparent d-spacings in detector.
    d, o, x,y = off_centre_strain(d0, deltad, orientation, x0, y0, wavelength, det_dist)
    
    
    azim = o#np.linspace(0, 360, num=361)
    # best fit paramets with fourier series.
    
    # fitting the data with fmin
    order = 2
    coeffs = np.ones(1+2*order)*.01
    coeffs[0] = 2
    if log == True:
        # print(np.log10(coeffs)) 
        p = fmin(e, np.log10(coeffs), args=(d,o),xtol=1e-06, ftol=1e-6)
        print('estimater parameters: ', 10**p)
    else:
        p = fmin(e, coeffs, args=(d,o),xtol=1e-06, ftol=1e-6, maxiter=100000, maxfun=100000,)
        print('estimater parameters: ', p)

    fig = plt.figure()
    spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[3, 1])
    ax0 = fig.add_subplot(spec[0])
    ax0.plot(d,o, 'o', fourier_series_tth(p,azim,wavelength,log=log),azim, '.')
    ax0.set_xlabel("two theta (degrees)")
    ax0.set_ylabel("azimuth (degrees)")
    ax1 = fig.add_subplot(spec[1])
    ax1.scatter((d-fourier_series_tth(p,azim,wavelength,log=log))/d*100,o,s=0.1)
    ax1.set_xlabel("differences (%)")
    
    
    

# %% Show works for centred diffraction on all geometries.
# keep X0=0 and Y0=0. then recover input values. 

geometries = cftc._sample_geometries(as_list=False)
strains = ["CompressionPositive", "ExtensionPositive"]
num_samples = 100

for SampleGeometry in list(geometries.keys()):#    list(["extension"]):#
    
    for SampleStrain in strains:
                        
        if SampleGeometry == "2d" or SampleGeometry == "shear" or SampleGeometry == "axial":
            d0 = [0.2, 6]
            deltad = [0, .3]
            orientation = [-180,180]
            
        elif SampleGeometry == "compression" or SampleGeometry == "3d": # 3d is assumed to be a compression geometry
            d0 = [0.2, 6]
            deltad = [0, .3]
            orientation = [-180,180]
            
        elif SampleGeometry == "extension":
            d0 = [0.2, 6]
            deltad = [-0.3, 0]
            orientation = [-180,180]
        
        elif SampleGeometry == "vertical" or SampleGeometry == "n-s" :
            d0 = [0.2, 6]
            deltad = [-.3, .3]
            orientation = [45,135]
        
        elif SampleGeometry == "horizontal" or SampleGeometry == "e-w":
            d0 = [0.2, 6]
            deltad = [-.3, .3]
            orientation = [-45,45]
        
        elif SampleGeometry == "ne-sw":
            d0 = [0.2, 6]
            deltad = [-.3, .3]
            orientation = [0,90]
        
        elif SampleGeometry == "nw-se":
            d0 = [0.2, 6]
            deltad = [-.3, .3]
            orientation = [-90,0]
            
        if 0:
            X0 = [0, 5]
            Y0 = [0, 5]
            deltad=[0,0]
            orientation = [0,0]
        else:
            X0 = [0, 0]
            Y0 = [0, 0]
    
        # solve for multiple instances -- then plot them all. 
        d0s          = np.random.uniform(low=d0[0], high=d0[1], size= num_samples)
        deltads      = np.random.uniform(low=deltad[0], high=deltad[1], size= num_samples)  
        orientations = np.random.uniform(low=orientation[0], high=orientation[1], size= num_samples) 
        x0s          = np.random.uniform(low=X0[0], high=X0[1], size= num_samples) 
        y0s          = np.random.uniform(low=Y0[0], high=Y0[1], size= num_samples) 
        input_as_dict = {"d0s": d0s,
                         "deltads": deltads,
                         "orientations": orientations,
                         "x0s": x0s,
                         "y0s": y0s
                        }
        
        all_properties = []
        for i in trange(num_samples, desc='Models', leave=True):
        
            #calculate new properties
            tth, o, x,y = off_centre_strain(d0s[i], deltads[i], orientations[i], x0s[i], y0s[i], wavelength, det_dist, geometry=SampleGeometry, compression=SampleStrain)
            d = wavelength/2/np.sin(np.deg2rad(tth/2))
            azim = o
        
            #use cpf to do the work
            # best fit paramets with fourier series.
            master_params = Parameters()  
            master_params = initiate_params(
                master_params,
                "peak_0",
                "d",
                coeff_type=0,
                trig_orders=2,
                limits=None,
                value=None,
                types=True,
            )
            out = coefficient_fit(d, o, master_params, param_str="peak_0_d")
            # recover values from fourier coefficients using cpf code. 
            p, p_err = gather_param_errs_to_list(out.params, "peak_0", "d")
        
            pp_coefs = {}
            for pps in range(len(p)):
                if "p"+str(pps) not in pp_coefs:
                    pp_coefs["p"+str(pps)] = []
                pp_coefs["p"+str(pps)].append(p[pps])
        
        
            #forge cpf style data structure for coefficients
            # should probably do this using cpf.lmfit_model methods
            cfp_data_structure = [{"peak":[{"d-space":p, 
                                            "d-space_type":"fourier",
                                            "d-space_err":np.abs(np.array(p)/100), 
                                           }]}]
            crystallographic_values = fourier_to_crystallographic(
                                cfp_data_structure,
                                SampleGeometry=SampleGeometry,
                                SampleStrain=SampleStrain,
                                subpattern=0,
                                peak=0,
                            )
            all_properties.append(crystallographic_values | pp_coefs)
            # properties = series_properties(
            #     cfp_data_structure, subpattern=0, peak=0, param="d-space",
            #     start_end = [0,360]
            # )
            # all_properties.append(properties)
            
        output_as_dict = {}
        for i in all_properties[0].keys():
            if i not in output_as_dict:
                output_as_dict[i] = []
            for ap in all_properties:
                output_as_dict[i].append(ap[i])
        
    
        #plot correlations  produced.
        
        plots = {"d mean": {"given": "d0s",
                         "calc": "dp",
                        "colours": "p3"},
                "Q/delta d": {"given": "deltads",
                         "calc": "Q",
                        "colours": "p4"},
                 "orientation": {"given": "orientations",
                         "calc": "orientation",
                        "colours": "p4"},
                }
        
        print("If all plots are 1-1 input vs. output then correct. Orientations might be a bit funny because of cycling round")
        fig, ax = plt.subplots(1,len(plots))
        for pp in range(len(plots.keys())):
            p = list(plots.keys())[pp]
            ax[pp].set_title(p)
            ax[pp].plot([np.min(input_as_dict[plots[p]["given"]]),np.max(input_as_dict[plots[p]["given"]]) ], [np.min(input_as_dict[plots[p]["given"]]),np.max(input_as_dict[plots[p]["given"]]) ], 'r:')
            if "colours" in list(plots[p].keys()):
                if plots[p]["colours"] in list(input_as_dict):
                    clrs = input_as_dict[plots[p]["colours"]]    
                else:
                    clrs = output_as_dict[plots[p]["colours"]]
                plttd = ax[pp].scatter(input_as_dict[plots[p]["given"]], output_as_dict[plots[p]["calc"]], s=2, c = clrs, cmap="jet")
                cbar = plt.colorbar(plttd, shrink=0.6, location='bottom')
                cbar.set_label(plots[p]["colours"], rotation=270)
            else:
                plttd = ax[pp].scatter(input_as_dict[plots[p]["given"]], output_as_dict[plots[p]["calc"]], s=2)
            ax[pp].set_xlabel(f"input: {plots[p]['given']}")
            ax[pp].set_ylabel(f"output: {plots[p]["calc"]}")
        fig.suptitle(f"{SampleGeometry} -- {SampleStrain}")
            
# also need to check some values are valid:
# differnetial = dmax-dmin