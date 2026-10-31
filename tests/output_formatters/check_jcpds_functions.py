# # Proofs and maths for stresses extracted from jcpds calculations.

# For continuous peak fit. 


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %% check stress is correct from jcpdc functions. 

"""
Q is related to the differential stess ($\sigma_d$) and the shear modulus (G) (e.g. Singh et al 1998, Eq4; Li et al 2004, J. Appl. Phys, Eq.6)
The relatinship in both references above rearrange to give:
    \sigma_d = dim * Q * (2*G) 
where dim is the dimensionality of the system (dim= 2 or 3). 
Formally each term (apart from dim) is dependent on the hkl but for testing purposeses the hkl dependece can be ignored.
"""
from cpf.output_formatters import jcpds

def stress_to_strain(jcpds_obj, pressure, temperature, stress, orientation, dimensionality=3, compression="CompressionPositive"):
    """
    Calculate strains from pressure, stress and elastic moduli.
    Then fit and get the same number back. 
    
    """
    
    if compression=="CompressionPositive":
        scl = -1
    else:
        scl = 1
        
    
    g0t = jcpds_obj.g0 + jcpds_obj.dg0dt * (temperature - 298.)
    g0pt = jcpds_obj.g0p0 + jcpds_obj.dg0pdt * (temperature - 298.)
    G = g0t + g0pt * (pressure)

    for refl in jcpds_obj.reflections:
        refl.orientationobs = orientation
        
        #FIXME -- change this to take stress and calculate differential from it
        
        #after Singh et al 1998, Eq4
        Q = stress /dimensionality / (G*2)
        
        d = refl.d
        dmin = d*(1 + scl * (1-dimensionality * 1) *Q) #cos^2(0) = 1
        
        differential = (dmin - d) / (dimensionality-1) * dimensionality
        refl.differentialobs = differential
        
    return jcpds_obj


# test once
if 1:
    jcpds_file="../../Example2-MgO/MgO_ICSD.jcpds"
    
    jcpds_obj = jcpds.jcpds()
    jcpds_obj.load_jcpds(jcpds_file)
    
    
    
    pressure = 10
    temperature = 298
    stress = 1
    orientation = 0


    jcpds_obj.compute_v0()
    # these next two lines should give the same answers.
    jcpds_obj.compute_volume(pressure=pressure, temperature=temperature)
    jcpds_obj = stress_to_strain(jcpds_obj, pressure=pressure, temperature=temperature, stress=stress, orientation=orientation)    
    jcpds_obj.compute_stress_per_hkl()
    
    jcpds_obj.compute_d(pressure=pressure, temperature=temperature)
    jcpds_obj = stress_to_strain(jcpds_obj, pressure=pressure, temperature=temperature, stress=stress, orientation=orientation)    
    jcpds_obj.compute_stress_per_hkl()



# %% test many times

# geometries = cftc._sample_geometries(as_list=False)
# strains = ["CompressionPositive", "ExtensionPositive"]
num_samples = 500


pressure_range = [0,30]
temperature_range = [25,1800]
stress_range = [0,.2]
orientation_range = [0,180]

    
# solve for multiple instances -- then plot them all. 
pressures      = np.random.uniform(low=pressure_range[0], high=pressure_range[1], size= num_samples)
tempaeratures  = np.random.uniform(low=temperature_range[0], high=temperature_range[1], size= num_samples)  
stresses       = np.random.uniform(low=stress_range[0], high=stress_range[1], size= num_samples) 
orientations   = np.random.uniform(low=orientation_range[0], high=orientation_range[1], size= num_samples) 


input_as_dict = {"pressures": pressures,
                 "temperatures": tempaeratures,
                 "stresses": stresses,
                 "orientations": orientations,
                }
        

jcpds_file="../../Example2-MgO/MgO_ICSD.jcpds"
    
jcpds_obj = jcpds.jcpds()
jcpds_obj.load_jcpds(jcpds_file)

for j in range(len(jcpds_obj.reflections)):
    if j != 0:
        jcpds_obj.remove_reflection(-1)
jcpds_obj.compute_v0()

all_recovered = []
for i in range(num_samples):

    jcpds_obj.compute_volume(pressure=input_as_dict["pressures"][i], temperature=input_as_dict["temperatures"][i])
    jcpds_obj = stress_to_strain(jcpds_obj, 
                                 pressure=input_as_dict["pressures"][i], 
                                 temperature=input_as_dict["temperatures"][i], 
                                 stress=input_as_dict["stresses"][i], 
                                 orientation=input_as_dict["orientations"][i])    
    jcpds_obj.compute_stress_per_hkl()

    recovered = {"pressure": jcpds_obj.pressure.nominal_value,
                 "pressure_err": jcpds_obj.pressure.std_dev,
                 "temperature": jcpds_obj.temperature,
                 "stress":jcpds_obj.compute_stress_per_hkl()[0].nominal_value,
                 "stress_err":jcpds_obj.compute_stress_per_hkl()[0].std_dev}    
    
    all_recovered.append(recovered)
    
    
all_recovered = pd.DataFrame(all_recovered)
        
# plot stresses



        
plots = {"stresses": {"given": "stresses",
                 "calc": "stress",
                "colours": "pressures"},
         "stresses2": {"given": "stresses",
                          "calc": "stress",
                         "colours": "temperatures"},
        }


fig, ax = plt.subplots(1,len(plots))
try: 
    len(ax)
except:
    ax = [ax]
for pp in range(len(plots.keys())):
    p = list(plots.keys())[pp]
    ax[pp].set_title(p)
    ax[pp].plot([np.min(input_as_dict[plots[p]["given"]]),np.max(input_as_dict[plots[p]["given"]]) ], [np.min(input_as_dict[plots[p]["given"]]),np.max(input_as_dict[plots[p]["given"]]) ], 'k:')
    if "colours" in list(plots[p].keys()):
        if plots[p]["colours"] in list(input_as_dict):
            clrs = input_as_dict[plots[p]["colours"]]    
        else:
            clrs = all_recovered[plots[p]["colours"]]
        plttd = ax[pp].scatter(input_as_dict[plots[p]["given"]], all_recovered[plots[p]["calc"]], s=2, c = clrs, cmap="jet")
        cbar = plt.colorbar(plttd, shrink=0.6, location='bottom')
        cbar.set_label(plots[p]["colours"])
    else:
        plttd = ax[pp].scatter(input_as_dict[plots[p]["given"]], all_recovered[plots[p]["calc"]], s=2)
    ax[pp].set_xlabel(f"input: {plots[p]['given']}")
    ax[pp].set_ylabel(f"output: {plots[p]["calc"]}")
# fig.suptitle(f"{SampleGeometry} -- {SampleStrain}")