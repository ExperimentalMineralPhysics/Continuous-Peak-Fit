
import cpf
from cpf.output_formatters.fits_io import ReadFits_to_dataframe

settings_class = cpf.XRD_FitPattern.initiate("CoSi22_MgO_input")

#get data from settings class
new_data = settings_class.data_class

settings_class.metadata = ['mean_start_time', 'mean_live_time', '6BMB_LVP:LVP_tc1_calcs.I', '6BMB_LVP:LVP_tc2_calcs.I', 'FILE_CREATION', 'frames start']

# read data from individual files
for n, i in enumerate(settings_class.image_list):
    
    settings_class.set_subpattern(n,0)
    new_data.fill_data(
        settings=settings_class,
    )
    meta = new_data.get_metadata()
    
    print("FILE:", settings_class.subfit_filename)
    print(meta)
    print("\n")
    

# read all fits with metadata (may require executing 'CoSi22_MgO_input' first. 
df = ReadFits_to_dataframe(settings_class)

print(df)
