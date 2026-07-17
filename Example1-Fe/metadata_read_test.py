
import cpf



settings_class = cpf.XRD_FitPattern.initiate("BCC1_Dioptas_input")


#get data from settings class
new_data = settings_class.data_class


for n, i in enumerate(settings_class.image_list):
    settings_class.set_subpattern(n,0)
    new_data.fill_data(
        settings=settings_class,
    )
    meta = new_data.get_metadata(settings_class=settings_class, metadata_values=['FILE_CREATION', 'FILE_MODIFIED', 'Exposure_time', 'Exposure_period', 'frames start', 'frames end', 'namename'])
    print(meta)

print("   ")

settings_class.metadata = ['FILE_CREATION', 'FILE_MODIFIED', 'Exposure_time', 'Exposure_period', 'frames start', 'frames end', 'namename']
# settings_class.metadata_settings = {'exposure': 'Exposure_period'}
for n, i in enumerate(settings_class.image_list):
    settings_class.set_subpattern(n,0)
    new_data.fill_data(
        settings=settings_class,
    )
    meta = new_data.get_metadata(settings_class=settings_class)
    print(meta)
