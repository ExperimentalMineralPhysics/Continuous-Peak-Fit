#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 20:03:57 2025

@author: j033794sh
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 06:39:10 2025

@author: j033794sh
"""

# import cpf

# settings_class = cpf.XRD_FitPattern.initiate("CoSi22_MgO_input")


# #get data from settings class
# new_data = settings_class.data_class


# from pathlib import Path
# if settings_class.calibration_data:
#     data_to_fill = Path(settings_class.calibration_data).resolve()
# else:
#     data_to_fill = settings_class.image_list[0]

# new_data.fill_data(
#     data_to_fill,
#     settings=settings_class,
# )


# for n, i in enumerate(settings_class.image_list):

#     settings_class.set_subpattern(n,0)
#     # meta = new_data.get_metadata(settings=settings_class)
#     meta = new_data.get_metadata(settings=settings_class, metadata_values=['mean_start_time', 'mean_live_time', '6BMB_LVP:LVP_tc1_calcs.I', '6BMB_LVP:LVP_tc2_calcs.I', 'FILE_CREATION', 'time_start'])

#     print(meta)
