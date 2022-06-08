__all__ = ["Requirements", "WriteOutput"]

import os

import numpy as np
# import cpf.PeakFunctions as ff
import json
import cpf.IO_functions as IO
import cpf.series_functions as sf


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
        ]
    OptionalParams = [
        #"Output_ElasticProperties",  # FIX ME: this needs to be included -- not needed or multifit. only for WritePolyDefix
        #"Output_directory",  # if no direcrtory is specified write to current directory. -- not wrapped in the settings class and set by default
        "Output_NumAziWrite", # now set by default here therefore optional
    ]

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, differential_only=False, **kwargs):
def WriteOutput(setting_class=None,setting_file=None,differential_only=False, debug=True, **kwargs):
    
    # writes output from multifit in the form of *.fit files required for polydefix.
    # writes a separate file for each diffraction pattern.
    # uses the parameters in the json files to do so.
    
    
    if setting_class is None and setting_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif setting_class is None:
        import cpf.XRD_FitPattern.initiate as initiate
        setting_class = initiate(setting_file)
        
    
    #FitParameters = dir(FitSettings)

    # Parse the required inputs.
    #base_file_name = FitSettings.datafile_Basename

    # diffraction patterns
    # diff_files, n_diff_files = IO.file_list(FitParameters, FitSettings)
    # if "Output_directory" in FitParameters:
    #     out_dir = FitSettings.Output_directory
    # else:
    #     out_dir = "./"

    # def WriteOutput(base_file_name, data_to_write, Num_Azi, wavelength):
    ## writes *.fit files required by polydefix.

    # force Num_Azi to be a float
    # Num_Azi = float(Num_Azi)
    # Num_Azi = FitSettings.Output_NumAziWrite
    Num_Azi = 90
    if "Output_NumAziWrite" in setting_class.output_settings:
        Num_Azi = setting_class.output_settings["Output_NumAziWrite"]

    # wavelength = parms_dict["conversion_constant"]
    wavelength = setting_class.data_class.calibration["conversion_constant"]

    for z in range(setting_class.datafile_number):

        # read file to write output for
        # filename = os.path.splitext(os.path.basename(diff_files[z]))[0]
        # filename = filename+'.json'
        setting_class.set_subpattern(z,0)   
        
        filename = IO.make_outfile_name(
            # diff_files[z],
            # directory=FitSettings.Output_directory,
            setting_class.subfit_filename,#diff_files[z],
            directory=setting_class.output_directory,#directory=FitSettings.Output_directory,
            extension=".json",
            overwrite=True,
        )  # overwrite =false to get the file name without incrlemeting it.

        # Read JSON data from file
        with open(filename) as json_data:
            data_to_write = json.load(json_data)

        # create output file name from passed name
        base = setting_class.subfit_filename 
        if base is None:
            print("No base filename, using input filename instead.")
            base = os.path.splitext(os.path.split(setting_class.settings_file)[1])[0]
        # if not base:
        #     print("No base filename, using input filename instead.")
        #     base = os.path.splitext(os.path.split(FitSettings.inputfile)[1])[0]
        if differential_only is not False:
            base = base + "_DiffOnly"

        out_file = IO.make_outfile_name(
            # diff_files[z],
            # directory=FitSettings.Output_directory,
            base,#diff_files[z],
            directory=setting_class.output_directory,#directory=FitSettings.Output_directory,
            extension=".fit",
            overwrite=True,
        )
        # path, filename = os.path.split(diff_files[z])
        # base, ext = os.path.splitext(filename)
        # if differential_only is not False:
        #     base = base+'_DiffOnly'
        # out_file = out_dir + base + '.fit'

        # base, ext = os.path.splitext(os.path.split(diff_files[z])[1])
        # if not base:
        #     print("No base filename, using input filename instead.")
        #     base =  os.path.splitext(os.path.split(FitSettings.inputfile)[1])[0]

        print("Writing:", out_file)

        text_file = open(out_file, "w")

        # headers. set file version to be 1.
        text_file.write("# File version\n")
        text_file.write("1\n")

        # header peak profile
        # I dont know how used this is.
        # Also need to write a single number.
        if "profile" in data_to_write:
            to_write = data_to_write["profile"]
        else:
            to_write = 0

        text_file.write("# Peak profile (0: gauss, 1: pseudo-voigt, 2: lorentz)\n")
        text_file.write("%8i\n" % to_write)
        # FIX ME: this is currently fixed. It should be read from the input constraints.

        # number of subpatterns
        # FIX ME: I dont know how the data will be input here.
        # num_subpatterns = 2

        num_subpatterns = len(data_to_write)

        text_file.write("# Number of sub-patterns\n")
        text_file.write("%8i\n" % num_subpatterns)

        # print type(num_subpatterns)
        # print num_subpatterns

        for j in range(num_subpatterns):

            text_file.write("# Sub-patterns        %i\n" % j)

            # number of peaks
            text_file.write("# Number of peaks\n")
            text_file.write("%8i\n" % len(data_to_write[j]["peak"]))

            # peak profile
            text_file.write("# Peak profile\n")
            pro = data_to_write[j]["peak"][0]["profile"]
            if pro == 0:
                profile_write = 0
            elif pro == 1:
                profile_write = 2
            else:
                profile_write = 1
            text_file.write("%8i\n" % profile_write)

            # width factor
            text_file.write("# Width factor\n")
            text_file.write("%13.5f\n" % 15)
            # [0.000000]

            # Range in azimuth of each window/subpattern. Preceeded by the number of azimuths
            text_file.write("# Range in azimuth: nvalues, and values\n")
            # [number of azimuth slices]
            text_file.write("%10i\n" % Num_Azi)
            # [azimuth1
            # ...
            # azimuth last]
            for k in range(int(Num_Azi)):
                az = k / Num_Azi * 360
                text_file.write("%10.1f\n" % az)

            # Start and end of the back ground in each subpattern
            text_file.write("# background positions, in 2 theta (left/right)\n")
            # [bg1left      bg1right
            # ... These are the minimum and maximum extent of the background in 2 theta
            # BGlastLeft    BGlastRight]
            for k in range(int(Num_Azi)):
                text_file.write(
                    "%13.4f %13.4f\n"
                    % (data_to_write[j]["range"][0][0], data_to_write[j]["range"][0][1])
                )

            # value of the background intensity and slope (assuming a linear background)
            text_file.write("# background coefficients\n")
            for k in range(int(Num_Azi)):
                az = np.array([k / Num_Azi * 360])
                inter = sf.coefficient_expand(
                    az,
                    data_to_write[j]["background"][0],
                    data_to_write[j]["background_type"],
                )
                if len(data_to_write[j]["background"]) > 1:
                    slop = sf.coefficient_expand(
                        az,
                        data_to_write[j]["background"][1],
                        data_to_write[j]["background_type"],
                    )
                else:
                    slop = 0

                text_file.write("%13.4f %13.4f\n" % (inter, slop))

            # [bg1intensity      bg1slope
            # ... These are the intensity at the min 2 theta and slope the background against 2 theta      ]
            # BGlastIntensity    BGlastSlope]                                                              ]

            # header for the peaks in the file
            text_file.write(
                "# Peak positions (2 theta), intensity, half-width, weight Gauss/Lorentz (for pseudo-voigt),\n"
            )
            # write the peak properties for each azimuth slice. and each peak.
            num_peaks = len(data_to_write[j]["peak"])
            for k in range(num_peaks):
                text_file.write("# peak number        %i\n" % k)
                # [position1     Intensity1  HW1            Weight1            ]  repeat for the
                # ...                                                          ]  number of peaks in
                # positionLast  IntensityLast  HWLast      WeightLast]         ]  subpattern
                # print data_to_write[j]['peak'][k]['d-space']
                if "symmetry" in data_to_write[j]["peak"][k]:
                    sym = data_to_write[j]["peak"][k]["symmetry"]
                else:
                    sym = 1
                for l in range(int(Num_Azi)):
                    az = np.array([l / Num_Azi * 360])
                    if (
                        differential_only is False
                        or len(data_to_write[j]["peak"][k]["d-space"]) <= 3
                    ):
                        d_coef = data_to_write[j]["peak"][k]["d-space"]
                    else:
                        d_coef = data_to_write[j]["peak"][k]["d-space"]
                        d_coef = d_coef[:5]  # truncate array to 5 numbers.
                        d_coef[1] = 0  #
                        d_coef[2] = 0
                    peak_d = sf.fourier_expand((az), d_coef)
                    peak_tth = 2.0 * np.degrees(np.arcsin(wavelength / 2 / peak_d))
                    # peak_i = ff.Fourier_expand((az)*sym, data_to_write[j]['peak'][k]['height'])  #FIX ME - Is this the height of the peak or the integral under it?
                    # peak_w = ff.Fourier_expand((az)*sym, data_to_write[j]['peak'][k]['width'])   #FIX ME - is this the correct half width?
                    # peak_p = ff.Fourier_expand((az)*sym, data_to_write[j]['peak'][k]['profile']) #FIX ME - is this 1 or 0 for Gaussian?
                    peak_i = np.max(
                        [
                            0,
                            sf.coefficient_expand(
                                (az) * sym,
                                param=data_to_write[j]["peak"][k]["height"],
                                coeff_type=data_to_write[j]["peak"][k]["height_type"],
                            ),
                        ]
                    )  # FIX ME - Is this the height of the peak or the integral under it?
                    peak_w = sf.coefficient_expand(
                        (az) * sym,
                        param=data_to_write[j]["peak"][k]["width"],
                        coeff_type=data_to_write[j]["peak"][k]["width_type"],
                    )  # FIX ME - is this the correct half width?
                    peak_p = sf.coefficient_expand(
                        (az) * sym,
                        param=data_to_write[j]["peak"][k]["profile"],
                        coeff_type=data_to_write[j]["peak"][k]["profile_type"],
                    )  # FIX ME - is this 1 or 0 for Gaussian?

                    text_file.write(
                        "%13.4f %13.4f %13.4f %13.4f\n"
                        % (peak_tth, peak_i, peak_w, peak_p)
                    )

        text_file.close()


def WriteTestCase(FitSettings, parms_dict, differential_only=False):
    # writes output from multifit in the form of *.fit files required for polydefix.
    # writes a separate file for each diffraction pattern.
    # uses the parameters in the json files to do so.

    Fouriers_to_write = [
        {
            "peak": [
                {
                    "profile": 0,
                    "width": [0.03296637],
                    "d-space": [
                        1.45521770e02,
                        3.15812130e-03,
                        3.63949255e-03,
                        -1.43496474e02,
                        -1.43500658e02,
                    ],
                    "height": [
                        -1.61854953e08,
                        -7.90551608e01,
                        -1.86863466e02,
                        1.84245125e08,
                        -5.01600833e05,
                        8.75674029e02,
                        3.71825217e03,
                        -5.17425137e06,
                        3.65458422e08,
                        -3.90556882e03,
                        -2.25773282e04,
                        1.86844720e08,
                        -2.28359051e08,
                        8.22152351e03,
                        5.67280685e04,
                        -3.07040812e08,
                        9.14094789e06,
                        -8.05594617e03,
                        -6.26194032e04,
                        1.24648536e08,
                        -5.65124620e06,
                        2.94882628e03,
                        2.51878566e04,
                        -2.16682847e07,
                        2.17678521e07,
                    ],
                }
            ],
            "background": [4.0],
            "range": [[11.700006, 12.079979], [1.9937644, 2.0582786]],
        },
        {
            "peak": [
                {
                    "profile": 0,
                    "width": [0.06473334],
                    "d-space": [
                        -1.48371656e02,
                        1.32237912e-03,
                        1.75463473e-03,
                        1.49805560e02,
                        1.49801480e02,
                    ],
                    "height": [
                        2.05724056e07,
                        -9.25144706e00,
                        -5.38376022e00,
                        -2.10599065e07,
                        -6.44172806e06,
                        1.50563973e02,
                        9.46366683e01,
                        -6.81122757e05,
                        -3.99125724e07,
                        -7.06638114e02,
                        -4.55423181e02,
                        -6.71592503e06,
                        3.26600155e07,
                        1.31413614e03,
                        9.36273739e02,
                        -1.63116330e06,
                        3.30713638e06,
                        -1.04632954e03,
                        -9.05402148e02,
                        1.44412561e07,
                        -1.51098537e07,
                        2.97575768e02,
                        3.39659190e02,
                        -4.92552756e06,
                        4.92459547e06,
                    ],
                }
            ],
            "background": [4.0],
            "range": [[16.600002, 17.099998], [1.4110823, 1.4532731]],
        },
        {
            "peak": [
                {
                    "profile": 0,
                    "width": [0.05999371],
                    "d-space": [
                        -1.38852649e02,
                        8.75848195e-04,
                        1.30742495e-03,
                        1.40024571e02,
                        1.40021371e02,
                    ],
                    "height": [
                        -2.94680617e05,
                        -7.51087348e00,
                        -5.01955745e-01,
                        1.71972262e05,
                        3.01886458e06,
                        1.47013354e02,
                        -3.97875774e01,
                        1.00812116e06,
                        1.31668009e06,
                        -8.64958125e02,
                        4.82478179e02,
                        -3.80495454e07,
                        2.73111970e06,
                        2.04088806e03,
                        -1.76521792e03,
                        9.43916850e07,
                        -4.37268233e07,
                        -2.05402403e03,
                        2.45598540e03,
                        -8.07779850e07,
                        6.04951412e07,
                        7.38560951e02,
                        -1.13116306e03,
                        2.35504567e07,
                        -2.35403011e07,
                    ],
                }
            ],
            "background": [4.0],
            "range": [[20.400015, 20.89999], [1.1566441, 1.1846806]],
        },
    ]

    file_name = "Output1.txt"

    Azimuths = 360.0

    WriteMultiFit.WriteOutput(file_name, Fouriers_to_write, Azimuths)

    # Write the fits to a temporary file
    TempFilename = open("Fit_previous_JSON.dat", "w")

    # Write a JSON string into the file.
    json_string = json.dumps(
        Fouriers_to_write,
        TempFilename,
        sort_keys=True,
        indent=4,
        separators=(",", ": "),
    )
    TempFilename.write(json_string)
