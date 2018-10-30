
import os


def WriteMultiFit(base_file_name, data_to_write):



    text_file = open("Output.txt", "w")

    #headers. set file version to be 1.
    text_file.write("# File version\n" )
    text_file.write("1\n" )

    #header peak profile
    # I dont know how used this is.
    #Also need to write a single number.
    if 'profile' in data_to_write:
        to_write = data_to_write['profile']
    else:
        to_write = 0

    text_file.write("# Peak profile (0: gauss, 1: pseudo-voigt, 2: lorentz)\n" )
    text_file.write("      %i\n" % to_write)
    
    #number of subpatterns
    # FIX ME: I dont know how the data will be input here.
    num_subpatterns = 2

    text_file.write("# Number of sub-patterns\n" )
    text_file.write("      %i\n" % num_subpatterns)

    print type(num_subpatterns)
    print num_subpatterns
    
    for j in range(num_subpatterns):

        text_file.write("# Sub-patterns        %i\n" % j)


        #number of peaks
        text_file.write("        # Number of peaks\n")
        #[1]

        #peak profile
        text_file.write("        # Peak profile\n")
        #[1]

        #width factor
        text_file.write("        # Width factor\n")
        #[0.000000]

        #Range in azimuth of each window/subpattern. Preceeded by the number of sub patterns
        text_file.write("        # Range in azimuth: nvalues, and values\n")
        #[number of azimuth slices]
        #[azimuth1
        #...
        #azimuth last]

        #Start and end of the back ground in each subpattern 
        text_file.write("        # background positions, in 2 theta (left/right)\n")
        #[bg1left      bg1right
        #... These are the minimum and maximum extent of the background in 2 theta
        #BGlastLeft    BGlastRight]

        # value of the background intensity and slope (assuming a linear background)
        text_file.write("        # background coefficients\n")
        #[bg1intensity      bg1slope
        #... These are the intensity at the min 2 theta and slope the background against 2 theta      ]
        #BGlastIntensity    BGlastSlope]                                                              ]

        #header for the peaks in the file
        text_file.write("        # Peak positions (2 theta), intensity, half-width, weight Gauss/Lorentz (for pseudo-voigt),\n")

        # write the peak properties for each azimuth slice. and each peak.
        num_peaks = 1
        for k in range(num_peaks):
            text_file.write("        # peak number        %i\n" % k)
            #[position1     Intensity1  HW1            Weight1            ]  repeat for the
            #...                                                          ]  number of peaks in
            #positionLast  IntensityLast  HWLast      WeightLast]         ]  subpattern




    text_file.close()