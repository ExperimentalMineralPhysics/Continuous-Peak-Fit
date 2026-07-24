
import glob

import fabio
import numpy as np
import numpy.ma as ma

from cpf.util.logging import get_logger

logger = get_logger("cpf.pattern_processing_methods.spots_scipy_peakfind")

import proglog
from pathlib import Path
from typing import Literal, Optional

from cpf.util.io import (
    has_value,
    make_outfile_name,
    numpy_to_json,
    peak_string,
    title_file_names,
)


from ImageD11.blobcorrector import perfect
from ImageD11.columnfile import columnfile
from ImageD11.labelimage import labelimage
from ImageD11.peaksearcher import peaksearch


"""
ths is probably a better way of finding peaks using ImageD111
Copied from ImageD11/test/peaksearchtiftest/scriptedpeaksearch.py
https://github.com/FABLE-3DXRD/ImageD11/blob/194d2fda453ee3e259e67651c3fcc1442dc3014b/test/peaksearchtiftest/scriptedpeaksearch.py
24th March 2023
"""

def method_defaults():
    """List non-universally required parameters for this image processing method."""
    # required parameters -- ones that cannot be guessed at in advance or given a pre-existing value
    # optional parameters -- values that can be guessed at in advance

    required_params = [
        #'apparently none!
    ]
    optional_params = {
        "threshhold": [1, 10, 100, 1000, 10000],
    }

    return required_params, optional_params




def spot_find(
        settings_class,
        data_class,
        **kwargs
        ):
    """
    Find bright spots in the diffration patterns using ImageD11. 

    Parameters
    ----------
    settings_class : TYPE
        DESCRIPTION.
    data_class : TYPE
        DESCRIPTION.
     : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # requies that submatter is set before calling this function -- so that can open the image again
    
    # Parse optional parameters
    threshold = settings_class.spot_find_settings.get("threshold", method_defaults()[1]["threshold"])

    # give somes minimal options
    corrector = perfect()  # no spatial disortion
    dims = data_class.intensity.shape
    label_ims = {
        t: labelimage(
            shape=dims,
            fileout=make_outfile_name(
                settings_class.subfit_filename,
                directory=settings_class.output_directory,
            )
            + "_t%d.flt" % (t),
            sptfile=make_outfile_name(
                settings_class.subfit_filename,
                directory=settings_class.output_directory,
            )
            + "_t%d.spt" % (t),
            spatial=corrector,
        )
        for t in threshold
    }

    # for filename, omega in lines:

    # open image with Fabio, required by peaksearch.
    # can't do this as a data class function because it cant pickle a Fabio instance.
    frame = fabio.open(settings_class.subfit_filename)
    frame.data = data_class.intensity

    frame.header["Omega"] = 0
    frame.data = frame.data.astype(np.float32)
    # corrections like dark/flat/normalise would be added here
    peaksearch(
        settings_class.subfit_filename.stem,
        frame,
        corrector,
        threshold,
        label_ims,
    )
    
    for t in threshold:
        label_ims[t].finalise()


def write_intermediate_files(filename, fits_to_write):
    """
    
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    fits_to_write : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    pass


def read_intermediate_files(
    settings_class,
    pattern = "all",
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs,
):
    """
    loads flt files for the selected subpattern.
    Makes an array of the peaks, where x,y are the centres and i is the size of the peak.
    Adds this array to the data class.

    Parameters
    ----------
    settings_file : TYPE, optional
        DESCRIPTION. The default is None.
    settings_class : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """

    print("Reading flt files")
    lgr = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger


    # find all the *.flt files for this image
    fnam = make_outfile_name(
        settings_class.subfit_filename, directory=settings_class.output_directory
    )
    fnams = glob.glob(fnam + "*.flt")
    print(fnams)


    # for f in range(setting_class.image_number):
    for f in lgr.iter_bar(iteration=range(settings_class.image_number)):
        settings_class.set_subpattern(f, 0)
        
        obj = []
        # obj = columnfile()
        for i in range(len(fnams)):
            try:
                obj.append(columnfile(fnams[i]))

                plt.scatter(-obj[i].dety, obj[i].detz, 1, c=(obj[i].sum_intensity))
                plt.colorbar()
                plt.title(fnams[i])
                plt.show()

                plt.scatter((obj[i].Number_of_pixels), (obj[i].sum_intensity))
                plt.title(fnams[i])
                plt.show()

            except:
                obj.append([])
        return obj

        
        
        filename = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            additional_text="chunks",
            extension=".json",
            overwrite=True,
        )
        if Path(filename).is_file():
            # Read JSON data from file
            with open(filename) as json_data:
                fit = json.load(json_data)
                all_azis.append(fit[1])
                all_fits.append(fit[0])

    print("Finished reading flt files")
    return all_azis, all_fits




    # if nothing has been set assume the first pattern
    if settings_for_fit.subfit_filename == None:
        settings_for_fit.set_subpattern(0, 0)

    # find all the *.flt files for this image
    fnam = make_outfile_name(
        settings_for_fit.subfit_filename, directory=settings_for_fit.output_directory
    )
    fnams = glob.glob(fnam + "*.flt")
    print(fnams)
    obj = []
    # obj = columnfile()
    for i in range(len(fnams)):
        try:
            obj.append(columnfile(fnams[i]))

            plt.scatter(-obj[i].dety, obj[i].detz, 1, c=(obj[i].sum_intensity))
            plt.colorbar()
            plt.title(fnams[i])
            plt.show()

            plt.scatter((obj[i].Number_of_pixels), (obj[i].sum_intensity))
            plt.title(fnams[i])
            plt.show()

        except:
            obj.append([])
    return obj



def write_grains_list(
    settings,
    
    prominence: int = 15,
    subpattern: str = "all",
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    
    ):
    """
    Write list of spots in the diffraction pattern, and their location and intensity.

    Parameters
    ----------
    settings : TYPE
        DESCRIPTION.
     : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """







# def load_spot_positions(
#     settings_file=None,
#     settings_class=None,
#     inputs=None,
#     report: Literal[
#         "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
#     ] = "INFO",
#     **kwargs,
# ):
#     """
#     loads flt files for the selected subpattern.
#     Makes an array of the peaks, where x,y are the centres and i is the size of the peak.
#     Adds this array to the data class.

#     Parameters
#     ----------
#     settings_file : TYPE, optional
#         DESCRIPTION. The default is None.
#     settings_class : TYPE, optional
#         DESCRIPTION. The default is None.

#     Returns
#     -------
#     None.

#     """

#     if settings_class is None:
#         settings_for_fit = initiate(settings_file, inputs=inputs, report=report)
#     else:
#         settings_for_fit = settings_class

#     # if nothing has been set assume the first pattern
#     if settings_for_fit.subfit_filename == None:
#         settings_for_fit.set_subpattern(0, 0)

#     # find all the *.flt files for this image
#     fnam = make_outfile_name(
#         settings_for_fit.subfit_filename, directory=settings_for_fit.output_directory
#     )
#     fnams = glob.glob(fnam + "*.flt")
#     print(fnams)
#     obj = []
#     # obj = columnfile()
#     for i in range(len(fnams)):
#         try:
#             obj.append(columnfile(fnams[i]))

#             plt.scatter(-obj[i].dety, obj[i].detz, 1, c=(obj[i].sum_intensity))
#             plt.colorbar()
#             plt.title(fnams[i])
#             plt.show()

#             plt.scatter((obj[i].Number_of_pixels), (obj[i].sum_intensity))
#             plt.title(fnams[i])
#             plt.show()

#         except:
#             obj.append([])
#     return obj


def make_im_from_flts(
    settings_file=None,
    settings_class=None,
    data_class=None,
    inputs=None,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    debug=False,
    **kwargs,
):
    if settings_class is None:
        settings_for_fit = initiate(settings_file, inputs=inputs, report=report)
    else:
        settings_for_fit = settings_class

    # if nothing has been set assume the first pattern
    if settings_for_fit.subfit_filename == None:
        settings_for_fit.set_subpattern(0, 0)

    if data_class == None:
        data_class = settings_for_fit.data_class
        data_class.fill_data(
            settings_for_fit.subfit_filename, settings=settings_for_fit
        )

    peaks_im = ma.zeros(data_class.intensity.shape)

    pks = load_flts(
        settings_file=settings_file,
        settings_class=settings_for_fit,
        data_class=data_class,
    )

    for i in range(len(pks)):
        try:
            x = -pks[i].dety
            y = pks[i].detz
            z = pks[i].sum_intensity

            for j in range(len(x)):
                peaks_im[int(y[j]), int(x[j])] = z[j]
        except:
            pass

    # peaks_im = peaks_im(mask=data_class.intensity.mask)
    peaks_im = ma.masked_where(ma.getmask(data_class.intensity), peaks_im)

    if debug:
        fig = plt.figure()
        plt.imshow(peaks_im, vmax=30)
        plt.colorbar()
        plt.show()

        fig = plt.figure()
        plt.imshow(np.log10(data_class.intensity))
        plt.colorbar()
        plt.show()

        fig = plt.figure()
        plt.scatter(data_class.tth, data_class.azm, s=0.1, c=peaks_im, vmax=30)
        plt.colorbar()
        plt.show()

    # print(pks[0].get_bigarray)
    # print(dir(pks[0]))

    data_class.peaks_image = peaks_im

    return data_class




# %% =================================================
# all dead code

# =================================================


def execute2(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug=False,
    save_all=False,
    parallel=True,
    subpattern="all",
    mode="cascade",
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    threshold=[1, 10, 100, 1000, 10000],
    **kwargs,
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param settings_file:
    :param parallel:
    :param report:
    :param mode:
    :param track:
    :param propagate:
    :param save_all:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """

    if settings_class is None:
        settings_for_fit = initiate(settings_file, inputs=inputs, report=report)
    else:
        settings_for_fit = settings_class
    new_data = settings_for_fit.data_class

    if settings_for_fit.calibration_data:
        data_to_fill = settings_for_fit.calibration_data.resolve()
    else:
        data_to_fill = settings_for_fit.image_list[0].resolve()

    new_data.fill_data(
        data_to_fill,
        settings=settings_for_fit,
        debug=debug,
    )

    # plot calibration file
    if debug and settings_for_fit.calibration_data is not None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        new_data.plot_collected(fig_plot=fig, axis_plot=ax)
        plt.title("Calibration data")
        plt.show()
        plt.close()

    if not isinstance(threshold, list):
        threshold = [threshold]

    # if parallel processing start the pool
    # if parallel is True:
    #     p = mp.Pool(processes=mp.cpu_count())
    # p = mp.Pool()

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

    # Process the diffraction patterns #
    for j in range(settings_for_fit.datafile_number):
        print("Process ", settings_for_fit.datafile_list[j])
        # Get diffraction pattern to process.
        new_data.import_image(settings_for_fit.datafile_list[j], debug=debug)

        if settings_for_fit.datafile_preprocess is not None:
            # needed because image preprocessing adds to the mask and is different for each image.
            new_data.mask_restore()
            if "cosmics" in settings_for_fit.datafile_preprocess:
                new_data = cosmicsimage_preprocess(new_data, settings_for_fit)
        else:
            # nothing is done here.
            pass

        # plot input file
        if debug:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax_o1 = plt.subplot(111)
            new_data.plot_calibrated(fig_plot=fig, axis_plot=ax, show="intensity")
            plt.title(settings_for_fit.datafile_list[j].stem)
            plt.show()
            plt.close()

        # give somes minimal options
        corrector = perfect()  # no spatial disortion
        dims = new_data.intensity.shape
        label_ims = {
            t: labelimage(
                shape=dims,
                fileout=make_outfile_name(
                    settings_for_fit.datafile_list[j],
                    directory=settings_for_fit.output_directory,
                )
                + "_t%d.flt" % (t),
                sptfile=make_outfile_name(
                    settings_for_fit.datafile_list[j],
                    directory=settings_for_fit.output_directory,
                )
                + "_t%d.spt" % (t),
                spatial=corrector,
            )
            for t in threshold
        }

        # for filename, omega in lines:

        # open image with Fabio, required by peaksearch.
        # can't do this as a data class function because it cant pickle a Fabio instance.
        frame = fabio.open(settings_for_fit.datafile_list[j])
        frame.data = new_data.intensity

        frame.header["Omega"] = 0
        frame.data = frame.data.astype(np.float32)
        # corrections like dark/flat/normalise would be added here
        peaksearch(
            settings_for_fit.datafile_list[j].stem,
            frame,
            corrector,
            threshold,
            label_ims,
        )
        for t in threshold:
            label_ims[t].finalise()

    # if debug:





"""
# ============================
# NOTES:
# These are a set of functions that attempted to find the spots in the raw data and count them.
# They were either very slow or didn't work.
# They are here to remind me to implement them at some time in the future because they are a 'better' way of making the measurements.
# ============================

def PeakCount(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, subpattern='all', **kwargs):


    FitSettings, FitParameters, new_data = XRDFit.initiate(settings_file, inputs=inputs, report=True)

    # Get list of diffraction patterns #
    diff_files, n_diff_files = IO.FileList(FitParameters, FitSettings)

    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60

    # restrict to subpatterns listed
    if subpattern=='all':
        subpats = list(range(0, len(FitSettings.fit_orders)))
    elif isinstance(subpattern,list):
        subpats = subpattern
    else:
        subpats = [int(x) for x in str(subpattern)]

    # make new order search list
    orders_tmp = []
    for i in range(len(subpats)):
        j = subpats[i]
        orders_tmp.append(FitSettings.fit_orders[j])

    FitSettings.fit_orders = orders_tmp


    FitSettings_files = copy.deepcopy(FitSettings.datafile_Files)
    num_peaks_all = []
    visible = []

    for h in range(n_diff_files):

        #logger.info(" ".join(map(str, [(FitParameters)])))
        #logger.info(" ".join(map(str, [(diff_files)])))

        #FitSettings_tmp = FitSettings
        FitSettings_tmp = FitSettings
        FitSettings_tmp.datafile_Files = [FitSettings_files[h]]
        #FitSettings_tmp.datafile_Files[0] = diff_files[h]

        #logger.info(" ".join(map(str, [(FitSettings_tmp.datafile_Files)])))
        #logger.info(" ".join(map(str, [(FitSettings_tmp)])))
        #logger.info(" ".join(map(str, [(FitParameters)])))

        outfname = IO.make_outfile_name(diff_files[h], directory=FitSettings_tmp.Output_directory, extension='csv', additional_text='peaks')
        logger.info(" ".join(map(str, [(outfname)])))

        if not os.path.exists(outfname):
            #This currently gets all the peaks in the image, not just those in the regios of interest
            intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data,
                    debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                    parallel=parallel,
                    mode='ReturnImage', report=True)

            # initialize with default parameters. The "denoise" parameter can be of use in your case
            # import 2D example dataset

            plt.imshow(twotheta)
            plt.colorbar()

            plt.imshow(np.array(twotheta))
            plt.colorbar()

            img = ma.array(intens)#, mask=False)
            plt.imshow(img)

            img.filled(0)
            img = ma.array(img, mask = None)
            plt.imshow(img)

            fp = findpeaks(limit=1, whitelist=['peak'], scale=False, togray=False, denoise=None)

            # make the fit
            fp.fit(img)
            # Make plot
            #fp.plot()
            #fp.plot_persistence()

            fp.results['persistence']
            logger.info(" ".join(map(str, [(fp.results['persistence'])])))
            logger.info(" ".join(map(str, [(type(fp.results['persistence']))])))

            peaks = fp.results
            #logger.info(" ".join(map(str, [(type(peaks))])))
            #logger.info(" ".join(map(str, [(peaks)])))
            #logger.info(" ".join(map(str, [(type(peaks))])))

            fp.results['persistence'].to_csv(outfname)
        else:

            peaks=pd.read_csv(outfname)
            #logger.info(" ".join(map(str, [(fp)])))
            #peaks = fp.results
            if h==0:
                #This currently gets all the peaks in the image, not just those in the regios of interest
                intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data,
                    debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                    parallel=parallel,
                    mode='ReturnImage', report=True)

        # get peaks in each diffraction peak window

        num_peaks = []
        for i in range(len(FitSettings.fit_orders)):

            tthRange = FitSettings.fit_orders[i]['range'][0]

            msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])


            # find and list all peaks in unmasked area.
            #FIXME: how does this interact with the region specific masks?
            pks=[]
            for j in range(len(peaks)):
                x = peaks['x'].loc[j]
                y = peaks['y'].loc[j]
                if twotheta.mask[y][x]==False:
                    if ((twotheta[y][x] > tthRange[0]) & (twotheta[y][x] < tthRange[1]) & (twotheta.mask[y][x] == False)):
                        pks.append(peaks.loc[j])
            num_peaks.append(len(pks))

            #approximate fraction of ring that is masked.
            #number of pixels in two theta range
            num_pix = ((np.array(twotheta) > tthRange[0]) & (np.array(twotheta) < tthRange[1])).sum()

            #number of unmasked pixles in two theta range
            num_vis_pix = ((tthRange[0] < twotheta) & (twotheta < tthRange[1])).sum()


            visible.append(num_vis_pix/num_pix)

            # # find and list all peaks in unmasked area.
            # #FIXME: how does this interact with the region specific masks?
            # peaks = []
            # for j in range(len(fp.results['persistence'])):
            #     x = fp.results['persistence']['x'].loc[j]
            #     y = fp.results['persistence']['y'].loc[j]
            #     if twotheta.mask[y][x]==False:
            #         if ((twotheta[y][x] > tthRange[0]) & (twotheta[y][x] < tthRange[1]) & (twotheta.mask[y][x] == False)):
            #             peaks.append(fp.results['persistence'].loc[j])
            # num_peaks.append(len(peaks))

            # #approximate fraction of ring that is masked.
            # #number of pixels in two theta range
            # num_pix = ((np.array(twotheta) > tthRange[0]) & (np.array(twotheta) < tthRange[1])).sum()

            # #number of unmasked pixles in two theta range
            # num_vis_pix = ((tthRange[0] < twotheta) & (twotheta < tthRange[1])).sum()


            # visible.append(num_vis_pix/num_pix)
        num_peaks_all.append(num_peaks)
        #logger.info(" ".join(map(str, [(num_peaks)])))
        #logger.info(" ".join(map(str, [(visible)])))

    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60
    num_peaks_all=np.array(num_peaks_all)
    logger.info(" ".join(map(str, [(num_peaks_all)])))
    logger.info(" ".join(map(str, [(modified_time)])))
    if 1:
        fig,ax = plt.subplots()
        for i in range(len(FitSettings.fit_orders)):
            logger.info(" ".join(map(str, [(num_peaks_all[:,i])])))
            plt.plot(modified_time, num_peaks_all[:,i], '.-', label=IO.peak_string(FitSettings.fit_orders[i]))

        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Intensity Maxima Count')
        plt.legend()
        plt.show()

        fig.savefig('IntensityMaximaTime2.png')
        fig.savefig('IntensityMaximaTime2.pdf')
        fig.savefig('IntensityMaximaTime2.eps')

    #plt.scatter([*range(len(num_peaks_all))], num_peaks_all)
    #plt.scatter

        # #fraction that is unmasked (as approximation for how much of ring is visible.)

        # tthRange = FitSettings.fit_orders[i]['range'][0]

        # msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])

        # img = ma.array(intens, mask=msk.mask)

        # plt.imshow(img)

        # img.filled(0)
        # img = ma.array(img, mask = None)
        # plt.imshow(img)

        # fp = findpeaks(limit=2, whitelist=['peak'])

        # # make the fit
        # fp.fit(img)
        # # Make plot

        # fp.plot()

        # fp.plot_persistence()

        # fp.results['persistence']

        # logger.info(" ".join(map(str, [(fp.results)])))










    for i in range(len(FitSettings.fit_orders)):


        tthRange = FitSettings.fit_orders[i]['range'][0]

        msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])

        img = ma.array(intens, mask=msk.mask)

        plt.imshow(img)

        img.filled(0)
        img = ma.array(img, mask = None)
        plt.imshow(img)

        fp = findpeaks(limit=2, whitelist=['peak'])

        # make the fit
        fp.fit(img)
        # Make plot

        fp.plot()

        fp.plot_persistence()

        fp.results['persistence']

        logger.info(" ".join(map(str, [(fp.results)])))







    if 1:
        fig,ax = plt.subplots()

        for i in range(len(FitSettings.fit_orders)):
            plt.plot(modified_time, peak_count[:][i], '.-', label=IO.peak_string(FitSettings.fit_orders[i]))

        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Intensity Maxima Count')
        plt.legend()
        plt.show()

        fig.savefig('IntensityMaximaTime.png')
        fig.savefig('IntensityMaximaTime.pdf')
        fig.savefig('IntensityMaximaTime.eps')



def WatershedCount(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, subpattern='all', **kwargs):


    FitSettings, FitParameters, new_data = XRDFit.initiate(settings_file, inputs=inputs, report=True)

    # Get list of diffraction patterns #
    diff_files, n_diff_files = IO.FileList(FitParameters, FitSettings)

    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60

    # restrict to subpatterns listed
    if subpattern=='all':
        subpats = list(range(0, len(FitSettings.fit_orders)))
    elif isinstance(subpattern,list):
        subpats = subpattern
    else:
        subpats = [int(x) for x in str(subpattern)]

    # make new order search list
    orders_tmp = []
    for i in range(len(subpats)):
        j = subpats[i]
        orders_tmp.append(FitSettings.fit_orders[j])

    FitSettings.fit_orders = orders_tmp


    FitSettings_files = copy.deepcopy(FitSettings.datafile_Files)
    num_peaks_all = []
    visible = []

    for h in range(n_diff_files):

        #This currently gets all the peaks in the image, not just those in the regios of interest
        intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data,
                debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                parallel=parallel,
                mode='ReturnImage', report=True)

        im = img_as_float(intens)

        # image_max is the dilation of im with a 20*20 structuring element
        # It is used within peak_local_max function
        image_max = ndi.maximum_filter(im, size=20, mode='constant')

        # Comparison between image_max and im to find the coordinates of local maxima
        coordinates = peak_local_max(im, min_distance=20)

        # display results
        fig, axes = plt.subplots(1, 3, figsize=(8, 3), sharex=True, sharey=True)
        ax = axes.ravel()
        ax[0].imshow(im, cmap=plt.cm.gray)
        ax[0].axis('off')
        ax[0].set_title('Original')

        ax[1].imshow(image_max, cmap=plt.cm.gray)
        ax[1].axis('off')
        ax[1].set_title('Maximum filter')

        ax[2].imshow(im, cmap=plt.cm.gray)
        ax[2].autoscale(False)
        ax[2].plot(coordinates[:, 1], coordinates[:, 0], 'r.')
        ax[2].axis('off')
        ax[2].set_title('Peak local max')

        fig.tight_layout()

        plt.show()





"""
