__all__ = ["outfile_version", "make_header", "write_csv", "csv_align_columns", "float_format"]


import numpy as np

from cpf.util.logging import get_logger
logger = get_logger("cpf.output_formatters.output_csv")


def outfile_version():
    """ returns outfile version
    
    v1: initial go
    v2: seprate version in all files
    v3: rewritten using underlying panda data frame - to force columns to line up.
    v4: writing of csvs is all from single function/file. 
    """
    f_version = 4
    return f_version


def make_header(settings_class, fits=False, derived=False, calc_options=None, additional=None, **kwargs):
    """
    Makes header for csv based output files

    if kwargs has "comment" this is treated as the comment denoter for the output file
    the default is "#". 

    Parameters
    ----------
    settings_class : cpf settings class
        Settings class for cpf.
    fits : bool, string, optional
        True or string describing types of fit being writtend. The default is False.
    derived : bool, string, optional
        True or string describing derived values being writtend. The default is False.
    additional : list of strings, optional
        Extra information to be written to the header. The default is None.
    calc_options : dict
        Dictionary of kay, value pair parameters to be written into the header.
    **kwargs :

    Returns
    -------
    file_header : list of strings
        Header to be written to the file.

    """
    comment = kwargs.pop("comment", "#")
        
    file_header = []
    
    if fits:
        if fits == True:
            fits = "Fits"
        file_header += [f"{comment} {fits} for input file/settings: {settings_class.settings_file}.\n"]
    elif derived:
        file_header += [f"{comment} {derived} derived from fits made using input file/settings: {settings_class.settings_file}.\n"]
    else:
        file_header += [f"{comment} Input file/settings: {settings_class.settings_file}.\n"]
        
    # FIXME: this line should include the version 
    file_header += [f"{comment} using *Continuous Peak Fit*. (https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit) \n",
                    f"{comment} outfile version: {outfile_version()} \n",
                    f"{comment} \n"]
    
    # loop over kawrgs and add a line per value.
    if calc_options:
        file_header += [f"{comment} *** Settings/options used: *** \n"]
        file_header += [
            ''.join(f"{comment} {key}: {value}\n" for key, value in calc_options.items()),
            f"{comment} \n",]
    #for pair in calc_options:
        # file_header += [
        #     "{comment} ",
        #     ', '.join(f"{pair}: {calc_options[pair]}" for key, value in pair.items()),
        #     "\n",]
        # file_header += ["{comment} \n"]
    
    if additional:
        for i in additional:
            if i[0] != comment:
                i = f"{comment} {i}\n"
            file_header += [i]
        file_header += [f"{comment} \n"]
    
    
    return file_header
        


def write_csv(out_file, df, column_headers, file_header=None, col_width=15, dp=5):
    """
    Writes a formatted csv file from panda data frame. 
    
    Forces all the column widths so that the csv files looks sensible. 

    Parameters
    ----------
    out_file : str, Path
        filename or Path to be written to.
    df : panda data frame
        DataFrame to be written to csv file.
    column_headers : list
        Which columns to write from dataframe to csv file.
    file_header : list, str, optional
        header to be written to the file. The default is None.
    col_width : int, optional
        column width for aligned columns in csv. The default is 15.
    dp : int, optional
        number of decimal places in written numbers. The default is 5.

    Returns
    -------
    None.

    """
    logger.info(" ".join(map(str, [("Writing %s" % out_file)])))
        
    ## format dateframe for writing to file neatly.
    # make strings in DateFile and Peak columns all the same length
    if "DataFile" in df:
        len_datafile = df["DataFile"].str.len().max()
        df_tmp = df["DataFile"].str.pad(
            np.max([len_datafile, col_width]), side="left", fillchar=" "
        )
        df["DataFile"] = df_tmp
    if "peak" in df:    
        len_peaks = df["peak"].str.len().max()
        df_tmp = df["peak"].str.pad(
            np.max([len_peaks, 5]), side="left", fillchar=" "
        )
        df["peak"] = df_tmp
    
    # make sure there are no commas in any entry 
    # some date formats might have them in.
    for col in list(df):
        df[col] = df[col].replace(',',' ', regex=True)
    
    #shorten hdf5 key names 
    for col in list(df):
        if "/" in col:
            df.rename(columns={col: col.split("/")[-1]})

    # make sure diferent columns are saved as desired.
    for i in df.columns:
        if ("date" in i.lower() or 
            "time" in i.lower() or 
            i.lower() == "FILE_CREATION".lower() or 
            i.lower() == "FILE_MODIFIED".lower() 
            ):
            # the a date or time so need to keep all precision.
            # convert to string so that float_format is passed over.
            # if treated as a float then unix time looses all the precision.
            if df[i].dtypes != "O":
                #then not object type and can assume is a number
                df[i] = df[i].apply(lambda x: f"{x: {np.max([18, col_width])}.{dp}f}")
            else:
                pass #can assume is string?
        elif 'residuals' in i.lower() and df[i].dtypes == "O":
            # make sure residual columns are saved as a single string with no line breaks.
            if df[i].dtypes == "O":
                #then object type column and can assume is a list
                df[i] = df[i].apply(lambda x: np.array2string(x, separator=";", max_line_width=np.inf, formatter={"float_kind": lambda x: float_format(x, np.min([12, col_width]), dp) }, sign=" "))
        

    # rename the columns so that the headers are the same width as the columns
    class NewClass(object):
        pass

    columns = NewClass()
    for i in range(len(column_headers)):
        if column_headers[i] == "DataFile":
            setattr(
                columns,
                column_headers[i],
                column_headers[i].rjust(np.max([len_datafile, col_width])),
            )
        elif column_headers[i] == "pdatafileeak":
            setattr(
                columns,
                column_headers[i],
                column_headers[i].rjust(np.max([len_peaks, col_width])),
            )
        else:
            setattr(columns, column_headers[i], column_headers[i].rjust(col_width))
    columns = columns.__dict__
    df.rename(columns=columns, inplace=True)

    # write data frame to csv file
    with open(out_file, "w") as f:
        if isinstance(file_header, list):
            for i in file_header:
                f.write(i)
        else:
            f.write(file_header)
        df.to_csv(
            f,
            index=False,
            header=True,
            na_rep="nan",#.rjust(col_width),
            float_format = lambda x: float_format(x, col_width, dp)
                )

    # rewrite the file adjusting the column widths to keep the data lined up. 
    csv_align_columns(out_file, col_width=col_width, dp=dp)


def csv_align_columns(file, col_width=15, dp=5):
    """
    Rewrite the file adjusting the column widths to keep the data lined up. 

    Parameters
    ----------
    out_file : str, Path
        File to be rewritten.
    col_width : int, optional
        Column width to work to. The default is 15.
    dp : int, optional
        Number of decimal places to write. The default is 5.
    """
    with open(file, 'r') as fl: 
        in_lines = fl.readlines() 
    with open(file, 'w') as f:
        for line in in_lines:
            split_line = line.split(",")
            split_line_out = []
            running_length = 0
            expected_length = 0
            for i in range(len(split_line)):
                
                split_line[i] = split_line[i].replace(f"{0:.{dp}e}", "0".rjust(len(f"{0:.{dp}e}")))
                
                if i==0: # first column. keep narrow
                    col_here = 3
                    if "num" in split_line[i]:
                        split_line_out.append(f"{split_line[i].replace(' ',''):>{col_here}}")
                    else:
                        split_line_out.append(f"{split_line[i]:>{col_here}}")
                elif i==1: # file names
                    col_here = col_width
                    split_line_out.append(f" {split_line[i]:>{col_here}}")
                else: # data values. adjust column width to line everything up.
                    col_here = col_width - (running_length-expected_length)
                    if col_here >= 0:
                        split_line_out.append(f"{split_line[i].strip():>{col_here}}")
                    else: 
                        split_line_out.append(f"{split_line[i].strip()}")
                    running_length += len(split_line_out[-1])
                    expected_length += col_width
                
            if split_line_out[-1][-1] != "\n":
                split_line_out[-1] += "\n"
            out_line = ",".join(split_line_out) 
            f.write(out_line)
    

def float_format(x, col_width=15, dp=5):
    """
    Formaat floats for writing to csv

    Parameters
    ----------
    x : number
        number to be formatted.
    col_width : int, optional
        Column width to work to. The default is 15.
    dp : int, optional
        Number of decimal places to write. The default is 5.

    Returns
    -------
    formatted number as string

    """
    if np.int32(x) == np.float64(x):
        out = f"{x:.0f}"
    elif np.abs(x) > 0.1 and np.log10(np.abs(x)) < col_width-dp-4:
        out = f"{x:{str(col_width)}.{str(dp)}f}"
    else:
        out = f"{x:{str(col_width)}.{str(dp)}e}"
    return out

    
    
    