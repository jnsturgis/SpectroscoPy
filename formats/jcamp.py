#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:08:22 2024

@author: James STURGIS

This file implements import and export of spectra with JCAMP DX format files
for my spectroscopy programmes and utilities.

Much of the code comes from Nathan Hagen and is used under the MIT/X11 Licence
the original code can be found here: https://github.com/nzhagen/jcamp

The original code uses as a memory structure a dictionary, this version is
modified to use the spectrum class for internal storage.

TODO: Modify return type and define it in separate module

The module implements 2 public functions:
    jcamp.read(filehandle, my_spectrum) and
    jcamp.write( filehandle, my_spectrum)
"""

import datetime
import re

import numpy as np
import spectroscopy as spc

## In SQZ_digits, '+' or '-' is for PAC, ',' for CSV.
SQZ_digits = {'@':'+0', 'A':'+1', 'B':'+2', 'C':'+3', 'D':'+4', 'E':'+5',
              'F':'+6', 'G':'+7', 'H':'+8', 'I':'+9', 'a':'-1', 'b':'-2',
              'c':'-3', 'd':'-4', 'e':'-5', 'f':'-6', 'g':'-7', 'h':'-8',
              'i':'-9', '+':'+',  '-':'-',  ',':' '}
DIF_digits = {'%': 0, 'J':1,  'K':2,  'L':3,  'M':4,  'N':5,  'O':6,  'P':7,
              'Q':8,  'R':9,  'j':-1, 'k':-2, 'l':-3, 'm':-4, 'n':-5, 'o':-6,
              'p':-7, 'q':-8, 'r':-9}
DUP_digits = {'S':1, 'T':2, 'U':3, 'V':4, 'W':5, 'X':6, 'Y':7, 'Z':8, 's':9}

## The specification allows multiple formats for representing LONGDATE.
## See `FRACTIONAL_SECONDS_PATTERN` below for the optional token representing
## fractional seconds. These fractional seconds are removed in advance. Thus
## `%N` is not referenced in the formats below.
DATE_FORMATS = ["%Y/%m/%d %H:%M:%S %z", "%Y/%m/%d %H:%M:%S", "%Y/%m/%d"]

## The optional token describing the fractional seconds is referenced in the
## specification as `.SSSS`. This number of digits (four) is rather unclear,
## since the usual presentation of a fraction of seconds would contain
## either 3, 6 or 9 digits.
FRACTIONAL_SECONDS_PATTERN = re.compile(
    r"^\d{4}/\d{2}/\d{2} +\d{2}:\d{2}\d{2}(?P<fractional_seconds>\d{1,9})"
)

##=============================================================================
def parse_longdate(date_string: str) -> datetime.datetime:
    """parse the "LONGDATE" field according to the JCAMP-DX specification

    raises ValueError in case of problems
    """
    fractional_seconds_match = FRACTIONAL_SECONDS_PATTERN.search(date_string)
    if fractional_seconds_match:
        ## Remove the fractional seconds string - to simplify `strptime`.
        date_string = FRACTIONAL_SECONDS_PATTERN.sub("", date_string)

        ## Try to interprete the fractional seconds. The JCAMP specification
        ## (v6.00) does not explain, how a string of arbitrary length is
        ## supposed to be interpreted. Thus we are just guessing based on the
        ## number of digits.
        fraction_secs = \
            fractional_seconds_match.group("fractional_seconds")
        if len(fraction_secs) in {7, 8, 9}:
            ## This is probably nanoseconds.
            microseconds = int(int(fraction_secs) / 1000)
        elif len(fraction_secs) in {4, 5, 6}:
            microseconds = int(fraction_secs)
        elif len(fraction_secs) in {1, 2, 3}:
            microseconds = 1000 * int(fraction_secs)
        else:
            ## We should never end up here.
            raise ValueError(f"Fractional seconds string could not be parsed: {fraction_secs}")
    else:
        microseconds = 0

    ## Parse the date and time.
    for fmt in DATE_FORMATS:
        try:
            parsed = datetime.datetime.strptime(date_string, fmt)
        except ValueError:
            pass
        else:
            ## Inject the previously parsed microseconds
            return parsed.replace(microsecond=microseconds)
    raise ValueError(f"Failed to parse the date string: {date_string}.")

##=============================================================================

def read(filehandle, my_spectrum ) -> None:
    '''
    Read a JDX-format file, and update the spectrum my_spectrum

    Parameters
    ----------
    filehandle : file object
        The object representing the JCAMP-DX filename to read.
    my_spectrum : a spectrum object
        Where to store the information

    Returns
    -------
    nothing.
    '''

    jcamp_dict = {}
    y = []
    x = []
    datastart = False
    is_compound = False
    in_compound_block = False
    compound_block_contents = []
    re_num = re.compile(r'\d+')
    lhs = None
    for line in filehandle:
        ## When parsing compound files, the input is an array of strings, so
        ## no need to decode it twice.
        if hasattr(line, 'decode'):
            line = line.decode('utf-8','ignore')

        if not line.strip():
            continue
        if line.startswith('$$'):
            continue

        ## Detect the start of a compound block
        if is_compound and line.upper().startswith('##TITLE'):
            in_compound_block = True
            compound_block_contents = [line]
            continue

        ## If we are reading a compound block, collect lines into an array to
        ## be processed by a recursive call this this function.
        if in_compound_block:
            ## Store this line.
            compound_block_contents.append(line)

            ## Detect the end of the compound block.
            if line.upper().startswith('##END'):
                ## Process the entire block and put it into the children array.
                child_spec = spc.Spectrum()
                read(compound_block_contents, child_spec)
                jcamp_dict['children'].append( child_spec )
                in_compound_block = False
                compound_block_contents = []
            continue

        ## Lines beginning with '##' are header lines.
        if line.startswith('##'):
            line = line.strip('##')
            (lhs,rhs) = line.split('=', 1)
            lhs = lhs.strip().lower()
            rhs = rhs.strip()
            #continuation = rhs.endswith('=')

            if rhs.isdigit():
                jcamp_dict[lhs] = int(rhs)
            elif is_float(rhs):
                jcamp_dict[lhs] = float(rhs)
            # processes numbers whose decimal separator is a comma.
            elif is_float(rhs.replace(",", ".", 1)):
                jcamp_dict[lhs] = float(rhs.replace(",", ".", 1))
            else:
                jcamp_dict[lhs] = rhs

            ## Detect compound files.
            ## See table XI in http://www.jcamp-dx.org/protocols/dxir01.pdf
            if (lhs in {'data type', 'datatype'}) and (rhs.lower() == 'link'):
                is_compound = True
                jcamp_dict['children'] = []

            if lhs in ('xydata', 'xypoints', 'peak table'):
                ## This is a new data entry, reset x and y.
                x = []
                y = []
                datastart = True
                datatype = rhs

                ## Calculate x-steps from mandatory metadata. If "xfactor" is
                ## not available in jcamp_dict, then use 1.0 as default.
                if ('lastx' in jcamp_dict) and ('firstx' in jcamp_dict) and   \
                   ('npoints' in jcamp_dict):
                    dx = (jcamp_dict["lastx"] - jcamp_dict["firstx"]) /       \
                         (jcamp_dict["npoints"] - 1)
                else:
                    dx = 1.0
                dx /= jcamp_dict.get("xfactor",1)
            elif lhs == 'end':
                bounds = [int(i) for i in re_num.findall(rhs)]
                datastart = True
                datatype = bounds
                datalist = []
            elif lhs == 'longdate':
                try:
                    parsed = parse_longdate(jcamp_dict[lhs])
                except ValueError:
                    ## Keep the original date string.
                    pass
                else:
                    ## Replace the string with the datetime object.
                    jcamp_dict[lhs] = parsed
            elif datastart:
                datastart = False
        elif lhs is not None and not datastart:  # multiline entry
            jcamp_dict[lhs] += f'\n{line.strip()}'

        if datastart and (datatype == '(X++(Y..Y))'):
            ## If the line does not start with '##' or '$$' then it should be
            ## a data line. The pair of lines below involve regex splitting
            ## on floating point numbers and integers. We can't just split on
            ## spaces because JCAMP allows minus signs to replace spaces in
            ## the case of negative numbers.

            ## Check the first data line only if ASDF format is implemented.
            if len(y) > 0:
                ## Check if the format is AFFN or ASDF:
                ASDF_format_detected = any(l in DIF_digits for l in line)
            datavals = jcamp_parse(line)

            ## X-check: Is the calculated x-value the same as in first value
            ## in line? Actual implementation checks whether difference is
            ## below 1. This threshold might require adjustment to higher
            ## values if needed (not encountered so far). The line_last pair
            ## will be generated after reading first line, see code below.
            ##
            ## TODO: Tidy up to remove this error
            if "line_last" in locals():
                next_x = line_last[0] + line_last[1] * dx
                if abs(datavals[0] - next_x) > 1:
                    print("X-Check failed. Expected value is "
                          f"{datavals[0]} but {next_x} has been calculated.")

            ## Only for ASDF format: Do y-checks (to ensure line integrity) and
            ##                       do y-value aggregation appropriately
            if ASDF_format_detected:
                if len(y) > 0:
                    line_last = (datavals[0], len(datavals[2:]))
                    ## Y-check: first y-value is used to check with last
                    ##          y-value to ensure integrity of all DIF
                    ##          operations done on previous line
                    if datavals[1] != y[-1]:
                        print("Y-Check failed. Last value of previous line "
                              f"is {y[-1]} but first value is {datavals[1]}.")
                    ## Aggregate y-values.
                    for dataval in datavals[2:]:
                        y.append(float(dataval))
                else:
                    ## Aggregate y-values; first line does not contain y-checks.
                    for dataval in datavals[1:]:
                        y.append(float(dataval))
                    ## Define last x and number of y-values for next x-check.
                    line_last = (datavals[0], len(y) - 1)
            else:
                line_last = (datavals[0], len(datavals[1:]))
                for dataval in datavals[1:]:
                    y.append(float(dataval))

        elif datastart and (('xypoints' in jcamp_dict) or
                            ('xydata' in jcamp_dict)) and                     \
                            (datatype == '(XY..XY)'):
                                ## be careful not to allow empty strings
            datavals = [v.strip() for v in re.split(r"[,;\s]", line) if v]
            if not all(is_float(datavals)):
                continue
            datavals = np.array(datavals)
                                ## every other data point starting at the 0'th
            x.extend(datavals[0::2])
                                ## every other data point starting at the 1'st
            y.extend(datavals[1::2])

        elif datastart and ('peak table' in jcamp_dict) and                   \
                            (datatype == '(XY..XY)'):
                                # See lines above
            datavals = [v.strip() for v in re.split(r"[,;\s]", line) if v]
            if not all(is_float(datavals)):
                continue
            datavals = np.array(datavals)
            x.extend(datavals[0::2])
            y.extend(datavals[1::2])

        elif datastart and isinstance(datatype,list):
            ## If the line does not start with '##' or '$$' then it should be
            ## a data line. The pair of lines below involve regex splitting
            ## on floating point numbers and integers. We can't just split on
            ## spaces because JCAMP allows minus signs to replace spaces in
            ## the case of negative numbers.
            datavals = jcamp_parse(line)
            datalist += datavals

    if ('xydata' in jcamp_dict) and (jcamp_dict['xydata'] == '(X++(Y..Y))'):
        ## You got all of the Y-values. Next you need to figure out how to
        ## generate the missing X's... According to JCAMP-DX specifications,
        ## the metadata contains actual x-values. X-values in the xydata-table
        ## are used for x-checks only. The variable "xfactor" is used to
        ## compress x-values, so decompression of actual x-values is not
        ## needed anymore.
        x = np.linspace(jcamp_dict["firstx"], jcamp_dict["lastx"],
                        jcamp_dict["npoints"])
        y = np.array([float(yval) for yval in y])
    else:
        x = np.array([float(xval) for xval in x])
        y = np.array([float(yval) for yval in y])
        ## The "xfactor" variables contain any scaling information that may
        ## need to be applied to the data. Go ahead and apply them.
        if 'xfactor' in jcamp_dict:
            x = x * jcamp_dict['xfactor']

    ## Check if arrays are the same length.
    if len(x) != len(y):
        print("Mismatch of array lengths found: "
              f"len(x) is {len(x)} and len(y) {len(y)}.")

    ## The "yfactor" variables contain any scaling information that may need
    ## to be applied to the data. Go ahead and apply them.

    if 'yfactor' in jcamp_dict:
        y *= jcamp_dict['yfactor']

    my_spectrum.x_data   = x
    my_spectrum.y_data   = y
    my_spectrum.name     = jcamp_dict['title']
    # TODO do a better job with x_label
    my_spectrum.x_label  = f"Wavenumber ({jcamp_dict['xunits'].lower()})"
    my_spectrum.y_label  = f"{jcamp_dict['yunits'].capitalize()}"
    if len(jcamp_dict['children']) > 0:
        my_spectrum.metadata['Children'] = jcamp_dict['children']

##=============================================================================
def is_float(s: str) -> bool:
    '''
    Test if a string, or list of strings, contains a numeric value(s).

    Parameters
    ----------
    s : str, or list of str
        The string or list of strings to test.

    Returns
    -------
    is_float_bool : bool or list of bool
        A single boolean or list of boolean values indicating whether each
        input can be converted into a float.
    '''

    if isinstance(s,(list,tuple)):
        if not all(isinstance(i, str) for i in s):
            raise TypeError(f"Input {s} is not a list of strings")
        if len(s) == 0:
            raise ValueError(f'Input {s} is empty')
        ret_bool = list(True for i in range(0,len(s)))
        for elem, i in s:
            try:
                float(elem)
            except ValueError:
                ret_bool[i] = False
    else:
        if not isinstance(s, str):
            raise TypeError(f"Input '{s}' is not a string")
        try:
            float(s)
            ret_bool = True
        except ValueError:
            ret_bool = False
    return ret_bool

##=============================================================================

def get_value(num, is_dif, vals) -> float :
    """
    Get a value either as absolute value or difference from previous one in vals.
    """
    n = float(num)
    if is_dif:
        n += vals[-1]
    return n

##=============================================================================
def jcamp_parse(line):
    """
    Parse a line of the jcamp file
    """
    line = line.strip()

    datavals = []
    num = ""

    ## Convert whitespace into single space by splitting the string then
    ## re-assembling with single spaces.
    line = ' '.join(line.split())

    ## If there are any coded digits, then replace the codes with the
    ## appropriate numbers.
    ## 'DUP_digits': ("duplicate suppression") replaces all but first value if
    ##               two or more adjacent y-values are identical
    ## 'DIF_digits': ("difference form") replace delimiter, leading digit and
    ##               sign of the difference between adjacent values
    ## 'SQZ_digits': ("squeezed form") replace delimiter, leading digit
    ##               and sign
    dup_set = set(DUP_digits)

    if any(c in dup_set for c in line):
        ## Split the line into individual characters so that you can check for
        ## coded characters one-by-one.
        newline = ''
        for (i,c) in enumerate(line):
            if c in DUP_digits:
                ## Check for last DIF_digit which is start of last y-value by
                ## default, so that all characters belonging to last value is
                ## fully decompressed by DUP compression.
                back = 1
                while line[i-back] not in DIF_digits:
                    back += 1
                prev_c = line[i-back:i]
                mul = DUP_digits[c]
                newline += prev_c * (mul-1)
            else:
                mul = ''
                newline += c
        line = "".join(newline)

    dif = False
    for c in line:
        if c.isdigit() or (c == "."):
            num += c
        elif c == ' ':
            dif = False
            if num:
                n = get_value(num, dif, datavals)
                datavals.append(n)
            num = ''
        elif c in SQZ_digits:
            dif = False
            if num:
                n = get_value(num, dif, datavals)
                datavals.append(n)
            num = SQZ_digits[c]
        elif c in DIF_digits:
            if num:
                n = get_value(num, dif, datavals)
                datavals.append(n)
            dif = True
            num = str(DIF_digits[c])
        else:
            raise SyntaxError(f"Unknown character ({c}) encountered while "
                            "parsing data")

    if num:
        n = get_value(num, dif, datavals)
        datavals.append(n)

    return datavals

##=============================================================================
def jcamp_write(jcamp_dict: dict, linewidth=75):
    '''
    Convert a dictionary into a JDX-format string for easy writing to a file.
    At a minimum, the input dictionary must contain the keys 'x' and 'y' for
    the data, and these two vectors must have the same length.

    Parameters
    ----------
    filename : str
        The name of the JCAMP-DX file to write.
    jcamp_dict : dict
        The input dictionary to write.
    linewidth : int, optional
        The maximum line width to allow when writing data lines.

    Returns
    ----------
    jcamp_str : str
        The JCAMP_DX formatted string containing the input dictionary's
        information.
    '''

    if 'x' not in jcamp_dict :
        raise ValueError('The input dictionary *must* have an "x" variable '
                         'for writing to a JCAMP file.')
    if 'y' not in jcamp_dict :
        raise ValueError('The input dictionary *must* have a "y" variable '
                         'for writing to a JCAMP file.')

    js = ''

    ## Write the first line.
    js += "##JCAMP-DX=5.01\n"

    ## First write out the header.
    for key in jcamp_dict:
        if key in ('x','y','xydata','end'):
            continue

        js += f"##{key.upper()}={str(jcamp_dict[key])}\n"

    ## Determine whether the spectra have a title and a datetime field in the
    ## labels, by default, the title if any will be is the first string; the
    ## timestamp will be the fist datetime.datetime.
    # x = jcamp_dict['x']
    # y = jcamp_dict['y']

    if 'firstx' not in jcamp_dict:
        js += f"##FIRSTX={jcamp_dict['x'][0]:.6f}\n"
    if 'lastx' not in jcamp_dict:
        js += f"##LASTX={jcamp_dict['x'][-1]:.6f}\n"
    if 'maxx' not in jcamp_dict:
        js += f"##MAXX={np.amax(jcamp_dict['x']):.6f}\n"
    if 'minx' not in jcamp_dict:
        js += f"##MINX={np.amin(jcamp_dict['x']):.6f}\n"

    if 'firsty' not in jcamp_dict:
        js += f"##FIRSTY={jcamp_dict['y'][0]:.4f}\n"
    if 'lasty' not in jcamp_dict:
        js += f"##LASTY={jcamp_dict['y'][-1]:.4f}\n"
    if 'maxy' not in jcamp_dict:
        js += f"##MAXY={np.amax(jcamp_dict['y']):.4f}\n"
    if 'miny' not in jcamp_dict:
        js += f"##MINY={np.amin(jcamp_dict['y']):.4f}\n"

    npts = jcamp_dict.get('npts', len(jcamp_dict['x']))
    js += f"##NPOINTS={npts}\n"
    js += f"##XFACTOR={jcamp_dict.get('xfactor', 1)}\n"
    yfactor = jcamp_dict.get('yfactor', 1)
    js += f"##YFACTOR={yfactor}\n"
    js += "##XYDATA=(X++(Y..Y))\n"

    line = f"{jcamp_dict['x'][0]:.6f} "
    for j in np.arange(npts):
        if np.isnan(jcamp_dict['y'][j]):
            line += '? '
        else:
            line += f"{jcamp_dict['y'][j] / yfactor:.4f} "
#        print(npts-1, j, linewidth, len(line), f'"{line}"')
        if (len(line) >= linewidth) or (j == npts-1):
            js += line + '\n'
            if j < npts-1:
                line = f"{jcamp_dict['x'][j+1]:.6f} "

    js += '##END=\n'
    return js

## ============================================================================

def main():
    """
    A main routine to do more or less nothing.
    """
    print("This file provides routines for reading and writing jcamp files")
    return True

if __name__ == '__main__':
    main()
