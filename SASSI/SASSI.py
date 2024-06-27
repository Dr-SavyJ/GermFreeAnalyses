# -*- coding: utf-8 -*-
"""
SASSI
Signal Analysis Selection, Segmentation, and Integration

formerly SASSI
Breathing Analysis Selection and Segmentation 
for Plethysmography and Respiratory Observations
***
built as part of the Russell Ray Lab Breathing And Physiology Analysis Pipeline
***
Breathe Easy - an automated waveform analysis pipeline
Copyright (C) 2022  
Savannah Lusk, Andersen Chang, 
Avery Twitchell-Heyne, Shaun Fattig, 
Christopher Scott Ward, Russell Ray.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

***

-History-
SASSI created by Christopher S Ward (C) 2020
updates include contributions from Avery Twitchell-Heyne and feedback from
Russell Ray, Savannah Lusk, and Andersen Chang

'Breath identification and annotation functionality adapted from prior efforts 
by Christopher S Ward
    'WardBreathcaller' for python 3 (C) 2015,
    'Breath Caller' for python 2 (C) 2014
    'Breathing Analysis' for Matlab (C) 2011

***

Command Line Arguments
    '-i', help='Path containing signal files for input'
    '-f', action='append', help='names of signal files,\
            declare multiple times to build a list of files'
    '-o', help='Path to directory for output files'
    '-p', help='Path to Analysis Parameters File'
    '-a', help='Path to Animal MetaData File'
    '-m', help='Path to Manual Selection File'
        *optional, to explicitly indicate none provided, use NONE as the value
    '-c', help='Path to Automated Criteria File'
        *optional, to explicitly indicate none provided, use NONE as the value


"""


__version__ = '36.8.0'
#Update to python 3.12
#Updated HR detection function
#Updated local module import to make more dynamic whether launched via GUI or script.

"""
# v36.7.13 README

# v36.7.10 README
    *merge for tidal flow feature

# v36.7.9 README
    *update to doc string for calculate_irreg_score to reflect the *100% 
     formula update.
     
# v36.7.8 README
    *update to RUID MUID PlyUID PMID re parsing
# v36.7.7 README
    *update of Irreg score calc to add scaling factor
    *addition of I/E ratio

# v36.7.6 README
    *minor bugfix, if ecg channel present but no heart beats identified, will 
     populate an empty RR column in the basicRR function to bypass keyerror
     caused by assumption of heartbeats if an ecg channel was present

# v36.7.5 README
    *bugfixes for SASSI testing and integration for signal converters
    
# v36.7.4 README
    * updates to extract ruid muid and plyuid calls

# v36.7.3 README
    *adapted signal imprt to external modules and added channel_name remapping
    *cropped README updates to v36.4.0

# v36.7.2 README
    *added compatability for adicht signal files
    
# v36.7.1 README
    *incorporation of suggestions following code review

# v36.7.0 README
    *updated SASSI to provide multiprocessing capability and adjust logging

# v36.6.2 README
    *updated imports for compatibility with reorganized directory
    *refinements to WA lab data analysis compatibility

# v36.6.1 README
    *merged update to apnea detection rules (was v36.3.5 in python_module - 
     apparently lost as SASSI was updated seperately and missed in merge with
     autores branch)
    *revised calls to argparse args (use --flags instead of -f labels)

# v36.6.0 README
    *added compatibility for alilain lab signals 
     (flow, chamber_temp, chamber_hum)

# v36.5.1 README
    *merge with autores branch (pickle file and ruid regex bugfix)

# v36.4.1 README
    *modified arguments for main/run_SASSI

# v36.4.0 README
    *updated name of script (from python_module to SASSI)
    *parameterized main & created helper function to facilitate external calls

...
    
***

!!!
Style Suggestions

CONSTANT
name_of_function(local_variable_1, local_variable_2)
Main_Variable_1 = ...
ClassName()


!!!
Automated Selection of Breaths for Analysis - Heirarchy/Nomenclature
Timestamp_Intervals
Blocks
Sections
Breaths
Segments

"""

# %% import libraries
from tkinter import filedialog, simpledialog, Tk
import statistics
from datetime import timedelta, datetime
import csv
import re
import logging
from logging import Manager
import os
import pandas
import numpy
from scipy import signal, stats
from typing import List, Dict, Tuple, Union
import math
import pytemperature
import sys
import argparse
import multiprocessing
import signal as sig
from scipy.signal import find_peaks
import time


#%% import local modules
# Get the directory of the current script
script_directory = os.path.dirname(os.path.realpath(__file__))

# Add the script directory to sys.path if it's not already included
if script_directory not in sys.path:
    sys.path.append(script_directory)

# Attempt to import your modules, using a direct import first, 
# then falling back to a relative import if the direct import fails

try:
    from signal_converters import adi_extract
except ImportError:
    from .signal_converters import adi_extract
    
extractors = {
    'adi':{'module': adi_extract, 'ext':'.adicht'}
}
# %% define functions
def gui_open_filename(kwargs=None):
    """
    This function creates a temporary Tkinter instance that provides a GUI 
    dialog for selecting a filename.

    Parameters
    ----------
    kwargs : dict, optional
        The default is None.
        Function calls on tkinter.filedialog and uses those arguments.
        Declare as a dictionary with possible keys:
        {"defaultextension": "", "filetypes": "", "initialdir": "",
         "initialfile": "", "multiple": False, "message": "", "parent": "", "title": ""}
        
    Returns
    -------
    output_text : str
        String describing the path to the file selected by the GUI.
    """

    # It's a good practice to use `None` as a default argument instead of a mutable object like {}
    if kwargs is None:
        kwargs = {}

    # Temporarily hide the main Tkinter window
    root = Tk()
    root.withdraw()  # This line is added to hide the main window

    output_text = filedialog.askopenfilename(**kwargs)

    root.destroy()
    
    return output_text


def gui_open_filenames(kwargs=None):
    """
    This function creates a temporary Tkinter instance that provides a GUI
    dialog for selecting multiple filenames.

    Parameters
    ----------
    kwargs : dict, optional
        The default is None.
        Function calls on tkinter.filedialog and uses those arguments.
        Declare as a dictionary with possible keys:
        {"defaultextension": "", "filetypes": "", "initialdir": "",
         "initialfile": "", "multiple": False, "message": "", "parent": "", "title": ""}
        
    Returns
    -------
    output_text : list of str
        List of strings describing the paths to the files selected by the GUI.
    """

    # It's a good practice to use `None` as a default argument instead of a mutable object like {}
    if kwargs is None:
        kwargs = {}

    # Temporarily hide the main Tkinter window
    root = Tk()
    root.withdraw()  # This line is added to hide the main window

    output_text = filedialog.askopenfilenames(**kwargs)

    root.destroy()
    
    return output_text


def gui_directory(kwargs={}):
    """
    This function creates a temporary Tkinter instance that provides a GUI
    dialog for selecting a directory.

    Parameters
    ----------
    kwargs : dict, optional
        Dictionary of options for the file dialog. Common options might include
        'initialdir', 'title'. The default is {}.

    Returns
    -------
    output_text : str
        Returns the directory path selected by the GUI.
    """

    # Ensure the root window doesn't appear
    root = Tk()
    root.withdraw()

    # Open the directory selection dialog
    output_text = filedialog.askdirectory(**kwargs)

    # Properly close the Tkinter instance
    root.destroy()
    
    return output_text


def gui_get_int(title, text, default_if_canceled):
    """
    Opens a GUI dialog asking the user to enter an integer value.

    Parameters
    ----------
    title : str
        The title of the dialog window.
    text : str
        The prompt message in the dialog window.
    default_if_canceled : int
        The default value to return if the dialog is canceled or if no input is provided.

    Returns
    -------
    int
        The integer entered by the user, or the specified default value if the dialog is canceled.
    """
    
    root = Tk()
    root.withdraw()  # Hide the root window
    output_int = simpledialog.askinteger(title, text)
    root.destroy()  # Ensure the Tkinter root window is destroyed after dialog closure

    return output_int if output_int is not None else default_if_canceled


def get_avg(input_list):
    """
    Returns the average of the values in a list if it can be calculated, 
    otherwise returns 'NAN'.

    Parameters
    ----------
    input_list : list of numerical values (i.e. int, float)

    Returns
    -------
    float, string (NAN)
        arithmatic average of the values in a list 
    
    """

    try:
        return sum(input_list)/len(input_list)
    except:
        return 'NAN'


def get_med(inputlist):
    """
    Returns the median of the values in a list if it can be calculated,
    otherwise returns 'NAN'.

    Parameters
    ----------
    inputlist : list of numerical values (i.e., int, float)

    Returns
    -------
    float, int, string (NAN)
        The numeric value if median can be calculated, otherwise 'NAN'.
    """

    try:
        return statistics.median(inputlist)
    except statistics.StatisticsError:  # Catching specific exception related to calculating median
        return 'NAN'


def log_info_from_dict(local_logger, input_dict, log_prefix=''):
    """
    Creates log entries from dict input.

    Parameters
    ----------
    local_logger : instance of a logging.Logger
        Logger object used to log messages.

    input_dict : dict
        Dictionary contents will be used to populate log entries. Expected format: {key: message, ...}.

    log_prefix : str, optional
        Prefix string to prepend to each log entry. Default is an empty string.

    Returns
    -------
    None
    """

    # Iterate over dictionary items directly
    for key, message in input_dict.items():
        local_logger.info(f'{log_prefix}{key} : {message}')

def convert_seconds_to_time_text(seconds: timedelta) -> str:
    """
    Creates a string summarizing the time in hrs, mins, secs based on an input
    of duration in seconds.

    Parameters
    ----------
    seconds : datetime.timedelta
        An interval of time.

    Returns
    -------
    string
        String summarizing the time in hrs, mins, secs.
    """

    # Extract hours, minutes, and seconds directly from the timedelta object
    hrs, remainder = divmod(seconds.total_seconds(), 3600)
    mins, secs = divmod(remainder, 60)

    return f'{int(hrs)} hrs {int(mins)} mins {int(secs)} secs'

def get_animal_metadata(csvpath):
    """
    Returns a dictionary extracted from an 'Animal Metadata' csv file.
    Note: 'PlyUID' is a CASE-SENSITIVE REQUIRED column.

    Parameters
    ----------
    csvpath : str
        Path to the file containing Animal_Metadata.

    Returns
    -------
    dict
        Dictionary containing metadata associated with animals (indexed by PlyUID).
    """

    animal_metadata = {}
    if not csvpath:
        return animal_metadata

    try:
        with open(csvpath, 'r', encoding='utf-8-sig') as file:
            data = csv.DictReader(file)
            for row in data:
                plyuid = row.get('PlyUID')
                if plyuid is None:
                    continue  # Skip rows where 'PlyUID' is missing

                animal_metadata[plyuid] = {}
                for key, value in row.items():
                    try:
                        animal_metadata[plyuid][key] = float(value) if value else value
                    except ValueError:
                        animal_metadata[plyuid][key] = value
    except FileNotFoundError:
        print(f"File not found: {csvpath}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return animal_metadata

def get_animal_metadata_pneumo(csvpath):
    """
    Returns a dictionary extracted from an 'Animal Metadata' csv file.
    Note: 'RUID' is a CASE-SENSITIVE REQUIRED column.

    Parameters
    ----------
    csvpath : str
        Path to the file containing Animal_Metadata.

    Returns
    -------
    dict
        Dictionary containing metadata associated with animals (indexed by RUID).
    """

    animal_metadata = {}
    if not csvpath:  # Checks if the path is None or empty
        return animal_metadata

    try:
        with open(csvpath, 'r', encoding='utf-8-sig') as file:
            data = csv.DictReader(file)
            for row in data:
                ruid = row.get('RUID')
                if ruid is None:
                    continue  # Skip rows where 'RUID' is missing

                animal_metadata[ruid] = {}
                for key, value in row.items():
                    # Attempt to convert each value to float; fallback to original value if conversion fails
                    try:
                        animal_metadata[ruid][key] = float(value) if value else value
                    except ValueError:  # Catches conversion errors
                        animal_metadata[ruid][key] = value
    except FileNotFoundError:
        print(f"File not found: {csvpath}")
    except Exception as e:
        print(f"An error occurred while processing the file: {e}")

    return animal_metadata


def calculate_time_remaining(
        cur_file_num: int,
        tot_file_num: int,
        first_file_start: datetime,
        cur_file_start: datetime
        ) -> timedelta:
    """
    Estimates time remaining based on time used to complete analyses so far.

    Parameters
    ----------
    cur_file_num : int
        Iteration number for the current file.
    tot_file_num : int
        Total number of files submitted for analysis.
    first_file_start : datetime.datetime
        Timestamp at the start of analysis of the first file.
    cur_file_start : datetime.datetime
        Timestamp at the start of analysis for the current file.

    Returns
    -------
    datetime.timedelta
        Estimated time remaining to complete analysis of all files.
    """
    
    # Calculate the time elapsed from the start of the first file to the start of the current file
    time_elapsed = cur_file_start - first_file_start
    # Estimate the total time for all files based on the current pace
    estimated_total_time = time_elapsed / cur_file_num * tot_file_num
    # Calculate the estimated remaining time
    estimated_time_remaining = estimated_total_time - time_elapsed
    
    return estimated_time_remaining


def extract_muid_plyuid(filename, animal_metadata=None, local_logger=None):
    """
    Extracts muid and plyuid information from a filename provided in
    MUID_PLYUID.txt format.

    Parameters
    ----------
    filename : str
        Filename, expected as MUID_PLYUID.txt format.

    animal_metadata : dict, optional
        Dictionary indexed by 'PlyUID' containing animal metadata.

    local_logger : logging.Logger, optional
        Logger object for logging information about the extraction process.

    Returns
    -------
    muid : str or None
        Extracted MUID if present in the filename.

    plyuid : str or None
        Extracted PLYUID if present in the filename.
    """
    
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        local_logger.addHandler(handler)
        local_logger.setLevel(logging.INFO)

    local_logger.info(f"Extracting MUID and PLYUID from filename: {filename}")

    # Simplified regex for expected MUID_PLYUID.txt format
    #r'(?P<muid>M[^_]+)_(?P<plyuid>Ply[^.]+)\.txt$'
    muid_plyuid_re = re.compile(r'.*(?P<muid>[M][0-9]*)_(?P<plyuid>Ply[0-9]*).*.txt$', re.IGNORECASE)
    
    match = muid_plyuid_re.match(filename)
    if match:
        muid, plyuid = match.group('muid'), match.group('plyuid')
        local_logger.info(f"Successfully extracted MUID: {muid} and PLYUID: {plyuid}")
        return muid, plyuid
    else:
        local_logger.error(f"Failed to extract MUID and PLYUID from filename: {filename}")
        return None, None


def extract_ruid(filename):
    """
    Extracts ruid information from a filename provided in
    YYMMDD_RUID.txt format.

    Parameters
    ----------
    filename : str
        Filename, expected as RUID.txt format with anything surrounding the RUID.

    Returns
    -------
    ruid : str or None
        Extracted RUID if present in the filename, None otherwise.
    """
    
    # Extracts RUID from anywhere in the filename to allow for maximum flexibility in possible naming conventions.
    ruid_re = re.compile(r'.*(?P<ruid>[rR][0-9]*).*')
    match = ruid_re.match(filename)
    if match:
        return match.group(1)  # Return the matched RUID
    else:
        # Log or handle cases where the filename does not match the expected format
        return None  # Indicates that RUID was not found or filename format is incorrect


class InsufficientInformationError(Exception):
    """Custom exception for insufficient information to link MUID to PLYUID."""
    pass

class AmbiguousEntryError(Exception):
    """Custom exception for when MUID to PLYUID mapping is ambiguous."""
    pass

def resolve_plyuid(muid, animal_metadata, local_logger=None):
    """
    Identifies a plyuid using muid and information present in an 
    animal_metadata file. Custom exceptions are raised if no plyuid can be 
    found, or if a plyuid cannot be uniquely determined.

    Parameters
    ----------
    muid : str
        Serial identifier for mouse/subject.
    animal_metadata : dict
        Dictionary indexed by 'PlyUID' containing animal metadata.
    local_logger : logging.Logger, optional
        Logger object for logging information about the resolution process.

    Raises
    ------
    InsufficientInformationError
        Raised if there is no entry in metadata for the given MUID.
    AmbiguousEntryError
        Raised if MUID mapping to PLYUID is ambiguous.

    Returns
    -------
    plyuid : str
        Serial identifier for data collection event (plethysmography session ID).
    """

    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            # Setup basic logging configuration
            logging.basicConfig(level=logging.INFO)

    # Validate input parameters
    if not isinstance(animal_metadata, dict):
        raise ValueError("animal_metadata must be a dictionary")

    # Create muid to PlyUID mapping
    muid_dict = {}
    for plyuid, details in animal_metadata.items():
        muid_value = details.get('MUID')
        if muid_value:
            muid_dict.setdefault(muid_value, []).append(plyuid)

    # Check if 1 and only 1 entry can be found for the given MUID
    if muid not in muid_dict:
        msg = f'No entry for MUID {muid} in Metadata'
        local_logger.error(msg)
        raise InsufficientInformationError(msg)
    elif len(muid_dict[muid]) != 1:
        msg = f'MUID {muid} is ambiguous without PLYUID'
        local_logger.error(msg)
        raise AmbiguousEntryError(msg)

    plyuid = muid_dict[muid][0]
    return plyuid


def load_signal_data(
        filename, 
        analysis_parameters = None,
        channel_match = None,
        output_path = None,
        local_logger = None
        ):
    """
    Creates a dataframe containing plethysmography signal data
    includes calls to several other functions that are required for proper
    parsing of an exported lab chart signal file.
    -extract_header_type()
    -extract_header_locations()
    -read_exported_labchart_file()
    -merge_signal_data_pieces()
    
    *A Pickled Signal File may be used instead (improving load time). A direct
    call to pandas.read_pickle is used (file is assumed to be gzip compressed)
    columns are renamed to the traditional names provided by the legacy 
    load_signal_data function. filename must use .pkl.gzip extension
    
    *A matlab file may be used instead (but must requires inclusion of optional
    analysis_parameters to confirm style of data layout)
    !!! need info on matlab file specs from W.A. group (sunshine and taylor)

    Parameters
    ----------
    filename : string
        path to file containing signal data
    analysis_parameters : dict
        dictionary containing settings for analysis (optional)
    local_logger : instance of logging.logger (optional)


    Returns
    -------
    signal_data_assembled : pandas.DataFrame
        dataframe containing contents of signal file (merged into a single
        dataframe if multiple blocks are present)

    """

    if channel_match is None:
        if local_logger: local_logger.debug('no channel match settings')
        channel_match = pandas.DataFrame(
            {
                'setting_type':[],
                'setting_entry':[],
                'value':[]
            }
        )
    channel_match.columns = ['setting_type','setting_entry','value']
    if channel_match.empty:
        if local_logger: local_logger.info('no channel_match settings found - defaults will be used')
        channel_match_dict = {
            'time':'ts',
            'breathing':'vol',
            'oxygen':'o2',
            'co2':'co2',
            'tchamber':'temp',
            'channel 5':'ch5',
            'channel 6':'ch6',
            'channel 7':'ch7',
            'channel 8':'ch8',
            'vent flow':'flow',
            'comments':'comment',
            'all_comments':'comment'
            }
        channel_settings_dict = {}

    else:
        if local_logger: local_logger.info('loading channel_match settings')
        
        channel_match_dict = dict(
            zip(
                channel_match[channel_match['setting_type'].str.startswith('match')
                    ]['setting_entry'],
                channel_match[
                    channel_match['setting_type'].str.startswith('match')
                    ]['value']

                )            )
        channel_settings_dict = dict(
            zip(
                channel_match[
                    channel_match['setting_type'].str.startswith('setting')
                    ]['setting_entry'],
                channel_match[
                    channel_match['setting_type'].str.startswith('setting')
                    ]['value']
                )

            )
        
    
    
    filetype = []
    if channel_settings_dict.get('setting_filetype'):
        filetype.append(channel_settings_dict.get('setting_filetype'))
    for k, v in extractors.items():
        if filename.endswith(v["ext"]):
            filetype.append(k)
            
    if len(filetype) == 0:
        raise Exception('appropriate file extractor not known for signal file')
    elif len(filetype) >1 and not channel_settings_dict.get('setting_filetype'):
        if local_logger: local_logger.info(
                f'extension is ambiguous, extractor may need to be explicitly set - {filetype}'
                )
    filetype = filetype[0]
    
    if not extractors[filetype]["module"].__working__:
        raise Exception('no working file extractor available for signal file')
        

    if local_logger: local_logger.info(f'preparing to extract file using module: {filetype}')
    signal_data_assembled = extractors[filetype]['module'].SASSI_extract(
        filename,
        logger = local_logger
        )
    signal_data_assembled = signal_data_assembled.rename(columns = channel_match_dict)
    if local_logger: local_logger.info('signal data loaded')
    
    if channel_match.set_index('setting_type').to_dict('index').get(
            'setting_export_pklgzip',{'value':'false'}
    )['value'].lower() == 'true' and filetype != 'pklgzip':
        if local_logger: local_logger.info('exporting pkl.gzip file')
        extractors[filetype]['module'].convert_to_pickle(filename,output_path)
    
    return signal_data_assembled

def update_channel_match_dicts(channel_match, channel_match_dict, channel_settings_dict, logger):
    """
    Update channel match and settings dictionaries based on user-specified DataFrame.

    Parameters
    ----------
    channel_match : pd.DataFrame
        DataFrame with user-specified channel match and settings.
    channel_match_dict : dict
        Dictionary to update with channel matches.
    channel_settings_dict : dict
        Dictionary to update with channel settings.
    logger : logging.Logger
        Logger for logging the process.
    """
    logger.info("Updating channel match and settings dictionaries with user-specified values.")
    for _, row in channel_match.iterrows():
        if row['setting_type'].startswith('match'):
            channel_match_dict[row['setting_entry']] = row['value']
        elif row['setting_type'].startswith('setting'):
            channel_settings_dict[row['setting_entry']] = row['value']


def extract_header_type(filename,rows_to_check = 20):
    """
    Gathers information regarding the header format/column contents present
    in an exported lab chart signal file (assumes set-up is in line with
    Ray Lab specifications)

    Parameters
    ----------
    filename : string
        path for file containing signal data

    Returns
    -------
    header_tuples : list of tuples
        list of tuples specifying ([column name],[datatype])

    """
    
    # Initialize default columns
    ts_columns = ['ts']
    header_columns = []
    
    with open(filename, 'r') as opfi:
        for _ in range(rows_to_check):
            cur_line = opfi.readline()
            if "DateFormat=	M/d/yyyy" in cur_line:
                ts_columns.append('date')
            elif "ChannelTitle=" in cur_line:
                header_columns = cur_line.lower().strip().split('\t')[1:]
                header_columns = [col for col in header_columns if col]

    # Define special and rename columns
    special_columns = {'date': str, 'comment': str}
    rename_columns = {
        'breathing': 'vol',
        'oxygen': 'o2',
        'oxygen ': 'o2',  # Duplicate key with trailing space?
        'co2': 'co2',
        'tchamber': 'temp',
        'channel 5': 'ch5',
        'channel 6': 'ch6',
        'channel 7': 'ch7',
        'channel 8': 'ch8',
        'vent flow': 'flow'
    }

    # Combine and rename header columns as necessary
    combined_columns = ts_columns + [rename_columns.get(col, col) for col in header_columns] + ['comment']

    # Generate header tuples with data types
    header_tuples = [(col.lower(), special_columns.get(col, float)) for col in combined_columns]

    return header_tuples


def extract_header_locations(
        filename,
        header_text_fragment='Interval=',
        local_logger=None):
    """
    Gathers information regarding the locations of header information
    throughout a signal file - needed if files may contain multiple recording
    blocks

    Parameters
    ----------
    filename : string
        path to file containing signal data
    header_text_fragment : string, optional
        string that is present in header lines. The default is 'Interval='.
    local_logger : instance of logging.logger, optional
        The default is None (i.e. no logging)

    Returns
    -------
    headers : list
        list of rows in the datafile that indicate header content present
        
    """
    
    # Ensure that a logger is available
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            local_logger.addHandler(handler)
            local_logger.setLevel(logging.INFO)

    headers = []
    with open(filename, 'r') as opfi:
        for i, line in enumerate(opfi):
            if header_text_fragment in line:
                local_logger.info(f'Signal File has HEADER AT LINE: {i}')
                headers.append(i)

    return headers


def read_exported_labchart_file(
        lc_filepath,
        header_locations,
        header_tuples,
        delim='\t',
        rows_to_skip=6
        ):
    """
    Collects data from an exported lab chart file and returns a list of
    dataframes (in order) containing the extracted contents of the signal file.

    Parameters
    ----------
    lc_filepath : string
        path to file containing signal data
    header_locations : list of integers
        list describing the locations of headers throughout a signal file
    header_tuples : list of tuples
        list of tuples specifying ([column name],[datatype])
    delim : string, optional
        delimiter used. The default is '\t'.
    rows_to_skip : integer, optional
        the number of rows present in the header that should be skipped
        to get to the location containing data. The default is 6.

    Returns
    -------
    df_list : list of pandas.DataFrames
        list of dataframes containing signal data

    """
    df_list = []
    column_names = [name for name, _ in header_tuples]
    column_types = dict(header_tuples)

    for idx, location in enumerate(header_locations):
        # Determine the number of rows to read for this section
        nrows = None  # Assume reading till the end if it's the last section
        if idx + 1 < len(header_locations):
            nrows = header_locations[idx + 1] - location - rows_to_skip

        # Adjust skiprows to account for header location and rows to skip
        skiprows = rows_to_skip + location

        # Read the section of the file as a DataFrame
        df = pandas.read_csv(
            lc_filepath,
            sep=delim,
            names=column_names,
            skiprows=skiprows,
            nrows=nrows,
            dtype=column_types
        )
        df_list.append(df)

    return df_list


def merge_signal_data_pieces(df_list):
    """
    Merges multiple blocks contained in a list of dataframes into a single
    dataframe (current behavior will override timestamp information to
    place subsequent blocks of data using the next sequential timestamp).

    Parameters
    ----------
    df_list : list of pandas.DataFrames
        list of dataframes containing signal data

    Returns
    -------
    merged_data : pandas.DataFrame
        dataframe containing signal data

    """
    if not df_list:  # If df_list is empty, return an empty DataFrame
        return pandas.DataFrame()

    # Initialize merged_data with the first piece to avoid an empty DataFrame in the loop
    df_list[0]['ts'] = df_list[0]['ts'] - df_list[0]['ts'].min()  # Normalize the first block's timestamps if not starting from 0
    merged_data = df_list[0].copy()

    # Loop through subsequent DataFrames and adjust their timestamps
    for i in range(1, len(df_list)):
        # Calculate timestamp adjustments
        ts_adjustment = merged_data['ts'].max() + (df_list[0]['ts'][1] - df_list[0]['ts'][0])

        # Adjust and append each subsequent DataFrame
        df_list[i]['ts'] = df_list[i]['ts'] - df_list[i]['ts'].min() + ts_adjustment
        merged_data = pandas.concat([merged_data, df_list[i]], ignore_index=True)

    return merged_data


def apply_voltage_corrections(
        input_df, 
        factor_dict, 
        column_tuples,
        local_logger=None
        ):
    """
    Applies simple scaling factors to series within a dataframe specified by 
    settings in the factor_dict and column_tuples.

    Parameters
    ----------
    input_df : pandas.DataFrame
        DataFrame to modify by appending new columns.
    factor_dict : dict
        Simple dict containing factors to use for modification.
    column_tuples : list of tuples
        Each tuple containing 3 strings: (source column in input_df, destination column in input_df, key in factor_dict for factor).
    local_logger : logging.Logger, optional
        Logger for logging messages. Default is None (i.e., no logging).

    Returns
    -------
    df : pandas.DataFrame
        DataFrame that can be used to replace input_df, with voltage corrections applied.
    """
    # Ensure that a logger is available
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            local_logger.addHandler(handler)
            local_logger.setLevel(logging.INFO)

    # Copy the input DataFrame to avoid modifying the original data
    df = input_df.copy()

    # Apply scaling factors as specified in column_tuples
    for source_column, dest_column, factor_key in column_tuples:
        if source_column in df.columns:
            try:
                df[dest_column] = df[source_column] * float(factor_dict.get(factor_key, 1))  # Default to 1 if factor_key is not found
            except KeyError:
                local_logger.warning(f'Factor for {factor_key} not found in factor_dict. No correction applied to {source_column}.')
        else:
            local_logger.warning(f'{source_column} not present in DataFrame - no voltage correction applied.')

    return df     
    

def check_animal_metadata_and_analysis_parameters(
        animal_metadata,
        plyuid,
        analysis_parameters,
        parameter,
        default,
        local_logger=None):
    """
    Checks the contents of animal_metadata and analysis_parameters to determine
    if an analysis_parameter is 'overridden' by a setting in the 
    animal_metadata. The value of the setting is returned. Priority is to use 
    the value from animal_metadata, otherwise the value from 
    analysis_parameters, and finally a default value if no value is present in
    animal_metadata nor analysis_parameters.

    Parameters
    ----------
    
    animal_metadata : dict
        dict indexed by 'PlyUID' containing animal metadata
    plyuid : string
        unique identifier for data collection event 
        (plethysmography session id)
    analysis_parameters : dict
        dictionary containing settings for analysis
    parameter : string
        name of the parameter being checked
    default : string
        default value for parameter
    local_logger : instance of logging.logger (optional)
    
    Returns
    -------
    checked_value : string
        The determined value for the parameter.

    """
    # Ensure that a logger is available
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            local_logger.addHandler(handler)
            local_logger.setLevel(logging.INFO)

    # Check animal_metadata first
    if plyuid in animal_metadata and parameter in animal_metadata[plyuid] and animal_metadata[plyuid][parameter]:
        checked_value = animal_metadata[plyuid][parameter]
        local_logger.info(f"Using animal-specific setting for '{parameter}' from metadata for PLYUID {plyuid}.")
    # Then check analysis_parameters
    elif parameter in analysis_parameters and analysis_parameters[parameter]:
        checked_value = analysis_parameters[parameter]
        local_logger.info(f"Using analysis parameter for '{parameter}'.")
    # Default to the provided default value
    else:
        checked_value = default
        local_logger.info(f"Using default value for '{parameter}': {default}.")

    return checked_value


def make_reverse_timestamp_dict(timestamp_dict, key_list=None):
    """
    Creates a reverse dictionary from a given timestamp dictionary, optionally filtering
    by a specified list of keys.

    Parameters
    ----------
    timestamp_dict : dict
        The original dictionary to reverse, with timestamp values that may start with '#* '.
    key_list : list, optional
        A list of keys to include in the reversed dictionary. If None, include all keys.

    Returns
    -------
    reverse_timestamp_dict : dict
        The reversed dictionary, with values as keys and keys as values.
    """
    reverse_timestamp_dict = {}
    
    # Create the reverse dictionary, processing values as described
    for k, v in timestamp_dict.items():
        new_key = v[3:].strip() if v.startswith('#* ') else v.strip()
        reverse_timestamp_dict[new_key] = k
    
    # If a key_list is provided, filter the dictionary to include only those keys
    if key_list is not None:
        reverse_timestamp_dict = {k: v for k, v in reverse_timestamp_dict.items() if k in key_list}

    return reverse_timestamp_dict


def repair_temperature(
        signal_data,
        plyuid,
        animal_metadata,
        analysis_parameters,
        local_logger=None
        ):
    """
    temperature signal data is checked for anomolous values. 
    
    If no chamber temperature is present then a default value from 
    analysis_parameters or animal_metadata is used.
    If chamber temperature is present the first check uses cutoff values to 
    identify out of range values. Outlier regions are replaced with imputed 
    values. An additional optional check repeats this process where outliers 
    are identified as being greater than the IQR (75% - 25%) away from the 
    median temperature. This function utilizes 
    filter_and_replace_temperature() to apply corrections.

    Parameters
    ----------
    signal_data : pandas.DataFrame
        dataframe containing signal data (including corrected_temp column)
    plyuid : string
        serial identifier for plethysmography recording session
    animal_metadata : dict
        dict indexed by 'PlyUID' containing animal metadata
    analysis_parameters : dict
        dictionary containing settings for analysis
    local_logger : instance of logging.logger (optional)    
    
    **Note - column in animal_metadata or analysis_parameters contain values 
    used by this function
        * 'corrected_temp_default'
        * 'chamber_temp_cutoffs'
        * 'chamber_temperature_units'
        * 'chamber_temperature_default'
        * 'chamber_temperature_trim_size'
        * 'chamber_temperature_narrow_fix'

    Returns
    -------
    pandas.Series
        series containing corrected temperature values

    """
    # Ensure a logger is available
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            local_logger.addHandler(handler)
            local_logger.setLevel(logging.INFO)

    df = signal_data.copy()
    if 'corrected_temp' not in df.columns:
        local_logger.warning('corrected_temp not in signal data, replacing with default')
        default_temp = float(check_animal_metadata_and_analysis_parameters(
            animal_metadata, plyuid, analysis_parameters, 'corrected_temp_default', '35', local_logger))
        df['corrected_temp'] = default_temp
    else:
        # Perform various checks and corrections on the temperature data
        cutoffs_raw = check_animal_metadata_and_analysis_parameters(
            animal_metadata, plyuid, analysis_parameters, 'chamber_temp_cutoffs', '18 35', local_logger)
        cutoffs = [float(c) for c in cutoffs_raw.split()]

        temp_units = check_animal_metadata_and_analysis_parameters(
            animal_metadata, plyuid, analysis_parameters, 'chamber_temperature_units', 'C', local_logger)

        default_temp = float(check_animal_metadata_and_analysis_parameters(
            animal_metadata, plyuid, analysis_parameters, 'chamber_temperature_default', '28', local_logger))

        temperature_trim_size = int(check_animal_metadata_and_analysis_parameters(
            animal_metadata, plyuid, analysis_parameters, 'chamber_temperature_trim_size', '1000', local_logger))

        temperature_narrow_fix = check_animal_metadata_and_analysis_parameters(
            animal_metadata, plyuid, analysis_parameters, 'chamber_temperature_narrow_fix', 'T', local_logger)

        if temp_units.upper() == 'F':
            local_logger.info('Converting chamber temperature from F to C')
            df['corrected_temp'] = (df['corrected_temp'] - 32) * (5 / 9)

        # Apply the initial filtering and replacing of temperature
        df['corrected_temp'] = filter_and_replace_temperature(
            df['corrected_temp'], cutoffs, default_temp, temperature_trim_size, local_logger)

        if temperature_narrow_fix.upper() == 'T':
            local_logger.info('Applying enhanced chamber temperature check')
            # This is a simplified approach to reapply filtering based on a narrowed IQR criterion
            median_temp = numpy.median(df['corrected_temp'])
            iqr = numpy.subtract(*numpy.percentile(df['corrected_temp'], [75, 25]))
            narrowed_cutoffs = [median_temp - iqr, median_temp + iqr]

            df['corrected_temp'] = filter_and_replace_temperature(
                df['corrected_temp'], narrowed_cutoffs, default_temp, temperature_trim_size, local_logger)

    return df['corrected_temp']


def filter_and_replace_temperature(
        input_series,
        cutoffs=(10, 100),
        default=26,
        trim_size=1000,
        local_logger=None
        ):
    """
    Checks a series of data for values outside of cutoff values (exclusive), 
    and trims back neighboring values. The gap is then filled in with linear
    interpolated values.

    Parameters
    ----------
    input_series : DataSeries of Floats
        list of temperature values being screened for filtering/replacement
    cutoffs : tuple of Floats
        Boundary temperatures to use, values outside of this range are 
        identified and replaced
    default : Float, optional
        Temperature Value to use if temperature cannot be imputed from 
        neighboring measurements. The default is 26.
    trim_size : int, optional
        Number of samples to trim the data from the regio with the out of range
        value. The default is 1000.
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    temperature_series : DataSeries of Floats
        list of 'repaired' temperatures

    """
    # Ensure that a logger is available
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            local_logger.addHandler(handler)
            local_logger.setLevel(logging.INFO)

    # Copy the input Series to avoid modifying the original data
    temperature_series = input_series.copy()

    # Identify out-of-range values
    out_of_range = (temperature_series < cutoffs[0]) | (temperature_series > cutoffs[1])
    
    # Expand the filter to include neighboring values based on trim_size
    expanded_filter = out_of_range.rolling(window=trim_size*2, center=True, min_periods=1).max().fillna(1)

    # Replace out-of-range values with NaN
    temperature_series[expanded_filter == 1] = numpy.nan
    
    # Interpolate and fill remaining gaps with default value
    temperature_series.interpolate(method='linear', inplace=True)
    temperature_series.fillna(method='backfill', inplace=True)
    temperature_series.fillna(method='ffill', inplace=True)
    temperature_series.fillna(default, inplace=True)

    # Log changes
    if not out_of_range.any():
        local_logger.info('No changes made to chamber temperature data.')
    else:
        local_logger.warning('Chamber temperature contains out-of-range data; fixes applied.')

    return temperature_series


def calculate_body_temperature(
        signal_data,
        plyuid,
        timestamp_dict,
        animal_metadata,
        manual_selection,
        auto_criteria,
        analysis_parameters,
        local_logger
        ):
    """
    Calculates a linear interpolation of body temperature based on temperatures
    measuered at the start, middle, and end of the experiment. This funciton 
    utilizes extract_body_temperature_from_metadata(), and 
    extract_body_temperature_timings() to collect temperature and time 
    information.

    Parameters
    ----------
    signal_data : Pandas.DataFrame
        DataFrame containing signal data
    plyuid : string
        serial identifier for data collection event 
    timestamp_dict : Dict
        Dictionary containing timestamps and text used to describe them. 
        Captured from the commends in the signal_data
    animal_metadata : dict
        dict indexed by 'PlyUID' containing animal metadata
    manual_selection : Dict
        Dict with info describing manually selected start and stop points for
        breathing of interest
    auto_criteria : Dict
        Dict with info describing criteria for automated selection of breathing
        for 'calm quiet' breathing
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    DataSeries of Floats
        Series of Floats describing body temperature

    """
    # Assume these functions are implemented and correctly return temperatures and timestamps
    start_body_temp, mid_body_temp, end_body_temp = extract_body_temperature_from_metadata(animal_metadata, plyuid, local_logger)
    startpoint, midpoint, endpoint = extract_body_temperature_timings(plyuid, manual_selection, auto_criteria, analysis_parameters, timestamp_dict, local_logger)

    # Validate or set default timestamps
    startpoint = startpoint if startpoint is not None else signal_data['ts'].min()
    endpoint = endpoint if endpoint is not None else signal_data['ts'].max()
    midpoint = midpoint if midpoint is not None else (startpoint + endpoint) / 2  # Correct calculation of midpoint

    # Prepare DataFrame for interpolation
    df = pandas.DataFrame({
        'ts': signal_data['ts'],
        'body_temperature': numpy.nan  # Initialize with NaN for interpolation
    })

    # Time segments for interpolation
    time_segments = [(startpoint, midpoint, start_body_temp, mid_body_temp), (midpoint, endpoint, mid_body_temp, end_body_temp)]

    # Apply linear interpolation for each segment
    for seg_start, seg_end, temp_start, temp_end in time_segments:
        segment_mask = (df['ts'] >= seg_start) & (df['ts'] <= seg_end)
        segment_length = seg_end - seg_start
        segment_ts_normalized = (df.loc[segment_mask, 'ts'] - seg_start) / segment_length
        df.loc[segment_mask, 'body_temperature'] = temp_start + (temp_end - temp_start) * segment_ts_normalized

    # Fill in the end temperature beyond the endpoint if needed
    df.loc[df['ts'] > endpoint, 'body_temperature'] = end_body_temp

    # Log the operation
    local_logger.info('Body temperature interpolated based on start, middle, and end measurements.')

    return df['body_temperature']


def extract_body_temperature_from_metadata(
        animal_metadata,
        plyuid,
        local_logger,
        default=37.5):
    """
    Extracts temperature data from columns in animal_metadata. If no data is
    present a default value is used.
    Columns in animal_metadata containing temperature data:
    * 'Start_body_temperature'
    * 'Mid_body_temperature'
    * 'End_body_temperature'

    Parameters
    ----------
    animal_metadata : dict
        dict indexed by 'PlyUID' containing animal metadata
    plyuid : string
        serial identifier for data collection event 
        (plethysmography session id)
    local_logger : instance of logging.logger (optional)
    default : Float, optional
        Default body temperature to use if no data available. 
        The default is 37.5.

    Returns
    -------
    start_body_temp : Float
        Body Temperature of the animal at the beginning of the recording
    mid_body_temp : Float
        Body Temperature of the animal at the 'middle' of the recording
    end_body_temp : Fload
        Body Temperature of the animal at the end of the recording

    """
    # Ensure that a logger is available
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            local_logger.addHandler(handler)
            local_logger.setLevel(logging.INFO)

    # Initialize temperatures with None
    start_body_temp = mid_body_temp = end_body_temp = None

    # Attempt to extract temperature data
    try:
        start_body_temp = float(animal_metadata[plyuid].get('Start_body_temperature', default))
    except KeyError:
        local_logger.warning('No starting body temperature present in metadata for PLYUID "{}"'.format(plyuid))
    
    try:
        mid_body_temp = float(animal_metadata[plyuid].get('Mid_body_temperature', default))
    except KeyError:
        local_logger.warning('No midpoint body temperature present in metadata for PLYUID "{}"'.format(plyuid))
    
    try:
        end_body_temp = float(animal_metadata[plyuid].get('End_body_temperature', default))
    except KeyError:
        local_logger.warning('No ending body temperature present in metadata for PLYUID "{}"'.format(plyuid))

    # Handle cases where any of the temperatures are None
    temps = [start_body_temp, mid_body_temp, end_body_temp]
    if all(temp is None for temp in temps):
        local_logger.warning(f'No body temperature data present for PLYUID "{plyuid}", default value ({default}) used.')
        start_body_temp = mid_body_temp = end_body_temp = default
    else:
        # Fill in missing temperatures
        if start_body_temp is None:
            start_body_temp = mid_body_temp if mid_body_temp is not None else end_body_temp
            local_logger.info('Filling missing starting body temperature.')

        if mid_body_temp is None:
            mid_body_temp = (start_body_temp + (end_body_temp if end_body_temp is not None else start_body_temp)) / 2
            local_logger.info('Filling missing midpoint body temperature.')

        if end_body_temp is None:
            end_body_temp = mid_body_temp if mid_body_temp is not None else start_body_temp
            local_logger.info('Filling missing ending body temperature.')

    return start_body_temp, mid_body_temp, end_body_temp


def extract_body_temperature_timings(
        plyuid,
        manual_selection,
        auto_criteria,
        analysis_parameters,
        timestamp_dict,
        local_logger
        ):
    """
    Extracts timing information of body temperature measurements for use in 
    creating a linear interpolation of the body temperature. Timepoints are
    pulled from manual_selection or auto_criteria. Note, timepoints from 
    auto_criteria are based on the placement of the 'key' timestamp, without 
    any offsets or auto_criteria filters applied.

    Parameters
    ----------
    plyuid : string
        serial identifier for data collection event
    manual_selection : Dict
        Dict with info describing manually selected start and stop points for
        breathing of interest
    auto_criteria : Dict
        Dict with info describing criteria for automated selection of breathing
        for 'calm quiet' breathing
    timestamp_dict : Dict
        Dictionary containing timestamps and text used to describe them. 
        Captured from the commends in the signal_data
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    startpoint : float of None
        timestamp key associated startpoint in autoc_criteria or start time
        associated with startpoint in manual_selection
    midpoint : float or None
        timestamp key associated midpoint in autoc_criteria or start time
        associated with midpoint in manual_selection
    endpoint : float or None
        timestamp key associated endpoint in autoc_criteria or stop time
        associated with endpoint in manual_selection

    """
    # Initialize default return values
    startpoint = midpoint = endpoint = None

    # Extract auto criteria timings if available
    if analysis_parameters.get('Pneumo_Mode') == '1':
        if auto_criteria:
            for key, value in timestamp_dict.items():
                label = value.strip()[3:]  # Assuming value format starts with "#* "
                if label == auto_criteria.get('Start_key'):
                    startpoint = key
                    local_logger.info(f'startpoint found in auto_criteria for PLYUID {plyuid}')
                elif label == auto_criteria.get('Mid_key'):
                    midpoint = key
                    local_logger.info(f'midpoint found in auto_criteria for PLYUID {plyuid}')
                elif label == auto_criteria.get('End_key'):
                    endpoint = key
                    local_logger.info(f'endpoint found in auto_criteria for PLYUID {plyuid}')
    else: 
        if auto_criteria.shape[0]>0:
            for key, value in timestamp_dict.items():
                label = value.strip()[3:]  # Assuming value format starts with "#* "
                if label == auto_criteria.get('Start_key'):
                    startpoint = key
                    local_logger.info(f'startpoint found in auto_criteria for PLYUID {plyuid}')
                elif label == auto_criteria.get('Mid_key'):
                    midpoint = key
                    local_logger.info(f'midpoint found in auto_criteria for PLYUID {plyuid}')
                elif label == auto_criteria.get('End_key'):
                    endpoint = key
                    local_logger.info(f'endpoint found in auto_criteria for PLYUID {plyuid}')

    # Override with manual selections if available
    if manual_selection:
        startpoint = manual_selection.get('Startpoint', startpoint)
        midpoint = manual_selection.get('Midpoint', midpoint)
        endpoint = manual_selection.get('Endpoint', endpoint)

    # Log if manual selections are used
    if startpoint and startpoint == manual_selection.get('Startpoint'):
        local_logger.info(f'Using manual startpoint for PLYUID {plyuid}')
    if midpoint and midpoint == manual_selection.get('Midpoint'):
        local_logger.info(f'Using manual midpoint for PLYUID {plyuid}')
    if endpoint and endpoint == manual_selection.get('Endpoint'):
        local_logger.info(f'Using manual endpoint for PLYUID {plyuid}')

    return startpoint, midpoint, endpoint


def apply_smoothing_filter(signal_data,
                           column,
                           high_pass=0.1,
                           high_pass_order=2,
                           low_pass=50,
                           low_pass_order=10,
                           ):
    """
    Applies a highpass and lowpass filter to data, useful for reducing 
    artifacts from electrical noise or bias flow pump vibrations. It is 
    intended for use on a plethysmography or pneumotacography 'flow' signal.
    
    Parameters
    ----------
    signal_data : pandas.DataFrame
        data to be smoothed
    column : 'String'
        The name of the 
    high_pass : Float, optional
        Frequency cutoff (Hz) for the highpass filter. The default is 0.1.
    high_pass_order : Integer, optional
        order value for the high pass filter. The default is 2.
    low_pass : Float, optional
        Frequency cutoff (Hz) for the low_pass filter. The default is 50.
    low_pass_order : Integer, optional
        order value for the low pass filter. The default is 10.
     
    Returns
    -------
    lpf_hpf_signal : List, DataSeries of Floats
        smoothed data as a numpy array

    """
    # Ensure the timestamp column 'ts' is sorted
    signal_data.sort_values(by='ts', inplace=True)

    # Calculate sampling frequency
    sample_interval = signal_data['ts'].diff().median()  # More robust than subtracting two points
    sampleHz = round(1 / sample_interval)

    # High-pass filter
    hpf_b, hpf_a = signal.butter(high_pass_order, high_pass / (sampleHz / 2), 'high')
    hpf_signal = signal.filtfilt(hpf_b, hpf_a, signal_data[column])

    # Low-pass filter
    lpf_b, lpf_a = signal.bessel(low_pass_order, low_pass / (sampleHz / 2), 'low')
    lpf_hpf_signal = signal.filtfilt(lpf_b, lpf_a, hpf_signal)

    return lpf_hpf_signal


def basicFilt(CT,sampleHz,f0,Q):
    """
    Applies a notch and butter filter to data, useful for reducing artifacts 
    from electrical noise or voltage offsets. It is intended for use on an 
    ECG signal.

    Parameters
    ----------
    CT : list or pandas.Series
        ecg voltage values
    sampleHz : Float
        the sampling rate of the data
    f0 : Float
        The target frequency to exclude
    Q : Float
        Quality factor. Dimensionless parameter that characterizes
        notch filter -3 dB bandwidth ``bw`` relative to its center
        frequency, ``Q = w0/bw``.

    Returns
    -------
    filtered : list or pandas.Series
        the filtered data

    """
    # Notch filter to remove the f0 frequency component
    b, a = signal.iirnotch(f0 / (sampleHz / 2), Q)
    notched = signal.filtfilt(b, a, CT)
    
    # High-pass Butterworth filter to remove low-frequency noise
    b, a = signal.butter(1, 1 / (sampleHz / 2), btype='highpass')
    filtered = signal.filtfilt(b, a, notched)
    
    # Return the same type as the input
    if isinstance(CT, pandas.Series):
        return pandas.Series(filtered, index=CT.index)
    else:
        return filtered


def smooth_gas_signals(signal_data, analysis_parameters, samplingHz=1000, window_sec=10, logger=None):
    """
    Optionally applies median smoothing to gas signal data ('corrected_o2' and 'corrected_co2')
    based on analysis parameters.

    Parameters
    ----------
    signal_data : pandas.DataFrame
        DataFrame containing gas signal data.
    analysis_parameters : dict
        Dictionary containing settings for analysis.
    samplingHz : int, optional
        Sampling frequency in Hz. The default is 1000.
    window_sec : int, optional
        Window size for the rolling median in seconds. The default is 10.
    logger : logging.Logger, optional
        Logger for logging messages. If None, a default logger is used.

    Returns
    -------
    tuple of pandas.Series
        Smoothed (or original, if smoothing not applied) O2 and CO2 signals.
    """
    # Ensure a logger is available
    if logger is None:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)

    # Check for 'g' in the analysis_parameters to apply smoothing
    if 'g' in str(analysis_parameters.get('apply_smoothing_filter', '')):
        try:
            smoothed_o2 = signal_data['corrected_o2'].rolling(int(samplingHz * window_sec), min_periods=1).median()
            smoothed_co2 = signal_data['corrected_co2'].rolling(int(samplingHz * window_sec), min_periods=1).median()
            logger.info("Applied smoothing to 'corrected_o2' and 'corrected_co2'.")
        except KeyError as e:
            logger.error(f"Missing column in signal_data: {e}.")
            raise
    else:
        smoothed_o2 = signal_data['corrected_o2']
        smoothed_co2 = signal_data['corrected_co2']
        logger.info("Smoothing not applied. Returning original signals.")

    return smoothed_o2, smoothed_co2
            

def getmaxbyindex(inputlist, indexlist):
    """
    Returns a list (numpy array) of maximum values between specified indexes.
    Bounds are set as [index[n] to index[n+1]).

    This function can be used for calculations like Peak Inspiratory Flow (PIF)
    using flow data and transition times. Note that a cross-reference of
    indexlist with Ti vs Te timestamps is needed for segregation of the data.

    Parameters
    ----------
    inputlist : list or numpy.array
        List or array to extract maximum values from.
    indexlist : list of int
        List of index values specifying the start of the unit to screen within,
        where n+1 specifies the end of the unit.

    Returns
    -------
    maxbyindexlist : numpy.array
        Array of the maximum values specified by the index boundaries.
    """
    maxbyindexlist = [max(inputlist[indexlist[i]:indexlist[i+1]])
                      for i in range(len(indexlist)-1)]
    return numpy.array(maxbyindexlist)


def getminbyindex(inputlist, indexlist):
    """
    Returns a numpy array of minimum values between specified indexes.
    Bounds are set as [index[n] to index[n+1]).
    
    This function can be used for calculations like Peak Expiratory Flow (PEF)
    using flow data and transition times. Note that cross-reference of
    indexlist with Ti vs Te timestamps is needed for segregation of the data.
    
    Parameters
    ----------
    inputlist : list or numpy.array
        List or array to extract minimum values from.
    indexlist : list of int
        List of index values specifying the start of the unit to screen within, 
        where n+1 specifies the end of the unit.

    Returns
    -------
    minbyindexlist : numpy.array
        Array of the minimum values specified by the index boundaries.
    """
    minbyindexlist = [min(inputlist[indexlist[i]:indexlist[i+1]])
                      for i in range(len(indexlist)-1)]
    return numpy.array(minbyindexlist)


def getsumby2index(inputlist, index1, index2):
    """
    Returns a list of sums between specified indexes, with each pair of
    indexes defining a range for summation. This can be useful for
    calculating total volumes (TV) using flow data and transition times,
    requiring segregation of the data based on time indexes.

    Parameters
    ----------
    inputlist : list of float or int
        List from which to extract sums.
    index1 : list of int
        List of lower bounds for ranges of entries to sum.
    index2 : list of int
        List of upper bounds for ranges of entries to sum.

    Returns
    -------
    sumbyindexlist : list of float
        List of sums for the values specified by the index boundaries.
    """
    # Ensure index1 and index2 have the same length to match pairs correctly
    if len(index1) != len(index2):
        raise ValueError("index1 and index2 must have the same length.")
    
    # Calculate sums for each range defined by corresponding elements in index1 and index2
    sumbyindexlist = [sum(inputlist[start:end]) for start, end in zip(index1, index2)]
    
    return sumbyindexlist


def getmaxby2index(inputlist: List[float], index1: List[int], index2: List[int]) -> List[float]:
    """
    Returns a list of maximum values between specified indexes.
    Bounds are set as [index1[n] to index2[n]).

    Parameters
    ----------
    inputlist : List[float]
        List from which to extract maximum values.
    index1 : List[int]
        List of lower bounds of indexes for ranges of entries to find maximums.
    index2 : List[int]
        List of upper bounds of indexes for ranges of entries to find maximums.

    Returns
    -------
    maxbyindexlist : List[float]
        List of maximum values specified by the index boundaries.
    """

    maxbyindexlist = [max(inputlist[index1[i]:index2[i]])
                      for i in range(len(index1))]
    return maxbyindexlist


def get_avg_by_2_index(inputlist, index1, index2):
    """
    Returns a list of average values between specified indexes.
    Bounds are set as [index1[n] to index2[n]).
    
    Parameters
    ----------
    inputlist : list
        List to extract averages from.
    index1 : list of int
        List of lower bounds of indexes for ranges of entries to average.
    index2 : list of int
        List of upper bounds of indexes for ranges of entries to average.

    Returns
    -------
    avgbyindexlist : list of float
        List of averages for the values specified by the index boundaries.
    """
    
    avgbyindexlist = [numpy.mean(inputlist[index1[i]:index2[i]])
                      for i in range(len(index1))]
    return avgbyindexlist


def extract_candidate_breath_transitions(signal_data):
    """
    Generates a DataFrame describing candidate breaths using signal above and
    below zero to define transitions between inhalation and exhalation, based
    on flow signal data.

    Parameters
    ----------
    signal_data : pd.DataFrame
        DataFrame containing signal data with 'flow' and 'ts' (timestamp) columns.

    Returns
    -------
    candidate_breath_transitions : pd.DataFrame
        DataFrame with columns describing candidate breath parameters: timestamps,
        inhale vs exhale, PIF, PEF, and index list.
    """
    # Ensure 'flow' and 'ts' columns exist in signal_data
    if 'flow' not in signal_data.columns or 'ts' not in signal_data.columns:
        raise ValueError("signal_data must contain 'flow' and 'ts' columns.")
    
    # Identifying transitions
    signal_above_zero = numpy.greater(signal_data['flow'], 0)
    dif_signal_above_zero = numpy.diff(signal_above_zero.astype(int))
    index_list = numpy.array(dif_signal_above_zero.nonzero())[0]

    # Preparing the output DataFrame
    candidate_breath_transitions = pandas.DataFrame({
        'ts': numpy.take(signal_data['ts'], index_list,axis=0),
        'inhale vs exhale': numpy.take(dif_signal_above_zero, index_list, axis=0),
        'pif': numpy.append(getmaxbyindex(signal_data['flow'], index_list), 0),
        'pef': numpy.append(getminbyindex(signal_data['flow'], index_list), 0),
        'index_list': index_list
    })

    return candidate_breath_transitions


def extract_filtered_breaths_from_candidates(candidate_breath_transitions, analysis_parameters, local_logger=None):
    """
    Generates a refined list of breaths from candidate breath transitions,
    filtering based on PIF, PEF, and TI criteria from analysis_parameters.

    Parameters
    ----------
    candidate_breath_transitions : pd.DataFrame
        DataFrame containing annotated transitions between breaths.
    analysis_parameters : dict
        Dictionary of settings to use for analysis.
    local_logger : logging.Logger, optional
        Logger for logging messages. If None, a default logger is used.

    Raises
    ------
    Exception
        Raised if fewer than 30 breaths meet PIF, PEF, and TI criteria.

    Returns
    -------
    filtered_breaths : pd.DataFrame
        DataFrame describing parameters of the filtered breaths.
    """
    # Ensure a logger is available
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            local_logger.addHandler(handler)
            local_logger.setLevel(logging.INFO)

    # Filter based on minimum PIF and PEF
    filtered_transitions = candidate_breath_transitions[
        (candidate_breath_transitions['pif'] > float(analysis_parameters['minimum_PIF'])) |
        (candidate_breath_transitions['pef'] < float(analysis_parameters['minimum_PEF']))
    ]
    
    toggle_inhale = 0
    toggle_exhale = 1

    ts_list = list(filtered_transitions['ts'])
    ie_list = list(filtered_transitions['inhale vs exhale'])
    pif_list = list(filtered_transitions['pif'])
    pef_list = list(filtered_transitions['pef'])
    index_list = list(filtered_transitions['index_list'])
    inhale_filter = [0 for i in ts_list]
    exhale_filter = [0 for i in ts_list]

    # collect breaths based on toggle between inhale/exhale
    for i in range(len(ts_list)):
        if \
                ie_list[i] == 1 and \
                toggle_exhale == 1 and \
                pif_list[i] > float(analysis_parameters['minimum_PIF']):
            toggle_inhale = 1
            toggle_exhale = 0
            inhale_filter[i] = 1
        if \
                ie_list[i] == -1 and \
                toggle_inhale == 1 and \
                pef_list[i] < float(analysis_parameters['minimum_PEF']):
            toggle_inhale = 0
            toggle_exhale = 1
            exhale_filter[i] = 1

    number_of_full_breaths = min(sum(inhale_filter), sum(exhale_filter))
    ts_inhale = numpy.take(ts_list, numpy.nonzero(inhale_filter)[0])
    ts_exhale = numpy.take(ts_list, numpy.nonzero(exhale_filter)[0])
    il_inhale = numpy.take(index_list, numpy.nonzero(inhale_filter)[0])
    il_exhale = numpy.take(index_list, numpy.nonzero(exhale_filter)[0])
    duration_list = ts_exhale[:number_of_full_breaths] - \
        ts_inhale[:number_of_full_breaths]
    # filter based on minimum TI from Analysis Parameters
    ti_filter = numpy.array(
        duration_list > float(analysis_parameters['minimum_TI'])
        ).astype(int)
    filtered_breath_index_list = numpy.nonzero(ti_filter)[0]
    if len(filtered_breath_index_list) < 30:
        local_logger.error('insufficient breaths to score')
        raise Exception('fewer than 30 breaths meet PIF PEF and TI criteria')

    filtered_breaths = pandas.DataFrame()
    filtered_breaths['ts_inhale'] = numpy.take(
        ts_inhale, filtered_breath_index_list[:-1]
        )
    filtered_breaths['ts_exhale'] = numpy.take(
        ts_exhale, filtered_breath_index_list[:-1]
        )
    filtered_breaths['ts_end'] = numpy.take(
        ts_inhale, filtered_breath_index_list[1:]
        )
    filtered_breaths['il_inhale'] = numpy.take(
        il_inhale, filtered_breath_index_list[:-1]
        )
    filtered_breaths['il_exhale'] = numpy.take(
        il_exhale, filtered_breath_index_list[:-1]
        )
    filtered_breaths['il_end'] = numpy.take(
        il_inhale, filtered_breath_index_list[1:]
        )
    local_logger.info('{} breaths found'.format(len(filtered_breaths)))

    return filtered_breaths


def calculate_irreg_score(input_series):
    """
    takes a numpy compatible series and calculates an irregularity score
    using the formula |x[n]-X[n-1]| / X[n-1] * 100%. 
    A series of the irregularity scores will be returned.
    First value will be zero as it has no comparison to change from.


    Parameters
    ----------
    input_series : Pandas.DataSeries of Floats
        Data to use for Irreg Score Calculation

    Returns
    -------
    output_series : Pandas.DataSeries of Floats
        Series of Irreg Scores, paired to input_series

        
    """
    output_series = numpy.insert(
        numpy.multiply(
            numpy.divide(
                numpy.abs(
                    numpy.subtract(
                        list(input_series[1:]), list(input_series[:-1])
                    )
                ),
                list(input_series[:-1])
            ),
            100
        ),
        0, numpy.nan
    )
    return output_series


def calculate_moving_average(input_series, window, include_current=True):
    """
    Calculates a moving average of a series, with an option to exclude the current value
    from the average calculation.

    Parameters
    ----------
    input_series : pd.Series
        Data to use for moving average calculation.
    window : int
        Number of samples to use for the moving average.
    include_current : bool, optional
        Whether to include the 'middle' value in the moving average calculation.
        The default is True.

    Returns
    -------
    moving_average : pd.Series
        Moving average smoothed series paired to the input_series.
    """
    if include_current:
        moving_average = input_series.rolling(window, center=True, min_periods=1).mean()
    else:
        # Adjust for excluding the current value:
        # Calculate the sum over the window, subtract the current value, and divide by the adjusted count.
        total_sum = input_series.rolling(window, center=True, min_periods=1).sum()
        count = input_series.rolling(window, center=True, min_periods=1).count() - 1
        # Ensure we do not divide by zero
        count = count.apply(lambda x: max(x, 1))
        adjusted_sum = total_sum - input_series
        moving_average = adjusted_sum / count

    return moving_average


def beatcaller(CT, TS, signal_data, minRR=100, ecg_filter='1', ecg_invert='0', 
               analysis_parameters=None, sampling_rate=1000):
    """
    Extract R-R intervals and calculate heart rate from an ECG signal using segmented processing.

    Parameters:
    - CT (array-like): The raw ECG data.
    - TS (array-like): Timestamps corresponding to the ECG data points.
    - sampling_rate (float): The sampling rate of the ECG signal in Hz.
    - minRR (float): Minimum RR interval in seconds to consider for heart rate calculation.

    Returns:
    - DataFrame: DataFrame containing timestamps, RR intervals, and heart rates.
    """

    

    # Parameter overrides from analysis_parameters
    if analysis_parameters is not None:
        minRR = float(analysis_parameters.get('ecg_minRR', minRR))
        ecg_filter = str(analysis_parameters.get('ecg_filter', ecg_filter))
        ecg_invert = str(analysis_parameters.get('ecg_invert', ecg_invert))

    # Convert CT to a numpy array if it isn't one already
    CT = numpy.array(CT)
    
    # Invert ECG signal if required
    if ecg_invert == '1':
        CT = -CT
    
    #Calculate isoelectric line
    isoline = statistics.median(signal_data['ecg'])
    absthresh = isoline + 0.15
    
    # Identify peaks in the ECG signal; adjust parameters as necessary for your data
    peaks,_ = find_peaks(CT, height=absthresh, distance=100)  # Adjust 'distance' as needed    
    
    # Extract timestamps for the detected peaks
    timestamps_peaks = numpy.take(TS, peaks, axis=0)
    
    # Calculate RR intervals in seconds
    rr_intervals = numpy.diff(timestamps_peaks)
    
    # Calculate heart rate from rr intervals
    heart_rates = [60/ri for ri in rr_intervals]
    
    # Prepare dataframe to return
    beat_df = pandas.DataFrame({
        'ts': timestamps_peaks[:-1], # Exclude the last timestamp
        'RR': rr_intervals,
        'HR': heart_rates
        })
   
    return beat_df    


def calculate_basic_breath_parameters(
    signal_data,
    filtered_breaths,
    analysis_parameters,
    local_logger
    ):
    """
    Generates a datafrome containing a list of breaths with additional columns
    describing parameters for each breath.

    Parameters
    ----------
    signal_data : Pandas.DataFrame
        DataFrame containing signal data
    filtered_breaths : Pandas.DataFrame
        DataFrame containing filtered list of candidate breaths
    analysis_parameters : dict
        dictionary containing settings for analysis
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    breath_parameters : Pandas.DataFrame
        DataFrame containing annotated parameters for candidate breaths

    """
    # Initialize logger if not provided
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            local_logger.setLevel(logging.INFO)
            ch = logging.StreamHandler()
            ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            local_logger.addHandler(ch)
    
    # Copy the filtered breaths DataFrame to work with
    breath_parameters = filtered_breaths.copy()
    
    # Calculate breath duration
    breath_parameters['TI'] = breath_parameters['ts_exhale'] - breath_parameters['ts_inhale']
    breath_parameters['TE'] = breath_parameters['ts_end'] - breath_parameters['ts_exhale']
    breath_parameters['TT'] = breath_parameters['ts_end'] - breath_parameters['ts_inhale']
    breath_parameters['IS_TT'] = calculate_irreg_score(breath_parameters['TT'])
    breath_parameters['BPM'] = 60 / breath_parameters['TT']
    breath_parameters['IE_ratio'] = (breath_parameters['TI'] /breath_parameters['TE'])
    
    # Calculate sample interval based on timestamps
    sample_interval = numpy.diff(signal_data['ts']).mean()

    # Inspiratory and Expiratory Tidal Volumes (iTV and eTV)
    breath_parameters['iTV'] = numpy.multiply(getsumby2index(signal_data['flow'], filtered_breaths['il_inhale'], filtered_breaths['il_exhale']), sample_interval)
    breath_parameters['eTV'] = numpy.multiply(getsumby2index(-signal_data['flow'], filtered_breaths['il_exhale'], filtered_breaths['il_end']), sample_interval)
    
    # Irregularity Score for iTV
    breath_parameters['IS_iTV'] = calculate_irreg_score(breath_parameters['iTV'])

    # Peak Inspiratory and Expiratory Flows (PIF and PEF)
    breath_parameters['PIF'] = getmaxby2index(signal_data['flow'], breath_parameters['il_inhale'], breath_parameters['il_exhale'])
    breath_parameters['PEF'] = getmaxby2index(-signal_data['flow'], breath_parameters['il_exhale'], breath_parameters['il_end'])

    # Difference in Tidal Volumes (DVTV)
    breath_parameters['DVTV'] = numpy.divide(
        numpy.abs(numpy.subtract(numpy.abs(breath_parameters['iTV']), numpy.abs(breath_parameters['eTV']))),
        numpy.divide(numpy.add(numpy.abs(breath_parameters['iTV']),numpy.abs(breath_parameters['eTV'])), 2))

    # Local Thresholds for Apnea and Sigh Detection
    breath_parameters['apnea_local_threshold'] = calculate_moving_average(breath_parameters['TT'], int(analysis_parameters['apnea_window']), include_current=False)
    breath_parameters['sigh_local_threshold'] = calculate_moving_average(breath_parameters['iTV'], int(analysis_parameters['sigh_window']), include_current=False)

    # Percent X Calculation for BPM exceeding a threshold
    breath_parameters['per_X_calculation'] = calculate_moving_average(
        (
            breath_parameters['BPM'] >
            float(analysis_parameters['percent_X_value'])
            ).astype(int),
        int(analysis_parameters['percent_X_window']),
        include_current=True
        )
    
    # Define default values for missing analysis parameters
    default_values = {
        'minimum_apnea_time': 0.5,
        'minimum_apnea_fold_change': 2
        }

    # Check and update missing parameters with default values
    missing_keys = [key for key in default_values if key not in analysis_parameters]
    if missing_keys:
        local_logger.warning('Basic settings are missing for apnea values: ' + ', '.join(missing_keys))
        for key in missing_keys:
            local_logger.warning(f'Default value used for {key}')
            analysis_parameters[key] = default_values[key]
    
    # Calculate apnea detection flag
    breath_parameters['apnea'] = (
        (breath_parameters['TT'] > breath_parameters['apnea_local_threshold'] * float(analysis_parameters['minimum_apnea_fold_change'])) &
        (breath_parameters['TT'] > float(analysis_parameters['minimum_apnea_time']))
    ).astype(int)

    # Calculate sigh detection flag
    breath_parameters['sigh'] = (
        breath_parameters['iTV'] > breath_parameters['sigh_local_threshold'] * float(analysis_parameters['minimum_sigh_amplitude_x_local_VT'])
    ).astype(int)

    # Calculate average values for various gas concentrations and temperatures across breath intervals
    for param in ['o2', 'co2', 'corrected_o2', 'corrected_co2', 'temp', 'corrected_temp', 'Body_Temperature', 'mov_avg_vol']:
        breath_parameters[param] = get_avg_by_2_index(
            signal_data[param],
            filtered_breaths['il_inhale'],
            filtered_breaths['il_end']
        )

    # Calculate the number of breaths that do not meet the specified DVTV and per_X_calculation conditions
    breaths_to_remove = breath_parameters[(breath_parameters['DVTV'] >= float(analysis_parameters['maximum_DVTV'])) |
                                      (breath_parameters['per_X_calculation'] >= float(analysis_parameters['maximum_percent_X']))]

    # Calculate the number of breaths removed
    breaths_removed = len(breaths_to_remove)

    # Log the number of breaths removed
    local_logger.info(f'{breaths_removed} breaths removed from initial quality filter')

    return breath_parameters[
    (breath_parameters['DVTV'] < float(analysis_parameters['maximum_DVTV'])) &
    (breath_parameters['per_X_calculation'] < float(analysis_parameters['maximum_percent_X']))
    ]


def breath_caller(
        signal_data, analysis_parameters, local_logger
        ):
    """
    Generates a list of breaths with accompanying parameters starting from 
    flow signal data. This function is a wrapper that utilizes 
    extract_candidate_breath_transitions(), extract_filtered_breaths_from_candidates(),
    and calculate_basic_breath_parameters() to process the signal data.

    Parameters
    ----------
    signal_data : pd.DataFrame
        DataFrame containing signal data.
    analysis_parameters : dict
        Dictionary containing settings for analysis.
    local_logger : logging.Logger, optional
        Logger for tracking and logging the process. If None, a default logger is used.

    Returns
    -------
    basic_breath_parameters : pd.DataFrame
        DataFrame containing annotated breaths with calculated parameters.
    """
    # Initialize a local logger if none is provided
    if local_logger is None:
        local_logger = logging.getLogger(__name__)
        if not local_logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            local_logger.addHandler(handler)
            local_logger.setLevel(logging.INFO)

    try:
        local_logger.info('Detection of breaths started.')

        # Extract candidate breath transitions
        candidate_breath_transitions = extract_candidate_breath_transitions(signal_data)

        # Filter breaths based on candidate transitions
        filtered_breaths = extract_filtered_breaths_from_candidates(
            candidate_breath_transitions, analysis_parameters, local_logger)

        # Calculate detailed parameters for each breath
        basic_breath_parameters = calculate_basic_breath_parameters(
            signal_data, filtered_breaths, analysis_parameters, local_logger)

        local_logger.info('Detection of breaths completed successfully.')
        return basic_breath_parameters
    except Exception as e:
        local_logger.error(f'Error in breath detection process: {e}')
        raise


def create_auto_calibration_dict(
        signal_data: pandas.DataFrame,
        breath_list: pandas.DataFrame,
        auto_criteria: Dict,
        timestamp_dict: Dict,
        local_logger: logging.Logger
        ) -> Dict:
    """
    Generates a dictionary describing the calibration values for O2 and CO2
    measurements so that calibration curves can be used to improve accuracy of
    calorimetry values if possible. Calibration periods are identified through
    auto_criteria.

    Parameters
    ----------
    signal_data : Pandas.DataFrame
        DataFrame containing signal data
    breath_list : Pandas.DataFrame
        DataFrame containing annotated breaths
    auto_criteria : Dict
        Dict with info describing criteria for automated selection of breathing
        for 'calm quiet' breathing
    timestamp_dict : Dict
        Dictionary containing timestamps and text used to describe them. 
        Captured from the commends in the signal_data
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    auto_calibration_dict : Dict
        Dict describing baseline measures and parameters for calibration

    """
    if auto_criteria is None:
        return None

    auto_calibration_dict = {
        'o2': {
            'nom': [],
            'voltage': []
        },
        'co2': {
            'nom': [],
            'voltage': []
        }
    }
    cal_list = breath_list.copy()
    all_keys = auto_criteria['Key'].unique()

    cal_keys = auto_criteria[
        auto_criteria['AUTO_IND_GAS_CAL'] == 1
        ][['Key', 'Alias', 'CAL_O2', 'CAL_CO2']]
    cal_keys['Key_and_Alias'] = \
        cal_keys['Key'].astype(str) + cal_keys['Alias'].astype(str)

    reverse_timestamp_dict = make_reverse_timestamp_dict(timestamp_dict)

    for key_and_alias in cal_keys['Key_and_Alias']:

        key = cal_keys[
            cal_keys['Key_and_Alias'] == key_and_alias
        ]['Key'].values[0]
        alias = cal_keys[
            cal_keys['Key_and_Alias'] == key_and_alias
        ]['Alias'].values[0]
        if key not in reverse_timestamp_dict:
            continue

        current_auto_criteria = \
            auto_criteria[
                (auto_criteria['Key'] == key) & (
                    auto_criteria['Alias'] == alias)
            ].to_dict('records')[0]

        block_start, block_stop = extract_block_bounds(
            signal_data, reverse_timestamp_dict, all_keys, key
        )
        verified_block_start = max(
            block_start + float(current_auto_criteria['after_start']),
            block_stop - float(current_auto_criteria['within_end'])
        )
        verified_block_stop = min(
            block_stop - float(current_auto_criteria['before_end']),
            block_start + float(current_auto_criteria['within_start'])
        )

        try:
            nom_o2, min_v_o2, max_v_o2 = cal_keys[
                cal_keys['Key_and_Alias'] == key_and_alias
            ]['CAL_O2'].values[0].split(',')
            o2_filter = \
                (cal_list['o2'] >= float(min_v_o2)) & \
                (cal_list['o2'] <= float(max_v_o2))

        except:
            nom_o2, min_v_o2, max_v_o2 = None, None, None
            o2_filter = cal_list['o2'] == cal_list['o2']
        try:
            nom_co2, min_v_co2, max_v_co2 = cal_keys[
                cal_keys['Key_and_Alias'] == key_and_alias
            ]['CAL_CO2'].values[0].split(',')
            co2_filter = \
                (cal_list['co2'] >= float(min_v_co2)) & \
                (cal_list['co2'] <= float(max_v_co2))

        except:
            nom_co2, min_v_co2, max_v_co2 = None, None, None
            co2_filter = cal_list['co2'] == cal_list['co2']

        cal_list['filter_{}'.format(key_and_alias)] = 0
        cal_list.loc[
            (cal_list['TT'] >= float(current_auto_criteria['min_TT'])) &
            (cal_list['TT'] <= float(current_auto_criteria['max_TT'])) &
            (cal_list['DVTV'] <= float(current_auto_criteria['max_DVTV'])) &
            (cal_list['per_X_calculation'] <=
             float(current_auto_criteria['max_pX'])) &
            (cal_list['iTV'] >= float(current_auto_criteria['min_TV'])) &
            (cal_list['eTV'] >= float(current_auto_criteria['min_TV'])) &
            (o2_filter) &
            (co2_filter) &
            (cal_list['mov_avg_vol'] <=
             float(current_auto_criteria['vol_mov_avg_drift'])) &
            (cal_list['ts_inhale'] >= verified_block_start) &
            (cal_list['ts_end'] <= verified_block_stop),
            'filter_{}'.format(key_and_alias)
        ] = 1
        
        if sum(cal_list['filter_{}'.format(key_and_alias)]) == 0:
            if sum(
                (cal_list['ts_inhale'] >= verified_block_start) & \
                (cal_list['ts_end'] <= verified_block_stop) & \
                (o2_filter)
            ) < int(current_auto_criteria['min_bout']):
                local_logger.warning(
                    'insufficient gas calibration breaths satisfying O2 filter'
                )
            if sum(
                (cal_list['ts_inhale'] >= verified_block_start) & \
                (cal_list['ts_end'] <= verified_block_stop) & \
                (co2_filter)
            ) < int(current_auto_criteria['min_bout']):
                local_logger.warning(
                    'insufficient gas calibration breaths satisfying CO2 filter'
                )
            if sum(
                (cal_list['ts_inhale'] >= verified_block_start) & \
                (cal_list['ts_end'] <= verified_block_stop) & \
                (cal_list['mov_avg_vol'] <=
                 float(current_auto_criteria['vol_mov_avg_drift']))
            ) < int(current_auto_criteria['min_bout']):
                local_logger.warning(
                    'insufficient gas calibration breaths satisfying MAV filter'
                )
            
        biggest_chunk = get_the_biggest_chunk(
            get_the_chunks(
                cal_list['filter_{}'.format(key_and_alias)],
                1,
                int(current_auto_criteria['min_bout'])
            )
        )
        if biggest_chunk is None or biggest_chunk == []:
            local_logger.warning(
                'expected gas calibration data in ' + \
                'AUTO {} | {} - '.format(key, alias) + \
                'no suitable data found'
            )
            continue
        else:
            chunk_filter = \
                (cal_list['ts_inhale'] >=
                 cal_list.iloc[biggest_chunk[0][0]]['ts_inhale']) & \
                (cal_list['ts_end'] <=
                 cal_list.iloc[biggest_chunk[0][1]]['ts_end'])

        if nom_o2 is not None:
            auto_calibration_dict['o2']['nom'].append(float(nom_o2))
            auto_calibration_dict['o2']['voltage'].append(
                sum(
                    cal_list[chunk_filter]['o2'] *
                    cal_list[chunk_filter]['TT']
                ) / cal_list[chunk_filter]['TT'].sum()
            )
        if nom_co2 is not None:
            auto_calibration_dict['co2']['nom'].append(float(nom_co2))
            # create weighted average (TT) to generate average over sampling
            # time rather than per breath
            auto_calibration_dict['co2']['voltage'].append(
                sum(
                    cal_list[chunk_filter]['co2'] *
                    cal_list[chunk_filter]['TT']
                ) / cal_list[chunk_filter]['TT'].sum()
            )
    return auto_calibration_dict


def extract_block_bounds(
        signal_data: pandas.DataFrame,
        reverse_timestamp_dict: Dict,
        all_keys: List,
        key: str
        ) -> Tuple[float, float]:
    """
    Identifies the boundaries of blocks in data using timestamp information.
    
    Parameters
    ----------
    signal_data : Pandas.DataFrame
        DataFrame containing signal data
    reverse_timestamp_dict : Dict
        Dictionary containing text used to describe timestamps and the 
        actual timestamps themselves. Built from timestamp_dict. 
    all_keys : List
        Keys from reverse_timestamp_dict (not used)
    key : str
        Key to test in reverse_timestamp_dict

    Returns
    -------
    block_start : float
        Timestamp describing the start of a block interval
    block_end : float
        Timestamp describing the end of a block interval

    """
    block_start = reverse_timestamp_dict[key]
    block_end = min(
        [
            reverse_timestamp_dict[k]
            for k in reverse_timestamp_dict
            if reverse_timestamp_dict[k] > block_start
        ] + [signal_data['ts'].max()]
    )
    return block_start, block_end


def create_man_calibration_dict(
        signal_data: pandas.DataFrame,
        manual_selection: Dict,
        plyuid: str
        ) -> Dict:
    """
    Generates a dictionary describing the calibration values for O2 and CO2
    measurements so that calibration curves can be used to improve accuracy of
    calorimetry values if possible. Calibration periods are identified through
    manual_selection.

    Parameters
    ----------
    signal_data : Pandas.DataFrame
        DataFrame containing signal data
    manual_selection : Dict
        Dict with info describing manually selected start and stop points for
        breathing of interest
    plyuid : str
        Serial identifier for data collection event

    Returns
    -------
    man_calibration_dict : Dict
        Dict describing baseline measures and parameters for calibration from
        manual selection

    """
    
    if manual_selection is None:
        return None

    man_calibration_dict = {
        'o2': {
            'nom': [],
            'voltage': []
            },
        'co2': {
            'nom': [],
            'voltage': []
            }
        }

    for segment_index in manual_selection[
            (manual_selection['MAN_IND_GAS_CAL'] == 1) &
            (manual_selection['PLYUID'].astype(str) == plyuid)
            ].index:

        nom_o2 = manual_selection.iloc[segment_index]['CAL_O2']
        nom_co2 = manual_selection.iloc[segment_index]['CAL_CO2']

        if not math.isnan(nom_o2):
            man_calibration_dict['o2']['nom'].append(float(nom_o2))
            voltage_o2 = signal_data[
                (signal_data['ts'] >= manual_selection.iloc[segment_index]['start']) &
                (signal_data['ts'] <= manual_selection.iloc[segment_index]['stop'])
            ]['o2'].mean()
            man_calibration_dict['o2']['voltage'].append(voltage_o2)

        if not math.isnan(nom_co2):
            man_calibration_dict['co2']['nom'].append(float(nom_co2))
            voltage_co2 = signal_data[
                (signal_data['ts'] >= manual_selection.iloc[segment_index]['start']) &
                (signal_data['ts'] <= manual_selection.iloc[segment_index]['stop'])
            ]['co2'].mean()
            man_calibration_dict['co2']['voltage'].append(voltage_co2)
    
    return man_calibration_dict


def get_the_chunks(chunkable_list: List, value, min_chunk: int) -> List[tuple]:
    """
    Generates a list of tuples describing indexes where values in an input list
    match a value. Only continuous chunks that are at least min_chunk in 
    length are included in the returned list.
    
    Parameters
    ----------
    chunkable_list : List
        List containing values to search through
    value : ...
        Value that is being searched for in the chunkable list
    min_chunk : int
        Integer minimum chunk size (consecutive list entries matching value)

    Returns
    -------
    chunks : List of tuples
        Returns a list of tuples describing the indexes of chunks with the 
        value that meet the minimum chunk size

    """
    status = 0
    chunks = []
    chunk_start = 0
    chunk_stop = 1
    if len(chunkable_list) < min_chunk:
        return chunks

    for i in range(len(chunkable_list)):
        if chunkable_list[i] == value:
            if status == 0:
                chunk_start = i
                status = 1
            else:
                chunk_stop = i
        else:
            if status == 1:
                if chunk_stop - chunk_start >= min_chunk:
                    chunks.append((chunk_start, chunk_stop))
                status = 0
    if len(chunks) == 0 and status == 1:
        chunks.append((chunk_start, chunk_stop))
    return chunks


def get_the_biggest_chunk(chunk_list: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Compares entries in a list of tuples and returns the tuple that has the 
    largest difference between the values of the tuple. A list containing the 
    tuple with the largest difference is returned.
    
    Parameters
    ----------
    chunk_list : List of Tuples of Ints (or Floats)
        List of Tuples of Ints describing indexes of selections

    Returns
    -------
    biggest_chunk : List of Tuples of Ints (or Floats)
        List of Tuples describing the boundary of the largest 'chunk' in the chunk_list

    """
    biggest_chunk = []
    if len(chunk_list) == 0:
        return biggest_chunk

    biggest_chunk_size = 0

    for c in chunk_list:
        if c[1] - c[0] > biggest_chunk_size:
            biggest_chunk = [c]
            biggest_chunk_size = c[1] - c[0]
    return biggest_chunk


def reconcile_calibration_selections(
        auto_calibration_dict: Dict[str, Dict[str, List[Union[float, int]]]],
        man_calibration_dict: Dict[str, Dict[str, List[Union[float, int]]]],
        local_logger
        ) -> Tuple[Dict[Union[int, float, str], Union[float, int]], 
                   Dict[Union[int, float, str], Union[float, int]]]:
    """
    Checks for the presence of calibration information derived from automated 
    or manual selections. If both are present, the manual selections take
    precedent.
    
    Parameters
    ----------
    auto_calibration_dict : Dict
        Dict describing baseline measures and parameters to use for calibration
        from automated selection
    man_calibration_dict : Dict
        Dict describing baseline measures and parameters to use for calibration
        from manual selection
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    o2_cal_dict : Dict
        Dict describing baseline measures and parameters to use for calibration
    co2_cal_dict : Dict
        Dict describing baseline measures and parameters to use for calibration

    """
    o2_cal_dict = {}
    temp_o2_cal_dict = {}
    co2_cal_dict = {}
    temp_co2_cal_dict = {}

    if auto_calibration_dict is not None and man_calibration_dict is not None:
        local_logger.warning('gas calibration data present in both ' +
            'automated and manual selections - manual selections ' +
            'will override automated selections if they correscpond ' +
            'to the same nominal concentrations')

    if auto_calibration_dict is not None:

        for i in range(len(auto_calibration_dict['o2']['nom'])):
            if auto_calibration_dict['o2']['nom'][i] not in o2_cal_dict:
                o2_cal_dict[auto_calibration_dict['o2']['nom'][i]] = [
                    auto_calibration_dict['o2']['voltage'][i]]
            else:
                o2_cal_dict[auto_calibration_dict['o2']['nom'][i]].append(
                    auto_calibration_dict['o2']['voltage'][i]
                    )
        for i in range(len(auto_calibration_dict['co2']['nom'])):
            if auto_calibration_dict['co2']['nom'][i] not in co2_cal_dict:
                co2_cal_dict[auto_calibration_dict['co2']['nom'][i]] = [
                    auto_calibration_dict['co2']['voltage'][i]]
            else:
                co2_cal_dict[auto_calibration_dict['co2']['nom'][i]].append(
                    auto_calibration_dict['co2']['voltage'][i]
                    )

    if man_calibration_dict is not None:

        for i in range(len(man_calibration_dict['o2']['nom'])):
            if man_calibration_dict['o2']['nom'][i] not in temp_o2_cal_dict:
                temp_o2_cal_dict[man_calibration_dict['o2']['nom'][i]] = [
                    man_calibration_dict['o2']['voltage'][i]]
            else:

                temp_o2_cal_dict[man_calibration_dict['o2']['nom'][i]].append(
                    man_calibration_dict['o2']['voltage'][i]
                    )
        for i in range(len(man_calibration_dict['co2']['nom'])):
            if man_calibration_dict['co2']['nom'][i] not in temp_co2_cal_dict:
                temp_co2_cal_dict[man_calibration_dict['co2']['nom'][i]] = [
                    man_calibration_dict['co2']['voltage'][i]]
            else:

                temp_co2_cal_dict[man_calibration_dict['co2']['nom'][i]].append(
                    man_calibration_dict['co2']['voltage'][i]
                    )
        for k in temp_o2_cal_dict:
            o2_cal_dict[k] = temp_o2_cal_dict[k]
        for k in temp_co2_cal_dict:
            co2_cal_dict[k] = temp_co2_cal_dict[k]

    o2_cal_dict['avg'] = {}
    co2_cal_dict['avg'] = {}

    for k in o2_cal_dict:
        if k != 'avg':

            o2_cal_dict['avg'][k] = (sum(o2_cal_dict[k]) /
                len(o2_cal_dict[k]))
    for k in co2_cal_dict:
        if k != 'avg':
            co2_cal_dict['avg'][k] = (sum(co2_cal_dict[k]) /
                len(co2_cal_dict[k]))

    return o2_cal_dict, co2_cal_dict


def apply_gas_calibration(
        signal_data: pandas.DataFrame,
        breath_list: pandas.DataFrame,
        o2_cal_dict: Dict[Union[int, float, str], Union[float, int]],
        co2_cal_dict: Dict[Union[int, float, str], Union[float, int]],
        local_logger
        ) -> Tuple[pandas.Series, pandas.Series, pandas.Series, pandas.Series]:
    """
    Creates a calibration curve for O2 and CO2 values and applies it to O2 and
    CO2 data to improve accuracy of calorimetry measures.

    Parameters
    ----------
    signal_data : pd.DataFrame
        DataFrame containing signal data
    breath_list : pd.DataFrame
        DataFrame containing annotated breaths
    o2_cal_dict : Dict
        Dict describing O2 Calibration Parameters
    co2_cal_dict : Dict
        Dict describing CO2 Calibration Parameters
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    signal_calibrated_o2 : pd.Series
        Calibrated O2 matched to signal_data
    signal_calibrated_co2 : pd.Series
        Calibrated CO2 matched to signal_data
    breath_list_calibrated_o2 : pd.Series
        Calibrated O2 matched to breath_list
    breath_list_calibrated_co2 : pd.Series
        Calibrated CO2 matched to breath_list

    """

    if 0 not in o2_cal_dict and 0 < len(o2_cal_dict) < 2:
        local_logger.warning(
            'only one O2 calibration data point present - ' +
            '2 point curve will be created with assumption of 0 volts = 0%'
            )
        O2_curve = {0: 0}
    else: O2_curve = {}
    if 0 not in co2_cal_dict and 0 < len(co2_cal_dict) < 2:
        local_logger.warning(
            'only one CO2 calibration data point present - ' +
            '2 point curve will be created with assumption of 0 volts = 0%'
            )
        CO2_curve = {0: 0}
    else: CO2_curve = {}

    for k in o2_cal_dict:
        O2_curve[k] = o2_cal_dict[k]
    for k in co2_cal_dict:
        CO2_curve[k] = co2_cal_dict[k]

    O2_standards = list(O2_curve['avg'].keys())
    CO2_standards = list(CO2_curve['avg'].keys())

    if len(O2_standards) > 1:
        O2_slope, O2_intercept, r, p, s = stats.linregress(
            O2_standards, [O2_curve['avg'][k] for k in O2_standards]
            )
        signal_calibrated_o2 = (signal_data['o2']-O2_intercept)/O2_slope
        breath_list_calibrated_o2 = (breath_list['o2']-O2_intercept)/O2_slope
    else:
        signal_calibrated_o2 = signal_data['corrected_o2']
        breath_list_calibrated_o2 = breath_list['corrected_o2']
        local_logger.warning(
            'unable to calibrate O2 - insufficient data - ' +
            'correction factor multiplier provided in analysis parameters used'
            )

    if len(CO2_standards) > 1:
        CO2_slope, CO2_intercept, r, p, s = stats.linregress(
            CO2_standards, [CO2_curve['avg'][k] for k in CO2_standards]
            )
        signal_calibrated_co2 = (signal_data['co2']-CO2_intercept)/CO2_slope
        breath_list_calibrated_co2 = \
            (breath_list['co2']-CO2_intercept)/CO2_slope
    else:
        signal_calibrated_co2 = signal_data['corrected_co2']
        breath_list_calibrated_co2 = breath_list['corrected_co2']
        local_logger.warning(
            'unable to calibrate CO2 - insufficient data - ' +
            'correction factor multiplier provided in analysis parameters used'
            )

    return signal_calibrated_o2, signal_calibrated_co2, \
        breath_list_calibrated_o2, breath_list_calibrated_co2


def calibrate_gas(
        signal_data: pandas.DataFrame,
        breath_list: pandas.DataFrame,
        auto_criteria: Dict,
        manual_selection: Dict,
        plyuid: str,
        timestamp_dict: Dict,
        local_logger
        ) -> Tuple[pandas.Series, pandas.Series, pandas.Series, pandas.Series]:
    """
    Collects calibration information from auto_criteria, manual_selection, and 
    signal_data to improve accuracy of calorimetry data (O2 and CO2) by 
    creating a standard curve of gas concentrations (if possible) and adjusting
    O2 and CO2 values accordingly.

    Parameters
    ----------
    signal_data : pd.DataFrame
        DataFrame containing signal data
    breath_list : pd.DataFrame
        DataFrame containing annotated breaths
    auto_criteria : Dict
        Dict with info describing criteria for automated selection of breathing
        for 'calm quiet' breathing
    manual_selection : Dict
        Dict with info describing manually selected start and stop points for
        breathing of interest
    plyuid : str
        Serial identifier for data collection event
    timestamp_dict : Dict
        Dictionary containing timestamps and text used to describe them. 
        Captured from the comments in the signal_data
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    signal_calibrated_o2 : pd.Series
        Calibrated O2 matched to signal_data
    signal_calibrated_co2 : pd.Series
        Calibrated CO2 matched to signal_data
    breath_list_calibrated_o2 : pd.Series
        Calibrated O2 matched to breath_list
    breath_list_calibrated_co2 : pd.Series
        Calibrated CO2 matched to breath_list

    """

    auto_calibration_dict = create_auto_calibration_dict(
        signal_data,
        breath_list,
        auto_criteria,
        timestamp_dict,
        local_logger
        )

    man_calibration_dict = create_man_calibration_dict(
        signal_data,
        manual_selection,
        plyuid
        )

    o2_cal_dict, co2_cal_dict = reconcile_calibration_selections(
        auto_calibration_dict,
        man_calibration_dict,
        local_logger
        )

    signal_calibrated_o2, signal_calibrated_co2, \
    breath_list_calibrated_o2, breath_list_calibrated_co2 = \
        apply_gas_calibration(
        signal_data,
        breath_list,
        o2_cal_dict,
        co2_cal_dict,
        local_logger
        )

    return signal_calibrated_o2, signal_calibrated_co2, \
        breath_list_calibrated_o2, breath_list_calibrated_co2


def make_filter_from_chunk_list(ts_inhale_and_ts_end, chunk_list):
    """
    Creates a series where contents are a binary filter where values are 1 if 
    contained within indexes specified in the chunk_list tuples, otherwise they 
    are 0. A second series where all indexes specified by a tuple are labeled 
    with the index of the first member.

    Parameters
    ----------
    ts_inhale_and_ts_end : pd.DataFrame
        DataFrame describing transition times for inhalation and end of breath
    chunk_list : List of tuples of Ints
        List of Tuples of Floats describing Indexes that should be included 
        in the filter

    Returns
    -------
    ts_inhale_and_ts_end['filter'] : pd.DataSeries
        Binary filter based on the input timelist
    
    ts_inhale_and_ts_end['selection_id'] : pd.DataSeries
        Selection ID metadata for the input timelist

    """
    df_inhale_and_end = pandas.DataFrame(ts_inhale_and_ts_end.copy())
    df_inhale_and_end['filter'] = 0
    df_inhale_and_end['selection_id'] = ""
    for chunk in chunk_list:
        df_inhale_and_end.loc[
            (df_inhale_and_end['ts_inhale'] >=
                 df_inhale_and_end.iloc[chunk[0]]['ts_inhale']) &
            (df_inhale_and_end['ts_end'] <=
                 df_inhale_and_end.iloc[chunk[1]]['ts_end']),
            'filter'
            ] = 1
        df_inhale_and_end.loc[
            (df_inhale_and_end['ts_inhale'] >=
                 df_inhale_and_end.iloc[chunk[0]]['ts_inhale']) &
            (df_inhale_and_end['ts_end'] <=
                 df_inhale_and_end.iloc[chunk[1]]['ts_end']),
            'selection_id'
            ] = chunk[0]
    return df_inhale_and_end['filter'], df_inhale_and_end['selection_id']


def make_filter_from_time_list(ts_inhale_and_ts_end, time_list):
    """
    Creates a series where contents are a binary filter where values are 1 if 
    contained within timestamps specified in the ts_inhale_and_ts_end dataframe,
    otherwise they are 0. A second series is also returned which labels
    indexes associated with member breaths according to their selection_id.

    Parameters
    ----------
    ts_inhale_and_ts_end : pd.DataFrame
        DataFrame describing transition times for inhalation and end of breath
    time_list : List of tuples of Floats
        List of Tuples of Floats describing beginning and ending of timestamps
        that should be included in the filter

    Returns
    -------
    df_inhale_and_end['filter'] : pd.Series
        Binary filter based on the input timestamp list
    
    df_inhale_and_end['selection_id'] : pd.Series
        Selection ID metadata for the input timestamp list

    """
    df_inhale_and_end = pandas.DataFrame(ts_inhale_and_ts_end.copy())
    df_inhale_and_end['filter'] = 0
    df_inhale_and_end['selection_id'] = ""
    for timestamps in time_list:
        df_inhale_and_end.loc[
            (df_inhale_and_end['ts_inhale'] >= timestamps[0]) &
            (df_inhale_and_end['ts_end'] <= timestamps[1]),
            'filter'
            ] = 1
        df_inhale_and_end.loc[
            (df_inhale_and_end['ts_inhale'] >= timestamps[0]) &
            (df_inhale_and_end['ts_end'] <= timestamps[1]),
            'selection_id'
            ] = timestamps[0]
    
    return df_inhale_and_end['filter'], df_inhale_and_end['selection_id']


def create_filters_for_automated_selections(
        signal_data,
        breath_list,
        auto_criteria,
        timestamp_dict,
        analysis_parameters,
        local_logger
        ):
    """
    Creates a dictionary which describes the boundaries of automated selection 
    intervals. Contents include a breathlist filter describing membership to
    particular intervals as well as related timestamps.

    Parameters
    ----------
    signal_data : pd.DataFrame
        DataFrame containing signal data
    breath_list : pd.DataFrame
        DataFrame containing annotated breaths
    auto_criteria : Dict
        Dict with info describing criteria for automated selection of breathing
        for 'calm quiet' breathing
    timestamp_dict : Dict
        Dictionary containing timestamps and text used to describe them. 
        Captured from the commends in the signal_data
    analysis_parameters : Dict
        Dict of settings to use for analysis
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    block_dict : Dict
        Nested Dict describing boundaries of selected blocks and related 
        information such as a filter indicating member breaths, timestamp 
        boundaries, etc.

    """
    block_dict = {}
    if auto_criteria is None:
        return block_dict
    
    automated_block = pandas.DataFrame(breath_list['ts_inhale'].copy())
    automated_selection = pandas.DataFrame(breath_list['ts_inhale'].copy())
    automated_biggest_block = pandas.DataFrame(breath_list['ts_inhale'].copy())

    all_keys = auto_criteria['Key'].unique()

    auto_keys = auto_criteria[['Key', 'Alias']]
    auto_keys['key_and_alias'] = \
        auto_keys['Key'].astype(str) + auto_keys['Alias'].astype(str)

    reverse_timestamp_dict = make_reverse_timestamp_dict(timestamp_dict, key_list=all_keys)

    for key_and_alias in auto_keys['key_and_alias']:

        key = auto_keys[
            auto_keys['key_and_alias'] == key_and_alias
            ]['Key'].values[0]
        alias = auto_keys[
            auto_keys['key_and_alias'] == key_and_alias
            ]['Alias'].values[0]
        if key not in reverse_timestamp_dict:
            continue

        block_dict[key_and_alias] = {}

        current_auto_criteria = \
            auto_criteria[
                (auto_criteria['Key'] == key) & (
                    auto_criteria['Alias'] == alias)
                ].to_dict('records')[0]

        block_start, block_stop = extract_block_bounds(
            signal_data, reverse_timestamp_dict, all_keys, key
            )
        verified_block_start = max(
            block_start + float(current_auto_criteria['after_start']),
            block_stop - float(current_auto_criteria['within_end'])
            )
        verified_block_stop = min(
            block_stop - float(current_auto_criteria['before_end']),
            block_start + float(current_auto_criteria['within_start'])
            )

        high_chamber_temp_filter = breath_list['Body_Temperature'] < \
            breath_list['corrected_temp']
        
        if 'calibrated_TV' in breath_list.columns and \
                current_auto_criteria['max_calibrated_TV'] != 'off':
            high_TV_filter = breath_list['calibrated_TV'] < \
                float(current_auto_criteria['max_calibrated_TV'])
        else:
            high_TV_filter = \
                breath_list['ts_inhale'] == breath_list['ts_inhale']
        
        if 'VE_per_VO2' in breath_list.columns and \
                current_auto_criteria['max_VEVO2'] != 'off':
            high_VEVO2_filter = breath_list['VE_per_VO2'] < \
                float(current_auto_criteria['max_VEVO2'])
        else:
            high_VEVO2_filter = \
                breath_list['ts_inhale'] == breath_list['ts_inhale']
        
        if 'VO2_per_gram' in breath_list.columns and \
                current_auto_criteria.get('min_VO2pg','off') != 'off':
            low_VO2pg_filter = breath_list['VO2_per_gram'] > \
                float(current_auto_criteria['min_VO2pg'])
        else: 
            low_VO2pg_filter = \
                breath_list['ts_inhale'] == breath_list['ts_inhale']
                
        if 'VCO2_per_gram' in breath_list.columns and \
                current_auto_criteria.get('min_VCO2pg','off') != 'off':
            low_VCO2pg_filter = breath_list['VCO2_per_gram'] > \
                float(current_auto_criteria['min_VCO2pg'])
        else:
            low_VCO2pg_filter = \
                breath_list['ts_inhale'] == breath_list['ts_inhale']
                
        if 'VO2' in breath_list.columns and \
                current_auto_criteria.get('min_VO2','off') != 'off':
            low_VO2_filter = breath_list['VO2'] > \
                float(current_auto_criteria['min_VO2'])
        else: 
            low_VO2_filter = \
                breath_list['ts_inhale'] == breath_list['ts_inhale']
                
        if 'VCO2' in breath_list.columns and \
                current_auto_criteria.get('min_VCO2','off') != 'off':
            low_VCO2_filter = breath_list['VCO2'] > \
                float(current_auto_criteria['min_VCO2'])
        else:
            low_VCO2_filter = \
                breath_list['ts_inhale'] == breath_list['ts_inhale']

        automated_block[key_and_alias] = 0
        automated_selection[key_and_alias] = 0
        automated_biggest_block[key_and_alias] = 0

        automated_block.loc[
            (current_auto_criteria['min_CO2'] <= breath_list['calibrated_co2']) &
            (current_auto_criteria['max_CO2'] >= breath_list['calibrated_co2']) &
            (current_auto_criteria['min_O2'] <= breath_list['calibrated_o2']) &
            (current_auto_criteria['max_O2'] >= breath_list['calibrated_o2']) &
            (breath_list['ts_inhale'] >= verified_block_start) &
            (breath_list['ts_end'] <= verified_block_stop),
            key_and_alias
            ] = 1
        automated_selection.loc[
            (current_auto_criteria['min_CO2'] <= breath_list['calibrated_co2']) &
            (current_auto_criteria['max_CO2'] >= breath_list['calibrated_co2']) &
            (current_auto_criteria['min_O2'] <= breath_list['calibrated_o2']) &
            (current_auto_criteria['max_O2'] >= breath_list['calibrated_o2']) &
            (breath_list['ts_inhale'] >= verified_block_start) &
            (breath_list['ts_end'] <= verified_block_stop) &
            (current_auto_criteria['min_TT'] <= breath_list['TT']) &
            (current_auto_criteria['max_TT'] >= breath_list['TT']) &
            (current_auto_criteria['max_DVTV'] >= breath_list['DVTV']) &
            (current_auto_criteria['max_pX'] >=
             breath_list['per_X_calculation']) &
            (current_auto_criteria['min_TV'] <= breath_list['iTV']) &
            (current_auto_criteria['min_TV'] <= breath_list['eTV']) &
            (current_auto_criteria['include_apnea'] >= breath_list['apnea']) &
            (current_auto_criteria['include_sigh'] >= breath_list['sigh']) &
            (current_auto_criteria['include_high_chamber_temp'] >=
             high_chamber_temp_filter) &
            (current_auto_criteria['vol_mov_avg_drift'] >=
             breath_list['mov_avg_vol']) &
            (high_VEVO2_filter == 1) &
            (high_TV_filter == 1) &
            (low_VO2_filter) &
            (low_VCO2_filter) &
            (low_VO2pg_filter) &
            (low_VCO2pg_filter),
            key_and_alias
            ] = 1

        blocks = get_the_chunks(
            automated_block[key_and_alias],
            1,
            int(current_auto_criteria['min_bout'])
            )
        biggest_block = get_the_biggest_chunk(blocks)
        selection_list = get_the_chunks(
            automated_selection[key_and_alias],
            1,
            int(current_auto_criteria['min_bout'])
            )
        if biggest_block == []:
            local_logger.warning(
                'AUTO {} | {} - '.format(key, alias) +
                'no suitable data found'
            )
            continue
        else:
            block_dict[key_and_alias]['biggest_block'], \
            block_dict[key_and_alias]['block_id'] = \
                make_filter_from_chunk_list(
                    breath_list[['ts_inhale', 'ts_end']],
                    biggest_block
                    )

            block_dict[key_and_alias]['selection_filter'], \
            block_dict[key_and_alias]['selection_id'] = \
                make_filter_from_chunk_list(
                    breath_list[['ts_inhale', 'ts_end']],
                    selection_list
                    )

    return block_dict


def create_filters_for_manual_selections(
        breath_list,
        manual_selection,
        plyuid,
        logger
        ):
    """
    Creates a dictionary which describes the boundaries of user selected 
    intervals. Contents include a breathlist filter describing membership to
    particular intervals as well as related timestamps.

    Parameters
    ----------
    breath_list : Pandas.DataFrame
        DataFrame containing annotated breaths
    manual_selection : Dict
        Dict with info describing manually selected start and stop points for
        breathing of interest
    plyuid : string
        serial identifier for data collection event 
        (plethysmography session id)
    logger : instance of logging.logger (optional)

    Returns
    -------
    block_dict : Nested Dict
        Nested Dict describing boundaries of selected blocks and related 
        information such as a filter indicating member breaths, timestamp 
        boundaries, etc.

    """

    block_dict = {}
    if manual_selection is None:
        return block_dict
    
    alias_dict = {}
    for segment_index in \
        manual_selection[
            (manual_selection['PLYUID'].astype(str) == plyuid)
            ].index:
        chunk_start = manual_selection.iloc[segment_index]['start']
        chunk_stop = manual_selection.iloc[segment_index]['stop']

        if manual_selection.iloc[segment_index]['Alias'] not in alias_dict:
            alias_dict[manual_selection.iloc[segment_index]['Alias']] = []

        alias_dict[manual_selection.iloc[segment_index]
            ['Alias']].append([chunk_start, chunk_stop])

    for alias in alias_dict:
        selection_filter, selection_id = make_filter_from_time_list(
            breath_list[['ts_inhale', 'ts_end']],
            alias_dict[alias])
        
        block_dict[alias] = {
            'selection_filter': selection_filter,
            'selection_id': selection_id
            }
    return block_dict


def collect_calibration_parameters(
        breath_list,
        analysis_parameters,
        animal_metadata,
        plyuid,
        auto_criteria,
        manual_selection,
        automated_selections_filters,
        manual_selections_filters,
        local_logger
        ):
    """
    Provides a dictionary containing parameters relevent to calibration and 
    normalization of volume and VO2,VCO2. 

    Parameters
    ----------
    breath_list : Pandas.DataFrame
        DataFrame containing annotated breaths
    analysis_parameters : dict
        dictionary containing settings for analysis
    animal_metadata : dict
        dict indexed by 'PlyUID' containing animal metadata
    plyuid : string
        serial identifier for data collection event 
        (plethysmography session id)
    auto_criteria : Dict
        Dict with info describing criteria for automated selection of breathing
        for 'calm quiet' breathing
    manual_selection : Dict
        Dict with info describing manually selected start and stop points for
        breathing of interest
    automated_selections_filters : Nested Dict with Lists
        Nested Dict with binary lists describing which breaths are associated
        with passing/failing particular filters
    manual_selections_filters : Nested Dict with Lists
        Nested Dict with binary lists describing which breaths are associated 
        with passing/failing particular filters
    local_logger : instance of logging.logger

    Returns
    -------
    calibration_dict : Nested Dict
        Nested Dict containing calibration parameters to be applied to 
        automated of manually selected breaths

    """
    
    calibration_dict = {}

    if 'Bar_Pres' in animal_metadata[plyuid]:
        calibration_dict['room_baro_mmHg'] = \
            animal_metadata[plyuid]['Bar_Pres']*25.4
        local_logger.info('Barometric Pressure field found [Bar_Pres] ' +
                          'units assumed to be inHg')
    elif 'Room Baro mmHg' in animal_metadata[plyuid]:
        calibration_dict['room_baro_mmHg'] = \
            animal_metadata[plyuid]['Room Baro mmHg']
        local_logger.info('Barometric Pressure field found [Room Baro mmHg] ' +
                          'units assumed to be mmHg')
    else:
        calibration_dict['room_baro_mmHg'] = 760
        local_logger.warning('No Barometric Pressure field found, set at ' +
                             'default of 760 mmHg')

    if 'Calibration_Volume' in animal_metadata[plyuid]:
        calibration_dict['cal_vol_mL'] = \
            animal_metadata[plyuid]['Calibration_Volume'] / 1000
        local_logger.info('Calibration Volume field found ' +
                          '[Calibration_Volume] units assumed to be uL')
    elif 'Cal vol mL' in animal_metadata[plyuid]:
        calibration_dict['cal_vol_mL'] = \
            animal_metadata[plyuid]['Cal vol mL']
        local_logger.info('Calibration Volume field found ' +
                          '[Cal vol mL] units assumed to be mL')
    else:
        calibration_dict['cal_vol_mL'] = 0.02
        local_logger.warning('No Calibration Volume field found, set at ' +
                             'default of 0.02 mL')

    if 'Rotometer_Flowrate' in animal_metadata[plyuid]:
        calibration_dict['Flowrate (SLPM)'] = \
            numpy.interp(
                animal_metadata[plyuid]['Rotometer_Flowrate'],
                [float(i) for i in
                 analysis_parameters[
                     'rotometer_standard_curve_readings'
                     ].split(' ')
                 ],
                [float(i) for i in
                 analysis_parameters[
                     'rotometer_standard_curve_flowrates'
                     ].split(' ')
                 ]
                )
        local_logger.info('Flowrate field found [Rotometer_Flowrate] ' +
                          'units converted to SLPM using std curve')
    elif 'Flowrate (SLPM)' in animal_metadata[plyuid]:
        calibration_dict['Flowrate (SLPM)'] = animal_metadata[plyuid][
            'Flowrate (SLPM)'
            ]
        local_logger.info('Flowrate field found [Flowrate (SLPM)] ' +
                          'units assumed to be SLPM')
    else:
        calibration_dict['Flowrate (SLPM)'] = 0.5
        local_logger.warning('No Flowrate field found, set at ' +
                             'default of 0.5 SLPM')

    if 'Weight' in animal_metadata[plyuid]:
        calibration_dict['Weight'] = animal_metadata[plyuid]['Weight']
        local_logger.info('Body Weight field found [Weight] ' +
                          'units assumed to be g')
    else:
        calibration_dict['Weight'] = 1
        local_logger.warning('No Body Weight field found, set at ' +
                             'default of 20g')

    if auto_criteria is not None or manual_selection is not None:
        calibration_dict['Cal_Segs'] = {}

    if auto_criteria is not None:
        calibration_dict['Auto_Cal_Match'] = {}
        for row in auto_criteria.index:
            key_and_seg = '{}{}'.format(
                auto_criteria.iloc[row]['Key'],
                auto_criteria.iloc[row]['Alias']
                )
            if key_and_seg not in automated_selections_filters:
                continue
            calibration_dict['Auto_Cal_Match'][
                '{}{}'.format(
                    auto_criteria.iloc[row]['Key'],
                    auto_criteria.iloc[row]['Alias']
                    )
                ] = {
                    'Key': auto_criteria.iloc[row]['Key'],
                    'Alias': auto_criteria.iloc[row]['Alias'],
                    'Cal Seg': auto_criteria.iloc[row]['Cal Seg']
                    }
        for row in auto_criteria[auto_criteria['AUTO_IND_CAL'] == 1].index:
            key_and_seg = '{}{}'.format(
                auto_criteria.iloc[row]['Key'],
                auto_criteria.iloc[row]['Alias']
                )
            if key_and_seg not in automated_selections_filters:
                continue
            calibration_dict['Cal_Segs'][
                auto_criteria.iloc[row]['Alias']] = {
                    'Key or segment': auto_criteria.iloc[row]['Key'],
                    'Alias': auto_criteria.iloc[row]['Alias'],
                    'Auto_or_Man': 'Auto',
                    'o2': breath_list[
                        (automated_selections_filters[key_and_seg]\
                         ['biggest_block'] == 1) &
                        (automated_selections_filters[key_and_seg]
                         ['selection_filter'] == 1)
                        ]['calibrated_o2'].mean(),
                    'co2': breath_list[
                        (automated_selections_filters[key_and_seg]\
                         ['biggest_block'] == 1) &
                        (automated_selections_filters[key_and_seg]
                         ['selection_filter'] == 1)
                        ]['calibrated_co2'].mean(),
                    'tv': breath_list[
                        (automated_selections_filters[key_and_seg]\
                         ['biggest_block'] == 1) &
                        (automated_selections_filters[key_and_seg]
                         ['selection_filter'] == 1)
                        ]['iTV'].mean()
                    }

    if manual_selection is not None:
        calibration_dict['Man_Cal_Match'] = {}
        for row in manual_selection[
                manual_selection['PLYUID'].astype(str) == plyuid
                ].index:
            calibration_dict['Man_Cal_Match'][
                manual_selection.iloc[row]['Alias']
                ] = {
                    'segment': manual_selection.iloc[row]['segment'],
                    'Alias': manual_selection.iloc[row]['Alias'],
                    'Cal Seg': manual_selection.iloc[row]['Cal Seg']
                    }
        for row in manual_selection[
                (manual_selection['PLYUID'].astype(str) == plyuid) &
                (manual_selection['MAN_IND_CAL'] == 1)
                ].index:
            if manual_selection.iloc[row]['Alias'] in \
                    calibration_dict['Cal_Segs']:
                local_logger.warning(
                    'Manual Selection defines a calibration segment ' +
                    'that was already added, settings may be overwritten'
                    )
            calibration_dict['Cal_Segs'][
                manual_selection.iloc[row]['Alias']] = {
                    'Key or segment': manual_selection.iloc[row]['segment'],
                    'Alias': manual_selection.iloc[row]['Alias'],
                    'Auto_or_Man': 'Man',
                    'o2': breath_list[
                        manual_selections_filters['{}'.format(
                            manual_selection.iloc[row]['Alias']
                            )]['selection_filter'] == 1
                        ]['calibrated_o2'].mean(),
                    'co2': breath_list[
                        manual_selections_filters['{}'.format(
                            manual_selection.iloc[row]['Alias']
                            )]['selection_filter'] == 1
                        ]['calibrated_co2'].mean(),
                    'tv': breath_list[
                        manual_selections_filters['{}'.format(
                            manual_selection.iloc[row]['Alias']
                            )]['selection_filter'] == 1
                        ]['iTV'].mean()
                    }

    return calibration_dict


def getvaporpressure(temperature):
    """
    Calculates vapor pressure of water at a given temperature.
    
    Parameters
    ----------
    temperature : Float
        Temperature in degrees Celsius

    Returns
    -------
    Float
        Vapor Pressure of Water

    """
    return \
        (
            1.142 + (0.8017 * temperature) -
            (0.012 * temperature ** 2) + (0.0006468 * temperature ** 3)
            )

# def getK(C):
#     """
#     Converts Centigrade to Kelvin.
    
#     Parameters
#     ----------
#     C : Float
#         Temperature in degrees Celsius

#     Returns
#     -------
#     Float
#         Temperature in Kelvin

#     """
#     return (C+273.15)


def get_BT_TV_K(tv, calv, act_calv, t_body, t_chamb, room_pressure):
    """
    Calculates Bartlett and Tenney based tidal volume values.
    
    Parameters
    ----------
    tv : Float
        uncorrected tidal volume (V)
    calv : Float
        uncorrected tidal volume from calibration period (V)
    act_calv : Float
        nominal calibration volume (mL)
    t_body : Float
        Body Temperature (Celsius)
    t_chamb : Float
        Chamber Temperature (Celsius)
    room_pressure : Float
        Room Barometric Pressure (mmHg)

    Returns
    -------
    Float
        corrected tidal volume (mL)

    """
    # Convert temperatures to Fahrenheit
    t_kelvin_body = pytemperature.c2k(t_body)
    t_kelvin_chamb = pytemperature.c2k(t_chamb)

    # Calculate corrected tidal volume
    corrected_tv = (tv / calv) * (act_calv) * \
        (t_kelvin_body * (room_pressure - getvaporpressure(t_chamb))) / \
        ((t_kelvin_body * (room_pressure - getvaporpressure(t_chamb))) -
         t_kelvin_chamb * (room_pressure - getvaporpressure(t_body)))

    return corrected_tv


def get_pneumo_TV(tv: float, calv: float, act_calv: float) -> float:
    """
    Calculates pneumotacograph based tidal volume.
    
    Parameters
    ----------
    tv : float
        uncorrected tidal volume (V)
    calv : float
        uncorrected tidal volume from calibration period (V)
    act_calv : float
        nominal calibration volume (mL)

    Returns
    -------
    float
        corrected tidal volume (mL)
    """
    return tv / calv * act_calv


def get_VO2(o2_in, o2_out, flowrate):
    """
    Calculates Volume of O2 (mL) consumed per minute.

    Parameters
    ----------
    o2_in : Float
        O2 concentration (%) into the chamber
    o2_out : Float
        O2 concentration (%) out of the chamber
    flowrate : Float
        flowrate (mL) of air through the chamber

    Returns
    -------
    Float
        Volume of O2 (mL) consumed per minute.

    """
    
    return (flowrate*1000*o2_in/100-flowrate*1000*o2_out/100)


def get_VCO2(co2_in, co2_out, flowrate):
    """
    Calculates Volume of CO2 (mL) produced per minute.

    Parameters
    ----------
    co2_in : Float
        CO2 concentration (%) into the chamber
    co2_out : Float
        CO2 concentration (%) out of the chamber
    flowrate : Float
        flowrate (mL) of air through the chamber

    Returns
    -------
    Float
        Volume of CO2 (mL) produced per minute.    
        
    """
    
    return (flowrate*1000*co2_out/100-flowrate*1000*co2_in/100)


def apply_volumetric_and_respiratory_calculations(
    breath_list: pandas.DataFrame,
    calibration_parameters: dict,
    automated_selections_filters: dict,
    manual_selections_filters: dict,
    analysis_parameters: dict,
    local_logger: logging.Logger
) -> pandas.DataFrame:
    local_logger.info('Applying volumetric and respiratory calculations')
    basic_volume_list = populate_baseline_values_for_calibration_calculations(
        breath_list, 
        calibration_parameters, 
        automated_selections_filters, 
        manual_selections_filters, 
        local_logger
    )
    if analysis_parameters.get('Pneumo_Mode') == '1':
        enhanced_volume_list = calculate_calibrated_pneumo_volume_and_respiration(
            breath_list,
            basic_volume_list,
            calibration_parameters
        )
    else:
        enhanced_volume_list = calculate_calibrated_volume_and_respiration(
            breath_list,
            basic_volume_list,
            calibration_parameters
        )
    return enhanced_volume_list

def calculate_calibrated_pneumo_volume_and_respiration(
    breath_list: pandas.DataFrame,
    volume_list: pandas.DataFrame,
    calibration_parameters: dict
) -> pandas.DataFrame:
    """
    Calculates tidal volumes and flows from a pneumotacograph signal as well
    as VO2 and VCO2 related measures.

    Parameters
    ----------
    breath_list : Pandas.DataFrame
        DataFrame containing annotated breaths
    volume_list : Pandas.DataFrame
        DataFrame related to breath_list, but with baseline volume data placed
        for later use in calculations
    calibration_parameters : Dict
        Dict containing parameters needed for calculating calibrated volumes

    Returns
    -------
    enhanced_volume_list : Pandas.DataFrame
        DataFrame describing calibrated volume and respiratory data

    """

    enhanced_volume_list = volume_list.copy()
    
    enhanced_volume_list['calibrated_TV'] = get_pneumo_TV(
        breath_list['iTV'],
        volume_list['base_tv'],
        calibration_parameters['cal_vol_mL'],
    )

    enhanced_volume_list['calibrated_eTV'] = get_pneumo_TV(
        breath_list['eTV'],
        volume_list['base_tv'],
        calibration_parameters['cal_vol_mL'],
    )

    enhanced_volume_list['calibrated_PIF'] = get_pneumo_TV(
        breath_list['PIF'],
        volume_list['base_tv'],
        calibration_parameters['cal_vol_mL'],
    )

    enhanced_volume_list['calibrated_PEF'] = get_pneumo_TV(
        breath_list['PEF'],
        volume_list['base_tv'],
        calibration_parameters['cal_vol_mL'],
    )

    enhanced_volume_list['VO2'] = get_VO2(
        volume_list['base_o2'],
        breath_list['calibrated_o2'],
        calibration_parameters['Flowrate (SLPM)'],
    )

    enhanced_volume_list['VCO2'] = get_VCO2(
        volume_list['base_co2'],
        breath_list['calibrated_co2'],
        calibration_parameters['Flowrate (SLPM)'],
    )
    
    enhanced_volume_list['VE'] = breath_list['BPM'] * \
        enhanced_volume_list['calibrated_TV']
    enhanced_volume_list['VE_per_VO2'] = enhanced_volume_list['VE'] / \
        enhanced_volume_list['VO2']
    enhanced_volume_list['RER'] = enhanced_volume_list['VCO2'] / \
        enhanced_volume_list['VO2']
    enhanced_volume_list['calibrated_TV_per_gram'] = enhanced_volume_list[
        'calibrated_TV'] / calibration_parameters['Weight']
    enhanced_volume_list['VE_per_gram'] = enhanced_volume_list[
        'VE'] / calibration_parameters['Weight']
    enhanced_volume_list['VO2_per_gram'] = enhanced_volume_list[
        'VO2'] / calibration_parameters['Weight']
    enhanced_volume_list['VCO2_per_gram'] = enhanced_volume_list[
        'VCO2'] / calibration_parameters['Weight']
    enhanced_volume_list['TT_per_TV'] = breath_list['TT'] / \
        enhanced_volume_list['calibrated_TV']
    enhanced_volume_list['TT_per_TVpg'] = breath_list['TT'] / \
        enhanced_volume_list['calibrated_TV_per_gram']
    enhanced_volume_list['O2_per_Air'] = enhanced_volume_list[
        'VO2'] * breath_list['TT'] / 60 / enhanced_volume_list['calibrated_TV']
    
    enhanced_volume_list['calibrated_VT_per_TI'] = enhanced_volume_list[
        'calibrated_TV'
    ] / breath_list['TI']
    enhanced_volume_list['calibrated_VT_per_TE'] = enhanced_volume_list[
        'calibrated_eTV'
    ] / breath_list['TE']
    enhanced_volume_list['uncalibrated_VT_per_TI'] = breath_list[
        'iTV'
    ] / breath_list['TI']
    enhanced_volume_list['uncalibrated_VT_per_TE'] = breath_list[
        'eTV'
    ] / breath_list['TE']

    return enhanced_volume_list

def populate_baseline_values_for_calibration_calculations(
        breath_list,
        calibration_parameters,
        automated_selections_filters,
        manual_selections_filters,
        local_logger
        ):
    """
    Creates a dataframe that contains columns for baseline TV, O2, and CO2 
    measurements for use in creating compensated/calibrated.

    Parameters
    ----------
    breath_list : Pandas.DataFrame
        DataFrame containing annotated breaths
    calibration_parameters : Dict
        Dict containing parameters needed for calculating calibrated volumes
    automated_selections_filters : Nested Dict with Lists
        Nested Dict with binary lists describing which breaths are associated
        with passing/failing particular filters
    manual_selections_filters : Nested Dict with Lists
        Nested Dict with binary lists describing which breaths are associated 
        with passing/failing particular filters
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    volume_list : Pandas.DataFrame
        DataFrame related to breath_list, but with baseline volume data placed
        for later use in calculations

    """
    
    volume_list = pandas.DataFrame(breath_list['ts_inhale'].copy())
    
    volume_list['base_selection'] = numpy.nan
    volume_list['base_o2'] = numpy.nan
    volume_list['base_co2'] = numpy.nan
    volume_list['base_tv'] = numpy.nan

    if 'Auto_Cal_Match' in calibration_parameters:
        if len(calibration_parameters['Auto_Cal_Match']) > 0:
            for key_and_alias in automated_selections_filters:
                if str(
                        calibration_parameters[
                            'Auto_Cal_Match'
                            ][
                                key_and_alias
                                ][
                                    'Cal Seg'
                                    ]
                        ).lower() in ['nan','na']:
                    continue
                
                volume_list.loc[
                    (automated_selections_filters[key_and_alias]
                         ['biggest_block'] == 1) &
                    (automated_selections_filters[key_and_alias]
                         ['selection_filter'] == 1),
                    'base_selection'] = calibration_parameters[
                        'Auto_Cal_Match'
                    ][key_and_alias]['Cal Seg']
                volume_list.loc[
                    (automated_selections_filters[key_and_alias]
                         ['biggest_block'] == 1) &
                    (automated_selections_filters[key_and_alias]
                         ['selection_filter'] == 1),
                    'base_o2'] = calibration_parameters[
                        'Cal_Segs'][calibration_parameters['Auto_Cal_Match']
                                    [key_and_alias]['Cal Seg']]['o2']
                volume_list.loc[
                    (automated_selections_filters[key_and_alias]
                         ['biggest_block'] == 1) &
                    (automated_selections_filters[key_and_alias]
                         ['selection_filter'] == 1),
                    'base_co2'] = calibration_parameters[
                        'Cal_Segs'][calibration_parameters['Auto_Cal_Match']
                                    [key_and_alias]['Cal Seg']]['co2']
                volume_list.loc[
                    (automated_selections_filters[key_and_alias]
                         ['biggest_block'] == 1) &
                    (automated_selections_filters[key_and_alias]
                         ['selection_filter'] == 1),
                    'base_tv'] = calibration_parameters[
                        'Cal_Segs'][calibration_parameters['Auto_Cal_Match']
                                    [key_and_alias]['Cal Seg']]['tv']

    if 'Man_Cal_Match' in calibration_parameters:
        if len(calibration_parameters['Man_Cal_Match']) > 0:
            for alias in manual_selections_filters:
                if str(
                        calibration_parameters[
                            'Man_Cal_Match'
                            ][
                                alias
                                ][
                                    'Cal Seg'
                                    ]
                        ).lower() in ['nan','na']:
                    continue
                
                volume_list.loc[
                    manual_selections_filters[alias]
                         ['selection_filter'] == 1,
                         'base_selection'] = calibration_parameters[
                             'Man_Cal_Match'
                             ][alias]['Cal Seg']
                volume_list.loc[
                    manual_selections_filters[alias]
                        ['selection_filter'] == 1,
                        'base_o2'] = calibration_parameters[
                            'Cal_Segs'][calibration_parameters['Man_Cal_Match']
                                        [alias]['Cal Seg']]['o2']
                volume_list.loc[
                    manual_selections_filters[alias]
                        ['selection_filter'] == 1,
                        'base_co2'] = calibration_parameters[
                            'Cal_Segs'][calibration_parameters['Man_Cal_Match']
                                        [alias]['Cal Seg']]['co2']
                volume_list.loc[
                    manual_selections_filters[alias]
                        ['selection_filter'] == 1,
                        'base_tv'] = calibration_parameters[
                            'Cal_Segs'][calibration_parameters['Man_Cal_Match']
                                        [alias]['Cal Seg']]['tv']

    return volume_list



def calculate_calibrated_volume_and_respiration(
        breath_list,
        volume_list,
        calibration_parameters
        ):
    """
    Calculates Bartlett and Tenney compensated tidal volumes and flows as well
    as VO2 and VCO2 related measures.

    Parameters
    ----------
    breath_list : Pandas.DataFrame
        DataFrame containing annotated breaths
    volume_list : Pandas.DataFrame
        DataFrame related to breath_list, but with baseline volume data placed
        for later use in calculations
    calibration_parameters : Dict
        Dict containing parameters needed for calculating calibrated volumes

    Returns
    -------
    enhanced_volume_list : Pandas.DataFrame
        DataFrame describing calibrated volume and respiratory data

    """
    
    enhanced_volume_list = volume_list.copy()
    enhanced_volume_list['calibrated_TV'] = get_BT_TV_K(
        breath_list['iTV'],
        volume_list['base_tv'],
        calibration_parameters['cal_vol_mL'],
        breath_list['Body_Temperature'],
        breath_list['corrected_temp'],
        calibration_parameters['room_baro_mmHg']
        )

    enhanced_volume_list['calibrated_eTV'] = get_BT_TV_K(
        breath_list['eTV'],
        volume_list['base_tv'],
        calibration_parameters['cal_vol_mL'],
        breath_list['Body_Temperature'],
        breath_list['corrected_temp'],
        calibration_parameters['room_baro_mmHg']
        )

    enhanced_volume_list['calibrated_PIF'] = get_BT_TV_K(
        breath_list['PIF'],
        volume_list['base_tv'],
        calibration_parameters['cal_vol_mL'],
        breath_list['Body_Temperature'],
        breath_list['corrected_temp'],
        calibration_parameters['room_baro_mmHg']
        )

    enhanced_volume_list['calibrated_PEF'] = get_BT_TV_K(
        breath_list['PEF'],
        volume_list['base_tv'],
        calibration_parameters['cal_vol_mL'],
        breath_list['Body_Temperature'],
        breath_list['corrected_temp'],
        calibration_parameters['room_baro_mmHg']
        )

    enhanced_volume_list['VO2'] = get_VO2(
        volume_list['base_o2'],
        breath_list['calibrated_o2'],
        calibration_parameters['Flowrate (SLPM)']
        )

    enhanced_volume_list['VCO2'] = get_VCO2(
        volume_list['base_co2'],
        breath_list['calibrated_co2'],
        calibration_parameters['Flowrate (SLPM)']
        )
    
    enhanced_volume_list['VE'] = breath_list['BPM'] * \
        enhanced_volume_list['calibrated_TV']
    enhanced_volume_list['VE_per_VO2'] = enhanced_volume_list['VE'] / \
        enhanced_volume_list['VO2']
    enhanced_volume_list['RER'] = enhanced_volume_list['VCO2'] / \
        enhanced_volume_list['VO2']
    enhanced_volume_list['calibrated_TV_per_gram'] = enhanced_volume_list[
        'calibrated_TV'] / calibration_parameters['Weight']
    enhanced_volume_list['VE_per_gram'] = enhanced_volume_list[
        'VE'] / calibration_parameters['Weight']
    enhanced_volume_list['VO2_per_gram'] = enhanced_volume_list[
        'VO2'] / calibration_parameters['Weight']
    enhanced_volume_list['VCO2_per_gram'] = enhanced_volume_list[
        'VCO2'] / calibration_parameters['Weight']
    enhanced_volume_list['TT_per_TV'] = breath_list['TT'] / \
        enhanced_volume_list['calibrated_TV']
    enhanced_volume_list['TT_per_TVpg'] = breath_list['TT'] / \
        enhanced_volume_list['calibrated_TV_per_gram']
    enhanced_volume_list['O2_per_Air'] = enhanced_volume_list[
        'VO2'] * breath_list['TT'] / 60 / enhanced_volume_list['calibrated_TV']

    enhanced_volume_list['calibrated_VT_per_TI'] = enhanced_volume_list[
        'calibrated_TV'
    ] / breath_list['TI']
    enhanced_volume_list['calibrated_VT_per_TE'] = enhanced_volume_list[
        'calibrated_eTV'
    ] / breath_list['TE']
    enhanced_volume_list['uncalibrated_VT_per_TI'] = breath_list[
        'iTV'
    ] / breath_list['TI']
    enhanced_volume_list['uncalibrated_VT_per_TE'] = breath_list[
        'eTV'
    ] / breath_list['TE']

    return enhanced_volume_list



def create_output(
        breath_list,
        analysis_parameters,
        animal_metadata,
        auto_criteria,
        manual_selection,
        plyuid_or_ruid,
        automated_selections_filters,
        manual_selections_filters,
        column_dictionary,
        local_logger
        ):                
    """
    Creates JSON output for use with STAGG and optionally a csv file with all
    identified breaths.

    Parameters
    ----------
    breath_list : Pandas.DataFrame
        DataFrame containing annotated breaths
    analysis_parameters : Dict
        Dict of settings to use for analysis
    animal_metadata : Dict
        Dict of metadata information for an animal
    auto_criteria : Dict
        Dict with info describing criteria for automated selection of breathing
        for 'calm quiet' breathing
    manual_selection : Dict
        Dict with info describing manually selected start and stop points for
        breathing of interest
    plyuid : string
        serial identifier for data collection event 
        (plethysmography session id)
    automated_selections_filters : Nested Dict with Lists
        Nested Dict with binary lists describing which breaths are associated
        with passing/failing particular filters
    manual_selections_filters : Nested Dict with Lists
        Nested Dict with binary lists describing which breaths are associated 
        with passing/failing particular filters
        DESCRIPTION.
    column_dictionary : Dict
        Dict that indicates translation of current 'column names' for output
        to match the current preferred style
    local_logger : instance of logging.logger (optional)

    Returns
    -------
    output_list : Pandas DataFrame
        STAGG ready output, contains a breathlist and associated parameters 
        and metadata for subsequent statistical analysis

    """
    
    output_list = breath_list.copy()
    
    output_list['Exp_Condition'] = ""
    output_list['Auto_Condition'] = ""
    output_list['Breath_Inclusion_Filter'] = 0
    output_list['AUTO_Inclusion_Filter'] = 0
    output_list['MAN_Inclusion_Filter'] = 0
    output_list['Man_Condition'] = ""
    output_list['Auto_Block_Id'] = ""
    output_list['Auto_Selection_Id'] = ""
    output_list['Man_Selection_Id'] = ""
    output_list['AUTO_IND_INCLUDE'] = 0
    output_list['MAN_IND_INCLUDE'] = 0
    
    if auto_criteria is not None:
        
        auto_criteria['key_and_alias'] = \
            auto_criteria.copy()['Key'].astype(str) + \
            auto_criteria.copy()['Alias'].astype(str)
        
                            
        for row in auto_criteria.index:
            if auto_criteria.iloc[row]['key_and_alias'] not in \
                    automated_selections_filters:
                continue
            for c in auto_criteria.columns:
                if 'AUTO_IND_' == c[0:9]:
                        
                    if c not in output_list.columns:
                        output_list[c] = ""
                        
                    if auto_criteria.groupby(['Alias',c]).ngroups > \
                            auto_criteria.groupby(['Alias']).ngroups:
                        local_logger.warning(
                            'Shared Aliases in autocriteria with variant ' +\
                            'entries for {}'.format(c)
                            )
        
                    output_list.loc[automated_selections_filters[
                        auto_criteria.iloc[row]['key_and_alias']][
                            'biggest_block']==1,c] = auto_criteria.iloc[
                                row][c]
                                
            if 'AUTO_Condition_{}'.format(
                    auto_criteria.iloc[row]['Alias']
                    ) not in output_list.columns:
                output_list[
                    'AUTO_Condition_{}'.format(
                        auto_criteria.iloc[row]['Alias']
                        )
                    ] = 0
                
            output_list.loc[automated_selections_filters[
                auto_criteria.iloc[row]['key_and_alias']][
                    'biggest_block']==1,
                    'AUTO_Condition_{}'.format(
                        auto_criteria.iloc[row]['Alias']
                        )
                    ] = 1
                    
            output_list.loc[automated_selections_filters[
                auto_criteria.iloc[row]['key_and_alias']][
                    'biggest_block']==1,
                    'Auto_Condition'] = \
                        auto_criteria.iloc[row]['Alias']
                        
            output_list.loc[automated_selections_filters[
                auto_criteria.iloc[row]['key_and_alias']][
                    'biggest_block']==1,
                    'Auto_Block_Id'] = \
                        automated_selections_filters[
                            auto_criteria.iloc[row]['key_and_alias']
                            ]['block_id'][automated_selections_filters[
                                auto_criteria.iloc[row]['key_and_alias']][
                                    'biggest_block']==1]
                                
            
            output_list.loc[automated_selections_filters[
                auto_criteria.iloc[row]['key_and_alias']][
                    'selection_filter']==1,
                    'Breath_Inclusion_Filter'] = 1
                    
            output_list.loc[automated_selections_filters[
                auto_criteria.iloc[row]['key_and_alias']][
                    'selection_filter']==1,
                    'AUTO_Inclusion_Filter'] = 1
            
            output_list.loc[automated_selections_filters[
                auto_criteria.iloc[row]['key_and_alias']][
                    'selection_filter']==1,
                    'Auto_Selection_Id'] = \
                        automated_selections_filters[
                            auto_criteria.iloc[row]['key_and_alias']
                            ]['selection_id'][automated_selections_filters[
                                auto_criteria.iloc[row]['key_and_alias']][
                                    'selection_filter']==1]
                    
            
                        
    if manual_selection is not None:
        for alias in manual_selections_filters:
            for c in manual_selection.columns:
                if 'MAN_IND_' == c[0:8]:
                    
                    if c not in output_list.columns:
                        breath_list[c] = ""
                        
                    if manual_selection[
                            manual_selection['PLYUID'].astype(str)==\
                            plyuid_or_ruid \
                            ].groupby(['Alias',c]).ngroups > \
                            manual_selection[
                                manual_selection['PLYUID'].astype(str)==\
                                    plyuid_or_ruid \
                                    ].groupby(['Alias']).ngroups:
                        local_logger.warning(
                            'Shared Aliases in manual select with variant '+\
                            'entries for {}'.format(c)
                            )
                    output_list.loc[manual_selections_filters[alias][
                        'selection_filter'
                        ]==1,c] = manual_selection[
                            (manual_selection['Alias']==alias)&
                            (manual_selection['PLYUID'].astype(str)==\
                            plyuid_or_ruid)
                            ][c].values[0]
            output_list['MAN_Condition_{}'.format(alias)] = 0
            output_list.loc[manual_selections_filters[alias]\
                            ['selection_filter']==1,'MAN_Condition_{}'.format(
                                alias
                                )
                            ] = 1
            output_list.loc[:,'Man_Selection_Id'] = \
                manual_selections_filters[alias]['selection_id']
                
            output_list.loc[manual_selections_filters[alias]\
                            ['selection_filter']==1,
                            'Man_Condition'] = alias
            
            output_list.loc[manual_selections_filters[alias]\
                            ['selection_filter']==1,
                            'Breath_Inclusion_Filter'] = 1
            
            output_list.loc[manual_selections_filters[alias]\
                            ['selection_filter']==1,
                            'MAN_Inclusion_Filter'] = 1

    for m in animal_metadata[plyuid_or_ruid]:
        output_list[m] = animal_metadata[plyuid_or_ruid][m]
    
    if analysis_parameters.get('Pneumo_Mode') != '1':
        output_list['Mouse_And_Session_ID'] = \
            output_list['MUID'].copy().astype(str) + '_' + \
            output_list['PlyUID'].copy().astype(str)
    
    for p in analysis_parameters:
        output_list[p] = analysis_parameters[p]
    
    reversed_column_dictionary = {}
    for k in column_dictionary:
        reversed_column_dictionary[column_dictionary[k]] = k
    
    output_list.rename(columns=reversed_column_dictionary, inplace = True)
    column_list = []
    for k in column_dictionary:
        column_list.append(k)
        if k in [
                'Mouse_And_Session_ID',
                'Breath_Inclusion_Filter',
                'Exp_Condition',
                'Auto_Condition',
                'Man_Condition',
                'Breath Number',
                'Auto_Block_Id',
                'Auto_Selection_Id',
                'Man_Selection_Id',
                'Timestamp_Inspiration',
                'Timestamp_Expiration',
                ]:
            continue
        output_list['Irreg_Score_{}'.format(k)] = calculate_irreg_score(
            output_list[k]
            )
        column_list.append('Irreg_Score_{}'.format(k))
    
    output_list.loc[:,'Exp_Condition'] = \
        output_list['Auto_Condition'].astype(str) + \
        output_list['Man_Condition'].astype(str)
    
    return output_list


def start_logger(logger_output_path, 
                  level: int = logging.DEBUG, 
                  level_console: int = None,
                  level_file: int = None,
                  logname=__name__
                  ):
    
    if not level_console: level_console = level
    if not level_file: level_file = level
    level = min(level_console,level_file)
    
    # create logger
    logger = logging.getLogger(logname)
    logger.setLevel(level)

    # create file and console handlers to receive logging
    console_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(logger_output_path)
    
    console_handler.setLevel(level_console)
    file_handler.setLevel(level_file)

    # create format for log and apply to handlers
    log_format = logging.Formatter(
            '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
            )
    console_handler.setFormatter(log_format)
    file_handler.setFormatter(log_format)

    # add handlers to the logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    # log initial inputs
    logger.info('Logger Started')
    return logger


#%% functions to run SASSI

def run_SASSI_from_data(
        file : str = None,
        signal_data : pandas.DataFrame = None,
        output_path : str = None,
        animal_metadata : dict = None,
        manual_selection : pandas.DataFrame = None,
        auto_criteria : pandas.DataFrame = None,
        channel_match : pandas.DataFrame = None,
        analysis_parameters : dict = None,
        loglevel : int = logging.DEBUG,
        status_queue = None,
        logger : logging.Logger = None,
        pid : int = 0
        ):
    
    egregious=[]
    
    if status_queue: status_queue.put_nowait(
            f'{pid} : working on {os.path.basename(file)}'
            )
    
    if logger is None:
        log_path = os.path.join(output_path, os.path.basename(file) + '.log')
        logger = start_logger(log_path, level=logging.DEBUG, logname=os.path.basename(file))
    
    logger.info('Function start.')
    
    if logger is None:
        logger = start_logger(
            os.path.join(output_path,os.path.basename(file)+'.log'),
            level = loglevel,
            logname = os.path.basename(file)
            )
    logger.info('loading data')
    try:    
        # load signal data if not provided as a dataframe
        if signal_data is None:
            signal_data = load_signal_data(
                file, 
                analysis_parameters = analysis_parameters,
                channel_match = channel_match,
                output_path = output_path,
                local_logger = logger
                )
        
    
        column_dictionary = {
            'Mouse_And_Session_ID':'Mouse_And_Session_ID',
            'Breath_Inclusion_Filter':'Breath_Inclusion_Filter',
            'Auto_Condition':'Auto_Condition',
            'Man_Condition':'Man_Condition',
            'Exp_Condition':'Exp_Condition',
            'Auto_Block_Id':'Auto_Block_Id',
            'Auto_Selection_Id':'Auto_Selection_Id',
            'Man_Selection_Id':'Man_Selection_Id',
            'Breath Number':'il_inhale',
            'Timestamp_Inspiration':'ts_inhale',
            'Timestamp_Expiration':'ts_exhale',
            'Inspiratory_Duration':'TI',
            'Expiratory_Duration':'TE',
            'Tidal_Volume_uncorrected':'iTV',
            'Tidal_Volume_exhale_uncorrected':'eTV',
            'VT__Tidal_Volume_corrected':'calibrated_TV',
            'VT_per_TI__Tidal_Volume_corr_per_Ins_Time':'calibrated_VT_per_TI',
            'VT_per_TE__Tidal_Volume_corr_per_Exp_Time':'calibrated_VT_per_TE',
            'uncorrected_VT_per_TI':'uncalibrated_VT_per_TI',
            'uncorrected_VT_per_TE':'uncalibrated_VT_per_TE',
            'Tidal_Volume_exhale_corrected':'calibrated_eTV',
            'VTpg__Tidal_Volume_per_gram_corrected':'calibrated_TV_per_gram',
            'Peak_Inspiratory_Flow':'PIF',
            'Peak_Inspiratory_Flow_corrected':'calibrated_PIF',
            'Peak_Expiratory_Flow':'PEF',
            'Peak_Expiratory_Flow_corrected':'calibrated_PEF',
            'Breath_Cycle_Duration':'TT',
            'VE__Ventilation':'VE',
            'VEpg__Ventilation_per_gram':'VE_per_gram',
            'VO2':'VO2',
            'VO2pg':'VO2_per_gram',
            'VCO2':'VCO2',
            'VCO2pg':'VCO2_per_gram',
            'VEVO2':'VE_per_VO2',
            'VF':'BPM',
            'TT_per_TV':'TT_per_TV',
            'TT_per_TVpg':'TT_per_TVpg',
            'O2_per_Air__VO2_x_TT_per_TV_':'O2_per_Air',
            'Apnea':'apnea',
            'Sigh':'sigh',
            'DVTV':'DVTV',
            'per500':'per_X_calculation',
            'mav':'mov_avg_vol',
            'O2_concentration':'calibrated_o2',
            'CO2_concentration':'calibrated_co2',
            'Chamber_Temperature':'corrected_temp',
            'O2_uncalibrated':'o2',
            'CO2_uncalibrated':'co2',
            'Chamber_Temp_uncalibrated':'temp',
            'Body_Temperature_Linear':'Body_Temperature',
            }
        

        #% confirm animal is present in metadata
        # extract MUID and PLYUID
        
        
        if str(analysis_parameters.get('Pneumo_Mode')) == '1':
            
            # pneumo mode sample id extract


            RUID = extract_ruid(file)


                                    
            PLYUID_or_RUID = RUID

        else:
            
            # normal mode sample id extract

        
            MUID, PLYUID = extract_muid_plyuid(
                file,
                animal_metadata,
                local_logger=logger
                )


            PLYUID_or_RUID = PLYUID

    
        
    
        # apply basic voltage corrections to gas and temperature values
        signal_data = apply_voltage_corrections(
            signal_data,
            analysis_parameters,
            [
                ('o2', 'corrected_o2', 'O2_calibration_factor'),
                ('co2', 'corrected_co2', 'CO2_calibration_factor'),
                (
                    'temp',
                    'corrected_temp',
                    'temperature_calibration_factor'
                    )
                ],
            logger
            )
        
    
        
        if analysis_parameters.get('Pneumo_Mode') == '1' or \
                analysis_parameters.get('signal_format') in ['wa']:
            purge_columns = []
            if 'corrected_temp' not in signal_data:
                signal_data['corrected_temp'] = -1
                purge_columns.append('corrected_temp')
                logger.warning('corrected temperature data not in signal data')
            if 'temp' not in signal_data:
                signal_data['temp'] = -1
                purge_columns.append('temp')
                logger.warning('temperature data not in signal data')
                
            if 'corrected_co2' not in signal_data:
                signal_data['corrected_co2'] = \
                  analysis_parameters.get('calibrated_co2_default') if\
                  analysis_parameters.get('calibrated_co2_default') \
                  else -1
                signal_data['co2'] = \
                    signal_data['corrected_co2']
                purge_columns.append('calibrated_co2_default')
                logger.warning('CO2 data not in signal data')
                
            if 'corrected_o2' not in signal_data:
                signal_data['corrected_o2'] = \
                  analysis_parameters.get('calibrated_o2_default') if\
                  analysis_parameters.get('calibrated_o2_default') \
                  else -1
                signal_data['o2'] = \
                    signal_data['corrected_o2']
                purge_columns.append('calibrated_o2_default')
                logger.warning('O2 data not in signal data')
    
        # apply smoothing to gas signals if indicated in settings
        signal_data.loc[:,'corrected_o2'], \
        signal_data.loc[:,'corrected_co2'] = smooth_gas_signals(
            signal_data,
            analysis_parameters,
            logger = logger
            )
 
        
        signal_data.loc[:, 'corrected_temp'] = repair_temperature(
            signal_data,
            PLYUID_or_RUID,
            animal_metadata,
            analysis_parameters,
            local_logger=logger
            )
    
        # collect time stamps
        timestamp_dict = dict(
            zip(
                signal_data.ts[signal_data.comment.dropna().index],
                signal_data.comment.dropna()
                )
            )
    

    
        # apply linear estimation of body temperature
        signal_data['Body_Temperature'] = calculate_body_temperature(
            signal_data,
            PLYUID_or_RUID,
            timestamp_dict,
            animal_metadata,
            manual_selection,
            auto_criteria,
            analysis_parameters,
            logger
            )
    
        # differentiate breathing signal to 'derived flow'
        if 'flow' not in signal_data:
            signal_data['flow'] = signal_data['vol'].diff().fillna(1.0)
            logger.info('calculating flow from vol data')
        else:
            logger.info('flow data present in signal data')
    
        if 'vol' not in signal_data:
            signal_data['vol'] = signal_data['flow'].cumsum()
            if '1' in str(
                    analysis_parameters['apply_smoothing_filter']
                    ) or 'f' in str(
                        analysis_parameters['apply_smoothing_filter']
                        ): 
                signal_data['vol'] = apply_smoothing_filter(
                    signal_data,
                    'vol'
                    )
            logger.info('calculating vol data from flow data')
        else:                
            logger.info('vol data present in signal data')
    
        # calculate sampling rate (expected to be 1000Hz)
        sampleHz = round(1/(signal_data['ts'][2]-signal_data['ts'][1]))
    
        # calculate moving average of vol signal to use for QC
        # (1sec centered mov avg)
            
        signal_data['mov_avg_vol'] = signal_data['vol'].rolling(
            int(sampleHz), center=True
            ).mean()
    
        # smooth data with high and low pass filter if option selected
        if '1' in str(
                analysis_parameters['apply_smoothing_filter']
                ) or 'f' in str(
                    analysis_parameters['apply_smoothing_filter']
                    ):
            signal_data['original_flow'] = signal_data['flow'].copy()
            signal_data['flow'] = apply_smoothing_filter(
                signal_data,
                'flow'
                )    
        
        
        # call breaths
        if status_queue: status_queue.put_nowait(f'{pid}: Calling Breaths')
        breath_list = breath_caller(
            signal_data,
            analysis_parameters,
            logger
            )
    
        if 'ecg' in signal_data:
            logger.info('ECG data detected. Processing with improvedRR.')
            if signal_data['ecg'].empty:
                logger.warning('ECG data column is empty.')
            else:
                # Add logging inside improvedRR to ensure it's being executed
                beat_list = beatcaller(signal_data['ecg'], signal_data['ts'], signal_data, analysis_parameters=analysis_parameters)
                if beat_list.empty:
                    logger.warning('No heartbeats detected by improvedRR.')
                else:
                    logger.info(f'Heartbeats detected: {len(beat_list)}')
                if analysis_parameters.get('Pneumo_Mode') != '1':
                    beat_list.to_csv(
                        os.path.join(
                            output_path,f'{MUID}_{PLYUID}_beats.csv'
                            ),
                        index = False
                        )
                else:
                    new_output=os.path.join(output_path,f'{RUID}_beats.csv')
                    if os.path.exists(new_output):
                        beat_list.to_csv(
                            os.path.join(
                                output_path,f'{RUID}_beats_2.csv'
                                ),index=False
                            )
                        logger.warning(f'Duplicate files exist for {RUID}. The beatlist_2 file is generated from {file}.')
                        egregious.append(f'Duplicate files exist for {RUID}. The beatlist_2 file is generated from {file}.')
                    else:
                        beat_list.to_csv(
                            os.path.join(
                                output_path,f'{RUID}_beats.csv'
                                ), index=False
                            )
        else:
            logger.warning('No ecg channel detected in signal data.')
        
        
        
        # recalibrate gas concentration
        signal_data['calibrated_o2'], signal_data['calibrated_co2'], \
        breath_list['calibrated_o2'], breath_list['calibrated_co2'] = \
            calibrate_gas(
                signal_data,
                breath_list,
                auto_criteria,
                manual_selection,
                PLYUID_or_RUID,
                timestamp_dict,
                logger
                )
    
        # assign breaths to blocks and selections
        automated_selections_filters = \
            create_filters_for_automated_selections(
                signal_data,
                breath_list,
                auto_criteria,
                timestamp_dict,
                analysis_parameters,
                logger
                )
    
        manual_selections_filters = \
            create_filters_for_manual_selections(
                breath_list,
                manual_selection,
                PLYUID_or_RUID,
                logger
                )
    
        # apply volumetric calibrations
        calibration_parameters = collect_calibration_parameters(
            breath_list,
            analysis_parameters,
            animal_metadata,
            PLYUID_or_RUID,
            auto_criteria,
            manual_selection,
            automated_selections_filters,
            manual_selections_filters,
            logger
            )
        
        advanced_breath_parameters = \
            apply_volumetric_and_respiratory_calculations(
                breath_list, 
                calibration_parameters, 
                automated_selections_filters, 
                manual_selections_filters,
                analysis_parameters,
                logger)
        
        # review/revise automated filters in case volumetric 
        # calibration changes inclusion                
        revised_automated_selections_filters = \
            create_filters_for_automated_selections(
                signal_data,
                breath_list.merge(
                    advanced_breath_parameters, 
                    on = 'ts_inhale',
                    how = 'left',
                    ),
                auto_criteria,
                timestamp_dict,
                analysis_parameters,
                logger                        
                )
            
        revised_calibration_parameters = collect_calibration_parameters(
                breath_list.merge(
                    advanced_breath_parameters, 
                    on = 'ts_inhale',
                    how = 'left',
                    ),
                analysis_parameters,
                animal_metadata,
                PLYUID_or_RUID,
                auto_criteria,
                manual_selection,
                revised_automated_selections_filters,
                manual_selections_filters,
                logger
                )
        
        revised_advanced_breath_parameters = \
            apply_volumetric_and_respiratory_calculations(
                breath_list.merge(
                    advanced_breath_parameters, 
                    on = 'ts_inhale',
                    how = 'left',
                    ), 
                revised_calibration_parameters, 
                revised_automated_selections_filters, 
                manual_selections_filters, 
                analysis_parameters,
                logger)
    
        # merge values for calibrated volumes into the Breath_List
        breath_list = breath_list.merge(
            revised_advanced_breath_parameters, 
            on = 'ts_inhale',
            how = 'left',
            )
        
        # Include the version of this script as an Analysis Parameter
        analysis_parameters['Version'] = __version__
        
        output_list = create_output(
            breath_list, 
            analysis_parameters, 
            animal_metadata, 
            auto_criteria, 
            manual_selection, 
            PLYUID_or_RUID, 
            revised_automated_selections_filters, 
            manual_selections_filters, 
            column_dictionary, 
            logger
            )
        
        # Prepare the output for STAGG - only generate if not empty
        if len(
                output_list[
                    (output_list['AUTO_IND_INCLUDE'] == 1) |
                    (output_list['MAN_IND_INCLUDE'] == 1)
                    ]
                )<1:
            if analysis_parameters.get('Pneumo_Mode') != '1':
                
                logger.exception(
                    "{}_{} does not have any includable breaths ".format(
                        MUID,PLYUID
                        )+
                    "- no JSON file will be produced."+
                    "This is probably due to problems during "+
                    "sample collection or settings that are not "+
                    "appropriate for the current file"
                    )
                if status_queue: status_queue.put_nowait(
                    f'{pid}: No STAGG output produced, no includable breaths found'+
                    ' - check of your settings recommended'
                    )
            else:
                logger.info(
                    f'{RUID} run in Pneumo Mode. No '+
                    'includable breaths defined by automated or '+
                    'manual selection settings identified. This is '+
                    'expected if no manual or automated settings were '+
                    'provided.'
                    )
                if status_queue: status_queue.put_nowait(
                    f'{pid}: No STAGG output produced, no includable breaths found'+
                    ' - this might be expected if no manual or automated '+
                    'settings were provided.'
                    )
        else:
            if analysis_parameters.get('Pneumo_Mode') != '1':
                logger.info(
                    "{}_{} has includable breaths {}".format(
                        MUID,PLYUID,len(output_list)
                        )+
                    "- JSON file will be produced."
                    )
                if status_queue: status_queue.put_nowait(
                    f'{pid}: {len(output_list)} includable breaths found - SASSI'+
                    f' Output should be available at {output_path}'
                    )
                        
                output_list[
                    (output_list['AUTO_IND_INCLUDE'] == 1) |
                    (output_list['MAN_IND_INCLUDE'] == 1)
                    ].to_json(
                        os.path.join(
                            output_path,
                            '{}_{}.json'.format(MUID,PLYUID)
                            )
                        )
            else:
                logger.info(
                    "{} has includable breaths {}".format(
                        RUID,len(output_list)
                        )+
                    "- JSON file will be produced."
                    )
                output_list[
                    (output_list['AUTO_IND_INCLUDE'] == 1) |
                    (output_list['MAN_IND_INCLUDE'] == 1)
                    ].to_json(
                        os.path.join(
                            output_path,
                            '{}.json'.format(RUID)
                            )
                        )
        
        # provide 'all breath' breathlist if requested
        if analysis_parameters.get('All_Breath_Output') == '1':
            if analysis_parameters.get('Pneumo_Mode') != '1':
                output_list.to_csv(
                    os.path.join(
                        output_path,f'{MUID}_{PLYUID}_all.csv'
                        )
                    )
                    
            else:
                new_breathlist=os.path.join(output_path,f'{RUID}_all_breathlist.csv')
                if os.path.exists(new_breathlist):
                    output_list.to_csv(
                        os.path.join(
                            output_path,f'{RUID}_all_breathlist_2.csv'
                            )
                        )
                    logger.warning(f'Duplicate files exist for {RUID}. The breathlist_2 file is from {file}.')
                    egregious.append(f'Duplicate files exist for {RUID}. The breathlist_2 file is from {file}.')
                else:
                    output_list.to_csv(
                        os.path.join(
                            output_path,f'{RUID}_all_breathlist.csv'
                            )
                        )
                    
        # provide 'aggregated breath data' output if requested
        if 'Aggregate_Output' in analysis_parameters:
            if auto_criteria is not None:
                if analysis_parameters['Aggregate_Output'] != '':
                    group_by_vars = \
                        analysis_parameters[
                            'Aggregate_Output'
                            ].split('@')
                    group_by_vars = [i for i in group_by_vars if i!=\
                                     'Man_Condition']
                    temp_df1 = output_list.groupby(
                        group_by_vars
                        ).first()
                    temp_df2 = output_list.groupby(
                        group_by_vars
                        ).mean()
                    temp_df3 = output_list.groupby(
                        group_by_vars
                        ).count()
                
                    for k in temp_df2.columns:
                        temp_df1.loc[:,k] = temp_df2[k]
                        
                    temp_df1['N'] = temp_df3['Mouse_And_Session_ID']
                    
                    temp_df1['AGG_VF'] = \
                        60/temp_df1['Breath_Cycle_Duration']
                    temp_df1['AGG_VE'] = \
                        temp_df1['AGG_VF'] * \
                            temp_df1['VT__Tidal_Volume_corrected']
                    temp_df1['AGG_VEpg'] = \
                        temp_df1['AGG_VF'] * \
                            temp_df1[
                                'VTpg__Tidal_Volume_per_gram_corrected'
                                ]
                    temp_df1['AGG_VEVO2'] = \
                        temp_df1['AGG_VE'] / temp_df1['VO2']
                    
                    if analysis_parameters.get('Pneumo_Mode') != '1':
                        temp_df1.to_csv(
                            os.path.join(
                                output_path,
                                f'{MUID}_{PLYUID}_agg_auto.csv'
                                )
                            )
                    else:
                        temp_df1.to_csv(
                            os.path.join(
                                output_path,
                                f'{RUID}_agg_auto.csv'
                                )
                            )
            if manual_selection is not None:
                if analysis_parameters['Aggregate_Output'] != '':
                    group_by_vars = \
                        analysis_parameters[
                            'Aggregate_Output'
                            ].split('@')
                    group_by_vars = [i for i in group_by_vars if i!=\
                                     'Auto_Condition']
                    temp_df1 = output_list.groupby(
                        group_by_vars
                        ).first()
                    temp_df2 = output_list.groupby(
                        group_by_vars
                        ).mean()
                    temp_df3 = output_list.groupby(
                        group_by_vars
                        ).count()
                
                    for k in temp_df2.columns:
                        temp_df1.loc[:,k] = temp_df2[k]
                        
                    temp_df1['N'] = temp_df3['Mouse_And_Session_ID']
                    
                    temp_df1['AGG_VF'] = \
                        60/temp_df1['Breath_Cycle_Duration']
                    temp_df1['AGG_VE'] = \
                        temp_df1['AGG_VF'] * \
                            temp_df1['VT__Tidal_Volume_corrected']
                    temp_df1['AGG_VEpg'] = \
                        temp_df1['AGG_VF'] * \
                            temp_df1[
                                'VTpg__Tidal_Volume_per_gram_corrected'
                                ]
                    temp_df1['AGG_VEVO2'] = \
                        temp_df1['AGG_VE'] / temp_df1['VO2']
                    
                    if analysis_parameters.get('Pneumo_Mode') != '1':
                        temp_df1.to_csv(
                            os.path.join(
                                output_path,
                                f'{MUID}_{PLYUID}_agg_man.csv'
                                )
                            )
                    else:
                        temp_df1.to_csv(
                            os.path.join(
                                output_path,
                                f'{RUID}_agg_man.csv'
                                )
                            )
            
        if not egregious:
            pass
        else:
            file_name=os.path.join(output_path,'EGREGIOUS_ERRORS.csv')
            with open(file_name,mode='w',newline='') as file:
                writer = csv.writer(file)
                
                for row in egregious:
                    writer.writerow([row])
            
           
        
        logger.info('completed file {}'.format(file))
           
    except Exception as e:
        logger.exception(
            'Breath Caller encountered an ERROR: {}'.format(e),
            exc_info=True
            )
    

def run_SASSI_from_paths(
        signal_files : list = None,
        output_path : str = None,
        animal_metadata_path : str = None,
        manual_selections_path : str = None,
        auto_criteria_path : str = None,
        channel_match_path : str = None,
        analysis_parameters_path : str = None,
        loglevel : int = logging.DEBUG,
        parallel_runs : int = 1,
        internal_queue = None,
        status_queue = None
        ):
    
    """
    ...
    
    Inputs 
    ------
        Signal Files : ...
        Output_Path : Filename for Output Data
        Animal_Metadata_Path : File containing Animal Metadata
        Manual_Selections_Path : File containing Manual Selections
        Auto_Criteria_Path : File containing Autmated Selection Criteria
        Channel_Match_Path : File containing channel matching information for 
          signal files
        Analysis_Parameters_Path : File containing Analysis Parameters   
    
    Raises
    ------
    Exception
        Exceptions may be raised if necessary input files are not provided or
        if no suitable breaths have been identified.

    Outputs
    -------
        JSON Output for STAGG
        CSV All Breath Output [optional]
        CSV Aggregate Breath Summary [optional]
        log file

    """
    
    
    if not internal_queue:
        internal_queue = multiprocessing.Manager().Queue()
    
    try:
        
        
        if status_queue: status_queue.put_nowait('main: Starting SASSI')
            
        Launch_Time = datetime.now()
        
        
        if not output_path:
            main_logger_output_path = 'error_log.log'
        else:
            main_logger_output_path = os.path.join(
                output_path, 'main_log.log'
                )
        
        main_logger = start_logger(
            main_logger_output_path,
            level = logging.DEBUG,
            logname = f'SASSI {__version__}'
            )
        main_logger.info(f'Version {__version__}')

        
        log_info_from_dict(
            main_logger,
            {
                'input signals': signal_files,
                'output path': output_path,
                'analysis parameters path': analysis_parameters_path,
                'animal metadata path': animal_metadata_path,
                'manual selections path': manual_selections_path,
                'auto criteria path': auto_criteria_path,
                'channel match path': channel_match_path,
                'parallel_runs': parallel_runs
                },
            'command line argument : '
            )
        
        # collect data from input files - raise exception or warning if needed
        if not analysis_parameters_path:
            raise Exception('No Analysis Parameters File Provided')
        else:
            analysis_parameters = pandas.read_csv(
                analysis_parameters_path,
                sep=',',
                encoding='UTF-8',
                index_col='Parameter'
                )['Setting'].to_dict()
        
        analysis_integers=['apnea_window','percent_X_window','sigh_window','chamber_temperature_trim_size']
        analysis_floats=['minimum_TI','minimum_PIF','minimum_PEF','flowrate','percent_X_value','maximum_percent_X','maximum_DVTV',
                         'temperature_calibration_factor','CO2_calibration_factor','O2_calibration_factor',
                         'minimum_sigh_amplitude_x_local_VT','chamber_temperature_default','ecg_threshfactor',
                         'ecg_filter','ecg_invert','ecg_minRR']
        
        for x in analysis_integers:
            try:
                analysis_parameters[x]=int(analysis_parameters[x])
            except:
                pass
            
        for x in analysis_floats:
            try:
                analysis_parameters[x]=float(analysis_parameters[x])
            except:
                pass
        
        if animal_metadata_path is None or \
                animal_metadata_path == '':
            raise Exception('No Animal Metadata File Provided')
        else:
            if analysis_parameters.get('Pneumo_Mode') != '1':
                animal_metadata = get_animal_metadata(animal_metadata_path)
            else:
                animal_metadata = get_animal_metadata_pneumo(animal_metadata_path)
            if animal_metadata == {}:
                raise Exception('Empty or unparsable Animal Metadata') 
    
        if manual_selections_path and manual_selections_path != 'NONE':
            manual_selection = pandas.read_csv(
                manual_selections_path,
                sep=',',
                encoding='UTF-8')
        else: manual_selection = None
    
        if auto_criteria_path and auto_criteria_path != 'NONE':
            auto_criteria = pandas.read_csv(auto_criteria_path, sep=',')
        else: auto_criteria = None
        
        if channel_match_path and channel_match_path != 'NONE':
            channel_match = pandas.read_csv(
                channel_match_path, 
                sep = ',',
                header = None
                )
        else: channel_match = None
                
        # iterate through files and call breaths (multiprocess if applicable)
        
        if parallel_runs == 0 or \
                parallel_runs > multiprocessing.cpu_count(): 
            parallel_runs = multiprocessing.cpu_count()
            status_queue.put_nowait(
                f'Parallel Runs set to match CPU count {parallel_runs}'
                )
            main_logger.info(
                f'Parallel Runs set to match CPU count {parallel_runs}'
                )
        
        pool = multiprocessing.Pool(parallel_runs)
        pool_tracker = {}
        
        for file_number,file in enumerate(signal_files):
            
                pool_tracker[file_number] = {}
                pool_tracker[file_number]['Current_File_Start_Time'] = \
                    datetime.now()
                
                pool_tracker[file_number]['tracker'] = pool.apply_async(
                    run_SASSI_from_data,
                    (),
                    {
                        'file' : file,
                        'signal_data' : None,
                        'output_path' : output_path,
                        'animal_metadata' : animal_metadata,
                        'manual_selection' : manual_selection,
                        'auto_criteria' : auto_criteria,
                        'channel_match' : channel_match,
                        'analysis_parameters' : analysis_parameters,
                        'loglevel' : loglevel,
                        'status_queue' : internal_queue,
                        'logger' : None,
                        'pid' : file_number
                        }
                    )

                main_logger.info(f'added {file} to pool - pid {file_number}')
        # consider launching and tracking seperately this should permit 
        # external python access to the status_queue
        print ('launch processes')

        def term_handler(signum, frame, logger:logging.Logger=None):
            """
            Handles process termination signal (SIGTERM).
            
            This function ensures that the multiprocessing pool is terminated gracefully,
            logging the termination process, and exits the main process safely.
            
            Args:
                signum: Signal number.
                frame: Current stack frame.
                output_path (str): Path to the output directory for logs.
                file (str): Name of the current file being processed.
                pool (multiprocessing.Pool): Reference to the multiprocessing pool.
                status_queue (multiprocessing.Queue): Queue for status messages.
                logger (logging.Logger, optional): Logger instance. If None, a new logger is created.
            """
            # Pools are not guaranteed to close properly with
            #   Python's garbage collector, so we must manually
            #   'close' the pool and 'terminate' the processes,
            #   'joining' until the processes finish terminating
            # See: https://docs.python.org/3/library/multiprocessing.html#module-multiprocessing.pool
            #   > especially the "Warning" message!
            
            if logger is None:
                log_path = os.path.join(output_path, f"{os.path.basename(file)}.log")
                logger = start_logger(log_path, level=logging.DEBUG, logname=os.path.basename(file))

            pool.close()
            pool.terminate()
            logger.info('Terminating pool...')
            status_queue.put_nowait((1, "Terminating pool..."))
            pool.join()
            logger.info('The pool is terminated... No more swimming.')
            status_queue.put_nowait((1, "Terminated pool!"))

            # NOTE: here we could also do something to close more gracefully, such as retaining
            #   any progress so far, or saving any in-progress processing, or outputting some
            #   kind of status file, or whatever you may want to do...

            # Exit main process
            sys.exit()

        # Assign SIGTERM handler to handle parent process termination
        sig.signal(sig.SIGTERM, term_handler)
        
        completed = 0
        while pool_tracker:        
            time.sleep(0.5) 
            
            for pid,track_dict in pool_tracker.copy().items():
                if track_dict['tracker'].ready():
                    if track_dict['tracker'].successful():
                        if status_queue: 
                            status_queue.put_nowait(
                                f'{pid}: {os.path.basename(signal_files[pid])} '+
                                'completed successfully'
                            )
                        
                        pool_tracker.pop(pid)
                        main_logger.info(
                            f'{pid} - {signal_files[pid]} completed successfully'
                            )
                        
                    else:
                        if status_queue: 
                            status_queue.put_nowait(
                                f'{pid}: {os.path.basename(signal_files[pid])} '+
                                ' returned with error'
                            )
                        
                        pool_tracker.pop(pid)
                        main_logger.error(
                            f'{pid} - {signal_files[pid]} returned with error'
                            )
                else: 
                    pass # process is still running
            
        
            # collect queue of status messages
            while not internal_queue.empty():
                current_entry = internal_queue.get_nowait()
                if status_queue: status_queue.put_nowait(f'{pid}: {current_entry}')
                print(f'---{current_entry}')
                
            # report on progress
            current_completed = len(signal_files) - len(pool_tracker)
            if current_completed > completed:
                if status_queue: status_queue.put_nowait(

                        (
                            1,
                            f'still running {len(pool_tracker)}'+
                            f' of {len(signal_files)}'
                            )

                    )
                completed = current_completed
            

        # clear queue in case last addition wasn't removed in prior loop
        while not internal_queue.empty():
            current_entry = internal_queue.get_nowait()
            if status_queue: status_queue.put_nowait(f'{pid}: {current_entry}')
            print(f'---{current_entry}')
        
        
        # close out finished processes
        pool.close()
        pool.join()
        
    
                
        print('all done')
                
    
        
        Finish_Time=datetime.now()
        if status_queue: status_queue.put_nowait(
                    'main: Finished - Duration of Run : '+
                    f'{convert_seconds_to_time_text(Finish_Time-Launch_Time)}'
                )
                
        # This message signals the watcher process to complete
        if status_queue: status_queue.put_nowait("DONE")

        
    except Exception as e:
        
        main_logger.exception(
            'Breath Caller encountered an ERROR: {}'.format(e),
            exc_info = True
            )
    
        


# %% main

def main():
    # collect command line arguments
    
    parser = argparse.ArgumentParser(description='Automated Breath Caller')
    parser.add_argument(
        '-i','--input_directory', help='Path containing signal files for input'
        )
    parser.add_argument(
        '-f','--filename', action='append', help='names of signal files,\
            declare multiple times to build a list of files'
        )
    parser.add_argument(
        '-o', '--output_path', help='Path to directory for output files'
        )
    parser.add_argument(
        '-m', '--manual_selection', help='Path to Manual Selection File'
        )
    parser.add_argument(
        '-a', '--animal_metadata', help='Path to Animal MetaData File'
        )
    parser.add_argument(
        '-c', '--automated_criteria', help='Path to Automated Criteria File'
        )
    parser.add_argument(
        '-n', '--channel_match', help='Path to Channel Match Settings File'
        )
    parser.add_argument(
        '-p', '--analysis_parameters', help='Path to Analysis Parameters File'
        )
    parser.add_argument(
        '-l', '--parallel_runs', 
        help='number of parallel runs to process, default = 1'
        )

    args, others = parser.parse_known_args()
    
    
    
    if args.input_directory is None or args.filename is None:
        Signal_Files = gui_open_filenames(
            {
                'title': 'Load Ascii Files With Trace Data',
                'filetypes': [('txt', '.txt'), ('all files', '.*')]
                }
            )
    else:
        Signal_Files = []
        for file in args.filename:
            Signal_Files.append(os.path.join(args.input_directory, file))

    if args.output_path is None:
        Output_Path = gui_directory(
            {'title': 'Choose Filename for Output Data'}
            )
    else:
        Output_Path = args.output_path

    if args.animal_metadata is None:
        Animal_Metadata_Path = gui_open_filename(
            {'title': 'Select file containing Animal Metadata'}
            )
    else:
        Animal_Metadata_Path = args.animal_metadata

    if args.manual_selection is None:
        Manual_Selections_Path = gui_open_filename(
                {'title': 'Select file containing Manual Selections'}
                )
    else:
        Manual_Selections_Path = args.manual_selection

    if args.automated_criteria is None:
        Auto_Criteria_Path = gui_open_filename(
            {'title': 'Select file containing Autmated Selection Criteria'}
            )
    else:
        Auto_Criteria_Path = args.automated_criteria
        
    if args.channel_match is None:
        Channel_Match_Path = gui_open_filename(
            {'title': 'Select file containing Channel Match Settings'}
            )

    if args.analysis_parameters is None:
        Analysis_Parameters_Path = gui_open_filename(
            {'title': 'Select file containing Analysis Parameters'}
            )
    else:
        Analysis_Parameters_Path = args.analysis_parameters

    if args.parallel_runs is None:
        try:
            Parallel_Runs = gui_get_int(
                'How many parallel runs?',
                'How many parallel runs? (default = 1)',
                1
                )
                
        except:
            Parallel_Runs = 1
    else:
        Parallel_Runs = int(args.parallel_runs)

    try:
        multiprocessing.set_start_method('spawn')
    except Exception:
        print('multiprocessing environment already set...')
    multiprocessing.freeze_support()
    
    
    
    run_SASSI_from_paths(
        signal_files = Signal_Files,
        output_path = Output_Path,
        animal_metadata_path = Animal_Metadata_Path,
        manual_selections_path = Manual_Selections_Path,
        auto_criteria_path = Auto_Criteria_Path,
        channel_match_path = Channel_Match_Path,
        analysis_parameters_path = Analysis_Parameters_Path,
        parallel_runs = Parallel_Runs
        )
    
    
       
# %% Run the Main Program
if __name__ == '__main__':
    main()



