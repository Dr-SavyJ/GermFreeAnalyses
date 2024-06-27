# -*- coding: utf-8 -*-
"""
Selection Helper for Autoresuscitation Measurements (SHAM)
""
or
Automated Autoresuscitation Assay Analysis Helper (AAAH!)
"surprisingly helpful for autoresuscitation measurements and experiments"
or
Autoresuscitation Selection Shortcuts (ASS)
"you can take that autoresuscitation data and stick it in your ASS"


For use with the Ray Lab automated autoresuscitation system and PCC software,
and Breathe Easy with SASSI and STAGG.

Copyright (C) 2022  Christopher Scott Ward
Additional contributions from Savanah Lusk, Dipak Patel, Russel Ray

***
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

=Features=
*intake breathlists and timestampes
*provide summary of basic parameters and challenge scoreboard
    *Baseline summary values
     (defined by baseline settings [any non challenge period])
        *Avg. VF
        *Avg. Breath Duration
        *Avg. VT
        *Avg. VT/g
        *Avg. VE
        *Avg. VE/g
        *Avg. HR
        *Avg. RR
    *Sub-challenge summary values
     (prefill, exposure minus delay, early recovery, late recovery)
     definitions needed for:
         'prefill',
         'exposure minus delay',
         'early recovery',
         'late recovery'
        *Avg. VF
        *Avg. Breath Duration
        *Avg. VT
        *Avg. VT/g
        *Avg. VE
        *Avg. VE/g
        *Avg. HR
        *Avg. RR
    *Timestamp Challenge Start
    *Timestamp Exposure Start
    *Timestamp Apnea (Beginning of Breath)
    *Timestamp Apnea (+Trigger Duration)
    *Timestamp Recovery Gas
    *Timestamp post apnea gasps
    *Timestamp Recovery Gasp 1
    *Timestamp Recovery Gasp 2
    *Timestamp 1st Accum Recovery Start
    *Timestamp Consecutive Recovery Met
    *Timestamp Consecutive Recovery Start
    *Exposure Time (s)
    *Latency from Apnea to Recovery Gas
    *Apnea to 1st Gasp Latency (s)
    *1st Gasp Volume
    *1st Gasp Volume/g
    *1st Gasp to 2nd Gasp Latency (s)
    *2nd Gasp Volume
    *2nd Gasp Volume/g
    *Episode (s)
    *live vs post apnea discrepancy
    *latency from apnea/gasp to HR recovery (tunable)
        (set threshold and which condition to use as baseline)
    *latency from apnea/gasp to VF recovery (tunable)
        (set threshold and which condition to use as baseline)
    *Accumulated Recovery Latency (s)
    *Consecutive Recovery Latency (s)
    *# challenges survived

    tunable settings:
    minimum apnea duration:
    minimum apnea to gasp duration:
    prefill delay for hypervent phase: (0 starts at timestamp of gas start)
    early recovery starting breath: (0 starts at 1st gasp)
    condition to use for baseline:
    VF recovery threshold:
    HR recovery threshold:
    Consecutive Recovery Threshold:
    Accumulated Recovery Threshold:


    timestamp harmonizing settings:
    regex for...
    *condition
    *challenge start
    *exposure start
    *apnea detection
    *recovery gas start

    definitions of:
    'prefill': challenge start before gas exposure
    'exposure': challenge gas
        [optionally trimmed from start of gas or before apnea]
    'early recovery': Nth breath after gasp until start of Accum/Consec Recov
    'late recovery': last 30sec before next challenge
    ''

"""


__version__ = "0.1.9"


# %% import libraries
import pandas
import numpy
import re
import ast
import tkinter.filedialog
import tkinter
import logging
import argparse
import sys
import os

#%% import local modules

# Get the directory of the current script
script_directory = os.path.dirname(os.path.realpath(__file__))

# Change the current working directory to the script directory
# os.chdir(script_directory)

# Add the script directory to sys.path if it's not already included
if script_directory not in sys.path:
    sys.path.append(script_directory)
    sys.path.append(os.path.join(script_directory,'..','SASSI'))

# Change directory to SASSI directory
os.chdir('../SASSI')

# Now, attempt to import your modules
try:
    from signal_converters import adi_extract
except ImportError:
    # If the direct import fails, attempt an import directly from wkdir as a fallback
    # Note: This will work if the script_directory houses your modules
    try:
        from .signal_converters import adi_extract
    except ImportError as e:
        print(f"Failed to import modules: {e}")

#%%Change directory for running SHAM.
# os.chdir('../autoresuscitation')

# %% define functions and classes
class SETTINGS:
    def __init__(self):

        self.challenge_pf_threshold_multiple = 2
        self.minimum_apnea_duration = 5
        self.minimum_PIF = 0.05
        self.minimum_PEF = 0.05
        self.minimum_apnea_to_gasp_duration = 15
        self.trim_for_hypervent_phase = "5,5"
        self.early_recovery_starting_breath = 10
        self.late_recovery_start_from_end = 30

        self.condition_to_use_for_baseline = ["Pre-Inject","Baseline"]

        self.Baseline_minimum_bout = 5
        self.Baseline_VF = 250
        self.Baseline_TT = 1
        self.Baseline_isTT = 1
        self.Baseline_HR = 700
        self.Baseline_RR = 999
        self.Baseline_isRR = 1
        self.Baseline_BSD = 0.5
        self.Baseline_DVTV = 0.75
        self.VF_recovery_threshold = 50
        self.HR_recovery_threshold = 63
        self.Consecutive_Recovery_Threshold = 30
        self.Accumulated_Recovery_Threshold = 30
        self.Accumulated_Recovery_Minimum_Bout = 5
        self.IE_recovery_threshold = 60
        self.IE_recovery_bout = 5
        self.cal_vol_mL = 0.02

        # timestamp harmonizing settings:
        # regex for...
        # *condition
        # *challenge start
        # *exposure start
        # *apnea detection
        # *recovery gas start

    def load_from_file(self, filepath, logger=None):

        analysis_parameters = pandas.read_csv(
            filepath, sep=",", encoding="UTF-8", index_col="Parameter"
        )["Setting"].to_dict()
        for k,v in analysis_parameters.items():
            if k.startswith("PM_"):
                attr_name = k[3:]
                try:
                    attr_val = float(v)
                except ValueError:
                    attr_val = v
                    
                if '[' in v:
                    attr_val = ast.literal_eval(v)

                setattr(self, attr_name, attr_val)

                if logger:
                    logger.info(f'"{attr_name}" set to "{attr_val}" from file')

    def save_to_file(self, filepath, logger=None):

        analysis_parameters = {}
        for k in self.__dict__:
            if k.startswith('__'):
                continue
            analysis_parameters[f"PM_{k}"] = self.__dict__[k]
        AP_df = pandas.DataFrame({'Parameter':analysis_parameters.keys(),'Setting':analysis_parameters.values()})
        AP_df.to_csv(filepath,index=False)
        if logger: logger.info(f"analysis parameters saved to file: {filepath}")



def gui_open_filename(kwargs={}):
    """
    This function creates a temporary Tkinter instance that provides a GUI
    dialog for selecting a filename.

    Parameters
    ----------
    kwargs : Dict, optional
        The default is {}.
        *Function calls on tkFileDialog and uses those arguments
      ......
      (declare as a dictionary)
      {"defaultextension":'',"filetypes":'',"initialdir":'',...
      "initialfile":'',"multiple":'',"message":'',"parent":'',"title":''}
      ......

    Returns
    -------
    output_text : String
        String describing the path to the file selected by the GUI.

    """

    root = tkinter.Tk()
    output_text = tkinter.filedialog.askopenfilename(**kwargs)
    root.destroy()
    return output_text


def gui_open_filenames(kwargs={}):
    """
    This function creates a temporary Tkinter instance that provides a GUI
    dialog for selecting [multiple] filenames.

    Parameters
    ----------
    kwargs : Dict, optional
        The default is {}.
        *Function calls on tkFileDialog and uses those arguments
      ......
      (declare as a dictionary)
      {"defaultextension":'',"filetypes":'',"initialdir":'',...
      "initialfile":'',"multiple":'',"message":'',"parent":'',"title":''}
      ......

    Returns
    -------
    output_text : List of Strings
        List of Strings describing the paths to the files selected by the GUI.

    """

    root = tkinter.Tk()
    output_text_raw = tkinter.filedialog.askopenfilenames(**kwargs)
    output_text = root.tk.splitlist(output_text_raw)
    root.destroy()
    return output_text


def gui_directory(kwargs={}):
    """
    This function creates a temporary Tkinter instance that provides a GUI
    dialog for selecting a directory.

    Parameters
    ----------
    kwargs : Dict, optional
        The default is {}.
        *Function calls on tkFileDialog and uses those arguments
          ......
          (declare as a dictionary)
          {"defaultextension":'',"filetypes":'',"initialdir":'',...
          "initialfile":'',"multiple":'',"message":'',"parent":'',"title":''}
          ......

    Returns
    -------
    output_text : String
        Returns the directory path selected by the GUI.

    """

    root = tkinter.Tk()
    output_text = tkinter.filedialog.askdirectory(**kwargs)
    root.destroy()
    return output_text


def initialize_logger(filepath):
    logger = logging.getLogger("SHAM a.k.a. AAAH!")
    logger.setLevel(logging.DEBUG)

    # create file and console handlers to receive logging
    console_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(filepath)
    console_handler.setLevel(logging.DEBUG)
    file_handler.setLevel(logging.DEBUG)

    # create format for log and apply to handlers
    log_format = logging.Formatter(
        "%(asctime)s | %(name)s | %(levelname)s | %(message)s"
    )
    console_handler.setFormatter(log_format)
    file_handler.setFormatter(log_format)

    # add handlers to the logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    # log initial inputs
    logger.info("Initializing Log")

    return logger


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


def load_signal_data(filename, local_logger):
    """
    Creates a dataframe containing plethysmography signal data

    Parameters
    ----------
    filename : string
        path to file containing signal data
    local_logger : instance of logging.logger


    Returns
    -------
    signal_data_assembled : pandas.DataFrame
        dataframe containing contents of signal file (merged into a single
        dataframe if multiple blocks are present)

    """

    if filename.endswith(".adicht"):
        signal_data_assembled = adi_extract.SASSI_extract(
            filename, logger=local_logger
        )
    else:
        print(f"Cannot process this filetype {filename.endswith}.")

    signal_data_assembled.rename(
        columns={"time": "ts", "vent flow": "flow"}, inplace=True
    )

    return signal_data_assembled


def fix_broken_timestamps(
    timestamp_dict,
    recognized_commands_re=None,
    capture_version=None,
    logger=None,
):

    if capture_version is not None:
        pass  # customize recognized_commands to reflect capture_version

    if recognized_commands_re is None:
        recognized_commands_re = [
            "Starting: Calibrating for 0s",
            "Finished: Calibrating",
            "Habituation-1",
            "Pre-Inject",
            "Baseline",
            "Challenge",
            "ChallengeFinished: On Position ,0, Prefilled for 0s,Finished: On Anoxic Air",
            "Starting: On Anoxic Air,0,0,Ongoing: On Position ,0, Prefilled for 0s",
            "Starting: On Anoxic Air,Prefill for 0s,Ongoing: On Position 0, Prefilled for 0s",
            "Finished: On Position ,0, Prefilled for 0s,Finished: On Anoxic Air",
            "Finished: On Position 0, Prefilled for 0s,Finished: On Anoxic Air",
            "Starting: On Room Air,0,0,Ongoing: On Position ,0, Gas Off",
            "Starting: On Room Air,Prefill for 0s,Ongoing: On Position 0, Gas Off",
            "Finished: On Room Air",
            "Finished: On Room Air,ABORTED",
            "0",
            "R0",
            "Cal 20 Room Air",
            "Room Air",
        ]

    broken_comment_index_1 = 0
    for k, v in timestamp_dict.items():

        re_v = re.sub("[0-9]+", "[0-9]+", v)

        if any(
            [
                re.compile("^" + re_v + "$").search(i)
                for i in recognized_commands_re
            ]
        ):

            if logger:
                logger.info(f"Found: {k}:{v}")
        else:
            if logger:
                logger.info(f"unrecognized arduino com:\n{k}:{v}")

            for i in recognized_commands_re:
                # starts with
                if re.compile("^" + re_v).search(i):
                    broken_comment_index_1 = k

                    if logger:
                        logger.info(
                            f"--likely beginning fragment of command {k}"
                        )
                # ends with
                elif re.compile(re_v + "$").search(i) and re.sub(
                    "[0-9]+",
                    "[0-9]+",
                    timestamp_dict[broken_comment_index_1] + v,
                ) == re.sub("[0-9]+", "[0-9]+", i):
                    if logger:
                        logger.info(f"--likely ending fragment of command {k}")
                    timestamp_dict[broken_comment_index_1] = (
                        timestamp_dict[broken_comment_index_1] + v
                    )
                    timestamp_dict[k] = ""
                    if logger:
                        logger.info(
                            f"FIXED:{timestamp_dict[broken_comment_index_1]}"
                        )
                # middle piece
                elif re.compile(
                    "^"
                    + re.sub(
                        "[0-9]+",
                        "[0-9]+",
                        timestamp_dict.get(broken_comment_index_1, ""),
                    )
                    + re_v
                ).search(i):
                    if logger:
                        logger.info(
                            f"--likely mid fragment of command {broken_comment_index_1}"
                        )
                    timestamp_dict[broken_comment_index_1] = (
                        timestamp_dict[broken_comment_index_1] + v
                    )
                    if logger:
                        logger.info(
                            f"attempting repair:{timestamp_dict[broken_comment_index_1]}"

                        )
                    timestamp_dict[k] = ""
                # full but with gap
                elif re.compile(
                    "^"
                    + re.sub(
                        "[0-9]+",
                        "[0-9]+",
                        timestamp_dict.get(broken_comment_index_1, ""),
                    )
                ).search(i) and re.compile(re_v + "$").search(i):
                    if logger:
                        logger.info(
                            f"likely ending of fragment with middle gap {broken_comment_index_1}"
                        )
                    gap_finder = re.compile(
                        "^(?P<frag1>"
                        + timestamp_dict.get(broken_comment_index_1, "")
                        + ")(?P<gap>.*)(?P<frag2>"
                        + re_v
                        + ")$"
                    )
                    gap_search = gap_finder.search(i)
                    gap_contents = gap_search.group("gap")
                    timestamp_dict[broken_comment_index_1] = (
                        timestamp_dict[broken_comment_index_1]
                        + gap_contents
                        + v
                    )
                    timestamp_dict[k] = ""
                    if logger:
                        logger.info(
                            f"Fixed: {timestamp_dict[broken_comment_index_1]}"
                        )

    for k in list(timestamp_dict.keys()):
        if timestamp_dict[k] == "":
            timestamp_dict.pop(k)
    return timestamp_dict


def harmonize_timestamps(
    timestamp_dict, capture_version=None, translation_dict=None, logger=None
):

    if capture_version is not None:
        pass  # customize the translation dict to reflect the capture version

    if translation_dict is None:
        # !!!consider moving to external settings, may need to address handling
        #     of timestamps to be omitted from the harmonized output
        translation_dict = {


            "Starting: Cal": "calibration",

            "Finished: Cal": "Signal Preview-1",
            "Habituation": "Habituation",
            "Habituation-1": "Habituation-1",
            "Habituation-2": "Habituation-2",
            "Pre-Inject": "Pre-Inject",
            "Baseline": "Baseline",
            "Challenge": "Challenge",
            "Finished: On Anoxic": "Anoxic",
            "Finished: On Position": "Anoxic",
            "Starting: On Anoxic": "Prefill",
            "Starting: On Room Air": "Apnea(Arduino)",
            "Finished: On Room Air": "Recovery",

            "Cal 20 Room Air": "calibration",
            "Cal 20": "calibration",

            "apnea starts": "Apnea(Arduino)",
            "Room Air": "Baseline",
            "Room Ait": "Baseline",
            "Room AIr": "Baseline",
            "[1-9]": "Anoxic",
            "[1-9][1-9]": "Anoxic",
            "R[1-9]": "Recovery",
            "R[1-9][1-9]": "Recovery"
        }

    translated_results = {}

    for k_orig, v_orig in timestamp_dict.items():
        for k_trans, v_trans in translation_dict.items():
            if re.compile("^" + k_trans).search(v_orig):
                print(1)
                # if k_trans in v_orig:

                translated_results[k_orig] = v_trans
                logger.info(f'{k_orig}: matched to "{v_trans}"')
                break
            elif re.compile(".*,+" + k_trans).search(v_orig):
                print(1)
                # if k_trans in v_orig:

                translated_results[k_orig] = v_trans
                logger.info(f'{k_orig}: matched to "{v_trans}" as part of complex comm i/o')
                break
        if k_orig not in translated_results.keys():
            logger.warning(f"unable to match : {k_orig}:{v_orig}")

    challenge_found = any(value == 'Challenge' for value in translated_results.values())
    prefill_found = any(value == 'Prefill' for value in translated_results.values())
    anoxic_timestamps = [k for k, v in translated_results.items() if v == 'Anoxic']
    
    if not challenge_found:  
        translated_results=handle_missing_challenge(translated_results, anoxic_timestamps, logger)
    if not prefill_found:
        translated_results=handle_missing_prefill(translated_results, anoxic_timestamps, logger)
    translated_results=dict(sorted(translated_results.items()))        
    return translated_results

def handle_missing_challenge(translated_results, anoxic_timestamps, logger):
    try:
        translated_results[min(anoxic_timestamps)-10]='Challenge'
    except:
        logger.warning("Unable to find Challenge comment.")
    return translated_results

def handle_missing_prefill(translated_results, anoxic_timestamps, logger):
    for k in anoxic_timestamps:
        translated_results[k-5]='Prefill'
    return translated_results
    
def get_pneumo_TV(tv, calv, act_calv):
    """
    Calculates pneumotacograph based tidal volume.

    Parameters
    ----------
    tv : Float
        uncorrected tidal volume (V)
    calv : Float
        uncorrected tidal volume from calibration period (V)
    act_calv : Float
        nominal calibration volume (mL)

    Returns
    -------
    Float
        corrected tidal volume (mL)

    """

    return tv / calv * act_calv


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


def multi_filter(input_df, index_col, filter_dict, logger=None):
    """
    generates a dataframe output with a timestamp column "ts" and a 
    boolean filter column "filt" that indicates if the "ts" is passing or not
    for all of the filter criteria

    Parameters
    ----------
    input_df : pandas.DataFrame
        time series data that will be filtered against
    index_col : str
        name of the column containing the timestamp information
    filter_dict : dict 
        containing keys describing the filter and values as a 3 item tuple 
            [column name to filter],
            [comparison e.g. < > = as g l e]
            [value to compare against]
    logger : logging object, optional
        logging object. The default is None.

    Returns
    -------
    filter_output : pandas.DataFrame
        two column dataframe with a 'ts' column indicating timestamps and
        a 'filt' column indicating the pass/fall of that timestamps data as a 
        boolean for all of the filters.

    """

    filter_df = pandas.DataFrame(input_df[index_col])
    for k, v in filter_dict.items():
        # k describes a comparison name
        # v should be a 3 part tuple indicatling in position [0]
        # the column name in input_df
        # in position [1] g,l,e,ge,le corresponding to >,<,==,>=,<=
        # and indicating in position [2] the value to test against
        if v[1] == "g":
            filter_df[k] = input_df[v[0]] > v[2]
        elif v[1] == "ge":
            filter_df[k] = input_df[v[0]] >= v[2]
        elif v[1] == "l":
            filter_df[k] = input_df[v[0]] < v[2]
        elif v[1] == "le":
            filter_df[k] = input_df[v[0]] <= v[2]
        elif v[1] == "e":
            filter_df[k] = input_df[v[0]] == v[2]
        elif logger:
            filter_df[k] = 1
            logger.warning(f"unable to apply filter {v} - filter skipped")

    filter_output = pandas.DataFrame(
        {
            "ts": input_df[index_col],
            "filt": filter_df[filter_dict.keys()].min(axis=1),
        }
    )

    return filter_output


def resample_and_merge_filters(
    filter1_orig, filter2_orig, sample_int, minimum_bout, logger=None
):
    """
    

    Parameters
    ----------
    filter1_orig : pandas.DataFrame with 'ts' and 'filt' columns
        filter similar to that produced by multi_filter()
    filter2_orig : pandas.DataFrame with 'ts' and 'filt' columns
        filter similar to that produced by multi_filter()
    sample_int : float
        interval between samples to use for populating timing between 
        filter1_orig and filter2_orig timestamps
    minimum_bout : float
        minimum duration to use as requirement for filter1 and filter2 to both 
        pass so the output filter is all passing 
    logger : logging object, optional
        a logging object. The default is None.

    Returns
    -------
    filter1
        pandas.DataFrame of Timestamps from filter1_orig that passed in both 
        filter1_orig and filter2_orig for the minimum_bout duration
    filter2
        pandas.DataFrame of Timestamps from filter2_orig that passed in both 
        filter1_orig and filter2_orig for the minimum_bout duration
    resample_df
        pandas.DataFrame of Timestamps and filter columns indicating a 
        Timestamps status for filter1, filter2, and a merger of both

    """


    filter1 = filter1_orig.copy().rename(
        {"filt": "f1", "ts": "ts1"}, axis="columns"
    )
    filter2 = filter2_orig.copy().rename(
        {"filt": "f2", "ts": "ts2"}, axis="columns"
    )
    precision = (
        int(numpy.format_float_scientific(sample_int).split("e")[1]) * -1
    )
    filter1.loc[:, "round_ts"] = round(filter1["ts1"], precision)
    filter2.loc[:, "round_ts"] = round(filter2["ts2"], precision)

    resample_df = filter1.merge(filter2, how="outer", on="round_ts")

    resample_df["round_ts_as_sec"] = pandas.to_datetime(
        resample_df["round_ts"], unit="s"
    )

    resample_df = resample_df.set_index("round_ts_as_sec")

    resample_df = resample_df.sort_index()


    resample_df["f1"] = resample_df["f1"].ffill()
    resample_df["f2"] = resample_df["f2"].ffill()


    resample_df.loc[:, "f1_and_f2"] = resample_df["f1"] & resample_df["f2"]

    # !!! unit testing would be ideal for this
    # rolling minimum window requires all samples in window to be true to
    # be output as true, 'or' to combine a left and right justified filter
    # provides an output that passes with the minimum bout duration
    rolling_right = (
        resample_df["f1_and_f2"]
        .rolling(f"{int(minimum_bout/2*1000)}ms")
        .min()
        .fillna(0)
        .astype(bool)
    )
    # series is inverted for rolling window and inverted again for output
    # (rolling function only works right justified, workaround provides a
    # left justified rolling window)
    rolling_left = (
        resample_df.loc[::-1, "f1_and_f2"]
        .rolling(f"{int(minimum_bout/2*1000)}ms")
        .min()
        .fillna(0)
        .astype(bool)
        .sort_index()
    )
    resample_df.loc[:, "f1_and_f2"] = rolling_right | rolling_left

    filter1 = (
        filter1_orig.merge(
            resample_df[["ts1", "f1_and_f2"]],
            how="left",
            left_on="ts",
            right_on="ts1",
        )
        .fillna(0)
        .astype(bool)
    )
    filter2 = (
        filter2_orig.merge(
            resample_df[["ts2", "f1_and_f2"]],
            how="left",
            left_on="ts",
            right_on="ts2",
        )
        .fillna(0)
        .astype(bool)
    )

    if logger:
        logger.info(
            f"resample_and_merge_filters: f1 pass={filter1.filt.sum()}"
        )
        logger.info(
            f"resample_and_merge_filters: f2 pass={filter2.filt.sum()}"
        )
        logger.info(
            f"resample_and_merge_filters: f1_and_f2 f1 pass={filter1.f1_and_f2.sum()}"
        )
        logger.info(

            f"resample_and_merge_filters: f1_and_f2 f2 pass={filter2.f1_and_f2.sum()}"

        )

        if filter1["f1_and_f2"].sum() == 0 or filter2["f1_and_f2"].sum() == 0:
            logger.warning("resample_and_merge_filters yeild 0 passing")
    return (
        filter1["f1_and_f2"],
        filter2["f1_and_f2"],
        resample_df[["round_ts", "f1_and_f2"]],
    )


def get_ts_for_accumulated_value(
    input_df_orig, filter_series_orig, accum_col, target, ts_col
):
    input_df = input_df_orig.copy()
    filter_series = filter_series_orig.copy()
    input_df.loc[filter_series, "passing"] = input_df[filter_series][accum_col]
    input_df.loc[:, "accum"] = input_df["passing"].cumsum()
    ts = input_df[input_df["accum"] >= target][ts_col].min()
    return ts

def remove_data_after_corruption(ordered_dict, target_values=['data corruption boundary >','< data corruption boundary']):
    new_dict={}
    for key, value in ordered_dict.items():
        if value in target_values:
            break
        new_dict[key]=value
    return new_dict
            


# %%


def run_sham(signal_paths, breathbeat_dir, settings_path, output_path):

    #%% set up logging
    logger = initialize_logger(os.path.join(output_path, "autores_log.log"))

    # %%
    # Prepare Aggregate DataFrame
    aggregate_baseline = {}
    aggregate_challenge = {}
    aggregate_timestamps = pandas.DataFrame()
    aggregate_summary = {}

    # gather settings
    settings = SETTINGS()
    settings.load_from_file(settings_path, logger=logger)


    # parse and match files

    file_groups = {}
    #%%
    for f in signal_paths:
        #%%
        try:
            # Extract RUID
            ruid = extract_ruid(f)
            
            # Create filegroups for each of the expected inputs.
            file_groups[ruid] = {
                "signal": f,
                "breathlist": os.path.join(
                    breathbeat_dir, f"{ruid}_all_breathlist.csv"
                ),
                "beatlist": os.path.join(breathbeat_dir, f"{ruid}_beats.csv"),
            }

            logger.info(f"Loading data (breaths,beats,signals) for {ruid}")
            
            # Loading signal file
            signals = load_signal_data(f, logger)
            
            # Loading breathlist
            breath_list = pandas.read_csv(file_groups[ruid]["breathlist"])
            
            # Add BSD to breath_list
            signals.loc[:, "mov_avg_flow"] = (
                signals["flow"]
                .rolling(
                    int(1 / (signals["ts"].iloc[1] - signals["ts"].iloc[0])),

                    center=True,
                )
                .mean()
            )
            
            # determine sample interval
            sample_interval = signals["ts"].iloc[1] - signals["ts"].iloc[0]
            breath_list.loc[:, "BSD"] = abs(
                breath_list.merge(
                    signals,

                    how="left",
                    left_on="Timestamp_Inspiration",
                    right_on="ts",
                )["mov_avg_flow"]
            )

            # gather beatlist
            try:
                beat_list = pandas.read_csv(file_groups[ruid]["beatlist"])
            except:
                beat_list = pandas.DataFrame({"ts": [], "RR": [],"HR":[]})
                if logger:
                    logger.warning("no ecg data")

            # add isRR to beat_list
            beat_list.loc[:, "isRR"] = calculate_irreg_score(beat_list["RR"])
            
            # collect timestamps
            logger.info("collecting timestamps")
            timestamp_dict = dict(
                zip(
                    signals.ts[signals.comment.dropna().index],
                    [i for i in signals.comment.dropna()],
                )
            )
            for key, value in timestamp_dict.items():
                if value.startswith("#"):
                    timestamp_dict[key] = value[3:]
            
            # Remove corrupted data          
            timestamp_dict = remove_data_after_corruption(timestamp_dict)
            
            # Fix broken timestamps
            timestamp_dict = fix_broken_timestamps(
                timestamp_dict,
                recognized_commands_re = None,
                logger = logger
                )
            
            # Harmonize timestamps per dictionary in function
            harmonized_timestamps = harmonize_timestamps(
                timestamp_dict, logger=logger
            )
            
            # Identify challenges
            try:
                # identify challenge rounds
                logger.info("extracting challenge round info")
                challenge_timestamp = None
                challenge_list = []
                for k, v in harmonized_timestamps.items():
                    if v == "Challenge":
                        challenge_timestamp = k
                        break
                if challenge_timestamp is None:
                    # search for first 'prefill timestamp'
                    for k, v in harmonized_timestamps.items():
                        if v == "Prefill":
                            challenge_timestamp = k
                            break
                    if challenge_timestamp is None:

                        raise Exception("Challenge Mode Timestamp not found")

                open_challenge = False
                trial_counter = 0
                current_challenge = {}

                for k, v in harmonized_timestamps.items():

                    if k <= challenge_timestamp:
                        continue
                    if v == "Prefill":
                        if open_challenge is True:
                            logger.warning(

                                f"challenge trial timestamp anomoly -Prefill-{trial_counter}"
                            )
                        else:
                            trial_counter += 1
                            current_challenge = {
                                "trial_number": trial_counter,
                                "Prefill": k,
                            }
                            open_challenge = True
                    elif v == "Recovery":
                        current_challenge[v] = k
                        open_challenge = False

                        # !!! would be better to move this list to a settings
                        #     (probably shared with the setting for the
                        #     harmonization dict)
                        if (
                            len(
                                set(current_challenge).intersection(
                                    {
                                        "Prefill",
                                        "Anoxic",
                                        "Apnea(Arduino)",
                                        "Recovery",
                                    }
                                )
                            )
                            != 4
                        ):
                            current_challenge["trial_number"] = trial_counter

                            logger.warning(
                                f"Incomplete challenge timestamps -{trial_counter} - {current_challenge}"
                            )
                        challenge_list.append(current_challenge)


                        current_challenge = {}

                    elif v != "Prefill" and not open_challenge:

                        logger.warning(

                            f"challenge anomoly - {v}-{trial_counter}"
                        )
                        continue
                    elif v in current_challenge:

                        logger.warning(
                            f"challenge trial timestamp anomoly - {v}-{trial_counter}"
                        )
                        current_challenge["trial_number"] = trial_counter
                        challenge_list.append(current_challenge)


                        trial_counter += 1
                        current_challenge = {
                            "trial_number": trial_counter,
                            v: k,
                        }

                    else:
                        current_challenge[v] = k
                # add last challenge if it was left open
                if current_challenge != {}:

                    challenge_list.append(current_challenge)



            except Exception:
                logger.warning(
                    f"unable to determine challenge info{os.path.basename(f)}",
                    exc_info=True,
                )


            #%
            # generate calibration and pre challenge values
            try:
                calibration_start = max(
                    [

                        k
                        for k, v in harmonized_timestamps.items()
                        if v == "calibration"
                    ]
                )
                calibration_end = min(
                    [
                        k
                        for k in harmonized_timestamps.keys()
                        if k > calibration_start
                    ]
                )

                calibration_VT = breath_list[
                    (breath_list["Timestamp_Inspiration"] > calibration_start)
                    & (breath_list["Timestamp_Inspiration"] < calibration_end)
                ]["Tidal_Volume_uncorrected"].mean()
                calibration_VF = breath_list[
                    (breath_list["Timestamp_Inspiration"] > calibration_start)
                    & (breath_list["Timestamp_Inspiration"] < calibration_end)
                ]["VF"].mean()
                calibration_Breath_Cycle_Duration = breath_list[
                    (breath_list["Timestamp_Inspiration"] > calibration_start)
                    & (breath_list["Timestamp_Inspiration"] < calibration_end)
                ]["Breath_Cycle_Duration"].mean()

                # correct calibrated volume columns in breath_list
                breath_list.loc[
                    :, "VT__Tidal_Volume_corrected"
                ] = get_pneumo_TV(
                    breath_list["Tidal_Volume_uncorrected"],
                    calibration_VT,
                    settings.cal_vol_mL
                )

                breath_list.loc[:,'VTpg__Tidal_Volume_per_gram_corrected'] = \
                    breath_list['VT__Tidal_Volume_corrected'] \
                        / breath_list['Weight']
    
                breath_list.loc[:,'VE__Ventilation'] = \
                    breath_list['VT__Tidal_Volume_corrected'] * \
                        breath_list['VF']
                breath_list['VEpg__Ventilation_per_gram'] = \
                    breath_list['VE__Ventilation'] / \
                        breath_list['Weight']
                        
                        
                #%
                # generate pre-challenge summary values
                if type(settings.condition_to_use_for_baseline) is not list:
                    settings.condition_to_use_for_baseline = [
                        settings.condition_to_use_for_baseline    
                    ]
                    
                pre_challenge_dict = {}
                recovery_baselines = {}
                for cond in settings.condition_to_use_for_baseline:

                    recovery_baselines[cond] = {}
                    
                    cond_start = max(
                        [
                            k for k,v in harmonized_timestamps.items()
                            if v == cond
                        ]    
                    )
                    cond_end = min(
                        [
                            k for k in list(harmonized_timestamps.keys())
                            + [signals.ts.iloc[-1]]
                            if k > cond_start
                        ]    
                    )
                    
                    # !!! note post hoc analysis is using 'is' instead of 'cv' for
                    # variation acceptablility for inclusion
                    time_filter_breaths = multi_filter(
                        breath_list,
                        "Timestamp_Inspiration",
                        {
                            "start": (
                                "Timestamp_Inspiration",
                                "ge",
                                cond_start
                            ),
                            "end": (
                                "Timestamp_Inspiration", 
                                "l", 
                                cond_end
                            )
                        },
                        logger=logger
                    )
                    
                    time_filter_beats = multi_filter(
                        beat_list,
                        "ts",
                        {
                            "start": (
                                "ts",
                                "ge",
                                cond_start
                            ),
                            "end": (
                                "ts", 
                                "l", 
                                cond_end
                            )
                        },
                        logger=logger
                    )
                    
                    breath_filter = multi_filter(
                        breath_list,

                        "Timestamp_Inspiration",
                        {
                            "start": (
                                "Timestamp_Inspiration",
                                "ge",

                                cond_start
                            ),
                            "end": (
                                "Timestamp_Inspiration", 
                                "l", 
                                cond_end
                            ),
                            "VF": ("VF", "le", settings.Baseline_VF),
                            "TT": (
                                "Breath_Cycle_Duration",
                                "le",
                                settings.Baseline_TT,
                            ),
                            "isTT": ("IS_TT", "le", settings.Baseline_isTT),
                            "BSD": ("BSD", "le", settings.Baseline_BSD),
                            "DVTV": ("DVTV", "le", settings.Baseline_DVTV),
                        },
                        logger=logger,
                    )
                    
                    
                    if beat_list.shape[0] > 1 and settings.Baseline_RR > 0:
                        beat_filter = multi_filter(
                            beat_list,
                            "ts",
                            {
                                "start": (
                                    "ts", 
                                    "ge", 
                                    cond_start
                                ),
                                "end": (
                                    "ts", 
                                    "l", 
                                    cond_end
                                ),
                                "HR": ("HR", "le", settings.Baseline_HR),
                                "RR": ("RR", "le", settings.Baseline_RR),
                                "isRR": ("isRR", "le", settings.Baseline_isRR),
                            },
                            logger=logger,
                        )
                        # %
                        logger.info(
                            f"Deriving Quality Baseline Data -{cond}- using Breathing and Heartbeat"
                            
                        )
                        quality_breath, quality_beat = resample_and_merge_filters(
    
                            breath_filter,
                            beat_filter,
                            sample_interval,
                            settings.Baseline_minimum_bout,
    
                            logger=logger,
                        )[0:2]
                    else:
                        logger.info(
                            f"Deriving quality Baseline Data -{cond}- (Breathing Only)"
                        )
                        quality_breath, quality_beat = resample_and_merge_filters(
    
                            breath_filter,
                            breath_filter,
                            sample_interval,
                            settings.Baseline_minimum_bout,
                            logger=logger,
                        )[0:2]

                    logger.info("Collecting Baseline Measures for Recovery Calculations - {cond}")
                    # populate recovery related baseline measures
                    recovery_baselines[cond]['TT'] = \
                        breath_list[quality_breath]['Breath_Cycle_Duration'].mean()
                    recovery_baselines[cond]['VF'] = 60 / recovery_baselines[cond]['TT']
                    recovery_baselines[cond]['VF_alternate'] = \
                        (60 / breath_list[quality_breath]['Breath_Cycle_Duration']).mean()
                    recovery_baselines[cond]['RR'] = \
                        beat_list[quality_beat]['RR'].mean()
                    recovery_baselines[cond]['HR'] = 60 / recovery_baselines[cond]['RR']
                    recovery_baselines[cond]['HR_alternate'] =\
                        (60 / beat_list[quality_beat]['RR']).mean()
                    if 'IE_ratio' not in breath_list.columns:
                        logger.info('old breathlist format - calculating IE ratio')
                        breath_list['IE_ratio'] = \
                            breath_list['Inspiratory_Duration'] / \
                                breath_list['Expiratory_Duration']
                    else:
                        pass
                    recovery_baselines[cond]['IE_duration_ratio'] = \
                        breath_list[quality_breath]['IE_ratio'].mean()

                    logger.info("Populating PreChallenge Summary - {cond}")
    
                    # !!! add functionality to generate additional summaries, with or
                    # without tunable filters (is this needed if SASSI can do it?)

                    pre_challenge_dict[cond] = {
                        "VF": breath_list[quality_breath]["VF"].mean(),
                        "Breath Duration": breath_list[quality_breath][
                            "Breath_Cycle_Duration"
                        ].mean(),
                        "VT": breath_list[quality_breath][
                            "VT__Tidal_Volume_corrected"
                        ].mean(),
                        "VT/g": breath_list[quality_breath][
                            "VTpg__Tidal_Volume_per_gram_corrected"
                        ].mean(),
                        "VE": breath_list[quality_breath][
                            "VE__Ventilation"
                        ].mean(),
                        "VE/g": breath_list[quality_breath][
                            "VEpg__Ventilation_per_gram"
                        ].mean(),
                        "IE_ratio": breath_list[quality_breath][
                            "IE_ratio"
                        ].mean(),
                        "HR": beat_list[quality_beat]["HR"].mean(),
                        "RR": beat_list[quality_beat]["RR"].mean(),
                        
                        "VF_unfilt": breath_list[time_filter_breaths.filt]["VF"].mean(),
                        "Breath Duration_unfilt": breath_list[time_filter_breaths.filt][
                            "Breath_Cycle_Duration"
                        ].mean(),
                        "VT_unfilt": breath_list[time_filter_breaths.filt][
                            "VT__Tidal_Volume_corrected"
                        ].mean(),
                        "VT/g_unfilt": breath_list[time_filter_breaths.filt][
                            "VTpg__Tidal_Volume_per_gram_corrected"
                        ].mean(),
                        "VE_unfilt": breath_list[time_filter_breaths.filt][
                            "VE__Ventilation"
                        ].mean(),
                        "VE/g_unfilt": breath_list[time_filter_breaths.filt][
                            "VEpg__Ventilation_per_gram"
                        ].mean(),
                        "IE_ratio_unfilt": breath_list[time_filter_breaths.filt][
                            "IE_ratio"
                        ].mean(),
                        "HR_unfilt": beat_list[time_filter_beats.filt]["HR"].mean(),
                        "RR_unfilt": beat_list[time_filter_beats.filt]["RR"].mean(),
                        
                        "calibration_VT_voltage": calibration_VT,
                        "calibration_VF": calibration_VF,
                        "calibration_Breath_Cycle_Duration": calibration_Breath_Cycle_Duration,
                    }
            except Exception:
                logger.warning(
                    "unable to process calibration/baseline", exc_info=True
                )
            # %

            #% identify derrived timestamps
            recovery_dict = {}
            for cond, values in recovery_baselines.items():
                for i, v in enumerate(challenge_list):
                    try:
                        #
                        logger.info(f"Deriving Timestamps for Challenge {cond}")
                        challenge_list[i]["Challenge_Start"] = v["Prefill"]
                        challenge_list[i]["Exposure_Start"] = v["Anoxic"]
                        challenge_list[i]["Apnea_BoB(Arduino)"] = (
                            v["Apnea(Arduino)"] - settings.minimum_apnea_duration
                        )
                        challenge_list[i]["Recovery_Gas"] = v["Recovery"]
                        if i + 1 == len(challenge_list):
                            challenge_list[i]["Challenge_End"] = signals[
                                "ts"
                            ].iloc[-1]
                        else:
                            challenge_list[i]["Challenge_End"] = challenge_list[
                                i + 1
                            ]["Prefill"]
    
                        # search for breathlist apnea, filter by time, pif
                        challenge_list[i]["Apnea_BoB(B)"] = breath_list[
                            (
                                breath_list["Timestamp_Inspiration"]
                                >= challenge_list[i]["Exposure_Start"]
                            )
                            & (
                                breath_list["Timestamp_Inspiration"]
                                < challenge_list[i]["Recovery_Gas"]
                            )
                            & (
                                breath_list["Peak_Inspiratory_Flow"]
                                >= settings.minimum_PIF
                                * settings.challenge_pf_threshold_multiple
                            )
                        ]["Timestamp_Inspiration"].iloc[-1]
    
                        # search for post apnea gasps
                        # !!! make number of gasps tracked a tunable setting
                        gasp_list = breath_list[
                            (
                                breath_list["Timestamp_Inspiration"]
                                > challenge_list[i]["Apnea_BoB(B)"]
                                + settings.minimum_apnea_to_gasp_duration
                            )
                            & (
                                breath_list["Timestamp_Inspiration"]
                                < challenge_list[i]["Challenge_End"]
                            )
                            & (
                                breath_list["Peak_Inspiratory_Flow"]
                                >= settings.minimum_PIF
                                * settings.challenge_pf_threshold_multiple
                            )
                        ]
                        for g in range(10):
                            if g >= len(gasp_list):
    
    
                                challenge_list[i][f"Gasp_{g+1}"] = numpy.nan
    
                            else:
                                challenge_list[i][f"Gasp_{g+1}"] = gasp_list.iloc[
                                    g
                                ]["Timestamp_Inspiration"]
                                challenge_list[i][
                                    f"Gasp_{g+1}_Volume"
                                ] = gasp_list.iloc[g]["VT__Tidal_Volume_corrected"]
                                challenge_list[i][
                                    f"Gasp_{g+1}_Volume_per_g"
                                ] = gasp_list.iloc[g][
                                    "VTpg__Tidal_Volume_per_gram_corrected"
                                ]
    
                        # generate summary for sub-challenge components
                        for subchallenge, timing in {
                            "prefill": [
                                challenge_list[i]["Challenge_Start"],
                                challenge_list[i]["Exposure_Start"],
                            ],
                            "exposure minus delay": [
                                challenge_list[i]["Exposure_Start"]
                                + float(
                                    settings.trim_for_hypervent_phase.split(",")[0]
                                ),
                                challenge_list[i]["Apnea_BoB(Arduino)"]
                                - float(
                                    settings.trim_for_hypervent_phase.split(",")[1]
                                ),
                            ],
                            "early recovery": [
                                challenge_list[i][
                                    f"Gasp_{int(settings.early_recovery_starting_breath)}"
                                ],
                                challenge_list[i]["Challenge_End"]
                                - settings.late_recovery_start_from_end,
                            ],
                            "late recovery": [
                                challenge_list[i]["Challenge_End"]
                                - settings.late_recovery_start_from_end,
                                challenge_list[i]["Challenge_End"],
                            ],
                        }.items():
                            try:
                                sc_breaths = breath_list[
                                    (
                                        breath_list["Timestamp_Inspiration"]
                                        >= timing[0]
                                    )
                                    & (breath_list["ts_end"] < timing[1])
                                ]
                                sc_beats = beat_list[
                                    (beat_list["ts"] >= timing[0])
                                    & (beat_list["ts"] < timing[1])
                                ]
    
                                challenge_list[i][
                                    f"{subchallenge}_VF"
                                ] = sc_breaths["VF"].mean()
                                challenge_list[i][
                                    f"{subchallenge}_Breath_Duration"
                                ] = sc_breaths["Breath_Cycle_Duration"].mean()
                                challenge_list[i][
                                    f"{subchallenge}_VT"
                                ] = sc_breaths["VT__Tidal_Volume_corrected"].mean()
                                challenge_list[i][
                                    f"{subchallenge}_VT/g"
                                ] = sc_breaths[
                                    "VTpg__Tidal_Volume_per_gram_corrected"
                                ].mean()
                                challenge_list[i][
                                    f"{subchallenge}_VE"
                                ] = sc_breaths["VE__Ventilation"].mean()
                                challenge_list[i][
                                    f"{subchallenge}_VE/g"
                                ] = sc_breaths["VEpg__Ventilation_per_gram"].mean()
                                challenge_list[i][f"{subchallenge}_HR"] = sc_beats[
                                    "HR"
                                ].mean()
                                challenge_list[i][f"{subchallenge}_RR"] = sc_beats[
                                    "RR"
                                ].mean()
    
                            except Exception:
                                logger.error(
                                    f'sub-challenge "{subchallenge}" undefined for round {i+1}',
                                    exc_info=True,
                                )
                                challenge_list[i][f"{subchallenge}_VF"] = numpy.nan
                                challenge_list[i][
                                    f"{subchallenge}_Breath_Duration"
                                ] = numpy.nan
                                challenge_list[i][f"{subchallenge}_VT"] = numpy.nan
                                challenge_list[i][f"{subchallenge}_VT/g"] = numpy.nan
                                challenge_list[i][f"{subchallenge}_VE"] = numpy.nan
                                challenge_list[i][f"{subchallenge}_VE/g"] = numpy.nan
                                challenge_list[i][f"{subchallenge}_HR"] = numpy.nan
                                challenge_list[i][f"{subchallenge}_RR"] = numpy.nan
    
                        # create filter for recovered breaths and recovered beats
                        
                        breath_vf_recovery_filter = multi_filter(
                            breath_list,
                            "Timestamp_Inspiration",
                            {
                                "start": (
                                    "Timestamp_Inspiration",
                                    "ge",
                                    challenge_list[i]["Recovery_Gas"],
                                ),
                                "end": (
                                    "Timestamp_Inspiration",
                                    "l",
                                    challenge_list[i]["Challenge_End"],
                                ),
                                "VF": (
                                    "VF",
                                    "ge",
                                    settings.VF_recovery_threshold
                                    / 100
                                    * values['VF'],
                                ),
                            },
                            logger=logger,
                        )
                        
                        breath_ie_recovery_filter = multi_filter(
                            breath_list,
                            "Timestamp_Inspiration",
                            {
                                "start": (
                                    "Timestamp_Inspiration",
                                    "ge",
                                    challenge_list[i]["Recovery_Gas"],
                                ),
                                "end": (
                                    "Timestamp_Inspiration",
                                    "l",
                                    challenge_list[i]["Challenge_End"],
                                ),
                                "IE_ratio": (
                                    "IE_ratio",
                                    "ge",
                                    settings.IE_recovery_threshold
                                    / 100
                                    * values['IE_duration_ratio'],
                                ),
                            },
                            logger=logger,
                        )
    
                        beat_recovery_filter = multi_filter(
                            beat_list,
                            "ts",
                            {
                                "start": (
                                    "ts",
                                    "ge",
                                    challenge_list[i]["Recovery_Gas"],
                                ),
                                "end": (
                                    "ts",
                                    "l",
                                    challenge_list[i]["Challenge_End"],
                                ),
                                "HR": (
                                    "HR",
                                    "ge",
                                    settings.HR_recovery_threshold
                                    / 100
                                    * values['HR'],
                                ),
                            },
                            logger=logger,
                        )
    
                        # merge breath and beat filters
                        logger.info('Checking for "accumulated" recovery...')
                        recovery_filt_accum = resample_and_merge_filters(
                            breath_vf_recovery_filter,
                            beat_recovery_filter,
                            sample_interval,
                            settings.Accumulated_Recovery_Minimum_Bout,
                            logger=logger,
                        )[0]
                        logger.info('Checking for "consecutive" recovery...')
                        recovery_filt_consec = resample_and_merge_filters(
                            breath_vf_recovery_filter,
                            beat_recovery_filter,
                            sample_interval,
                            settings.Consecutive_Recovery_Threshold,
                            logger=logger,
                        )[0]
    
                        logger.info('Checking for "Breathing Only" recovery...')
                        breathing_filt_accum = resample_and_merge_filters(
                            breath_vf_recovery_filter,
                            breath_vf_recovery_filter,
                            sample_interval,
                            settings.Accumulated_Recovery_Minimum_Bout,
                            logger=logger,
                        )[0]
    
                        logger.info('Checking for "Heartbeat Only" recovery...')
                        beating_filt_accum = resample_and_merge_filters(
                            beat_recovery_filter,
                            beat_recovery_filter,
                            sample_interval,
                            settings.Accumulated_Recovery_Minimum_Bout,
                            logger=logger,
                        )[0]
                        
                        logger.info('Checking for "IE ratio" recovery...')
                        ie_filt = breath_ie_recovery_filter.rolling(
                            int(settings.IE_recovery_bout), center=False,
                        ).min().fillna(0).astype(bool).filt
    
                        # identify recovey timestamp
                        challenge_list[i][
                            "Accum Recovery Met"
                        ] = get_ts_for_accumulated_value(
                            breath_list,
                            recovery_filt_accum,
                            "Breath_Cycle_Duration",
                            settings.Accumulated_Recovery_Threshold,
                            "ts_end",
                        )
                        challenge_list[i][
                            "Accum Recovery Start"
                        ] = get_ts_for_accumulated_value(
                            breath_list,
                            recovery_filt_accum,
                            "Breath_Cycle_Duration",
                            0,
                            "Timestamp_Inspiration",
                        )
                        # !!! does the column name argument in position 3 need to be there?
                        # !!! get_ts_for_accumulated_value should be revisited
                        if ie_filt.sum() <= settings.IE_recovery_bout:
                            challenge_list[i][
                                "IE Recovery Met"
                            ] = numpy.nan  
                            challenge_list[i][
                                "IE Recovery Start"    
                            ] = numpy.nan
                        else:
                            challenge_list[i][
                                "IE Recovery Met"
                            ] = breath_list[ie_filt].Timestamp_Inspiration.iloc[
                                int(settings.IE_recovery_bout)
                            ]
                            challenge_list[i][
                                "IE Recovery Start"    
                            ] = breath_list[ie_filt].Timestamp_Inspiration.iloc[
                                0    
                            ]
                            
                        challenge_list[i][
                            "Consec Recovery Met"
                        ] = get_ts_for_accumulated_value(
                            breath_list,
                            recovery_filt_consec,
                            "Breath_Cycle_Duration",
                            settings.Consecutive_Recovery_Threshold,
                            "ts_end",
                        )
                        challenge_list[i][
                            "Consec Recovery Start"
                        ] = get_ts_for_accumulated_value(
                            breath_list,
                            recovery_filt_consec,
                            "Breath_Cycle_Duration",
                            0,
                            "Timestamp_Inspiration",
                        )
                        challenge_list[i][
                            "Breathing Only Recovery Start"
                        ] = get_ts_for_accumulated_value(
                            breath_list,
                            breathing_filt_accum,
                            "Breath_Cycle_Duration",
                            0,
                            "Timestamp_Inspiration",
                        )
    
                        challenge_list[i][
                            "Heartbeat Only Recovery Start"
                        ] = get_ts_for_accumulated_value(
                            beat_list, beating_filt_accum, "RR", 0, "ts"
                        )
    
                        print(challenge_list[i]["Accum Recovery Met"])
    
                        challenge_list[i][
                            "Latency_Apnea_BoB(B)_to_Recovery_Gas"
                        ] = numpy.subtract(
                            challenge_list[i]["Recovery_Gas"],
                            challenge_list[i]["Apnea_BoB(B)"]
                        )
                        challenge_list[i][
                            "Latency_Apnea_BoB(B)_to_1st_Gasp"
                        ] = numpy.subtract(
                            challenge_list[i]["Gasp_1"],
                            challenge_list[i]["Apnea_BoB(B)"]
                        )
                        challenge_list[i][
                            "Latency_Exposure_to_Apnea_BoB(B)"
                            ] = numpy.subtract(
                            challenge_list[i]["Apnea_BoB(B)"],
                            challenge_list[i]["Exposure_Start"]
                        )
                        challenge_list[i][
                            "Duration_Challenge_Start_to_End"
                        ] = numpy.subtract(
                            challenge_list[i]["Challenge_End"],
                            challenge_list[i]["Challenge_Start"]
                        )
                        challenge_list[i][
                            "Discrep_Apnea_BoB_B_vs_Arduino"
                        ] = numpy.subtract(
                            challenge_list[i]["Apnea_BoB(B)"],
                            challenge_list[i]["Apnea_BoB(Arduino)"]
                        )
                        challenge_list[i][
                            "Latency_1st_to_2nd_Gasp"
                        ] = numpy.subtract(
                            challenge_list[i]["Gasp_2"],
                            challenge_list[i]["Gasp_1"]
                        )
    
                        # latency from apnea/gasp to HR recovery
                        # (1st instance at accum bout)
                        challenge_list[i][
                            "Latency_Apnea_BoB(B)_to_HR_Recovery"
                        ] = numpy.subtract(
                            challenge_list[i]["Heartbeat Only Recovery Start"],
                            challenge_list[i]["Apnea_BoB(B)"]
                        )
                        # latency from apnea/gasp to VF recovery
                        # (1st instance accum bout)
                        challenge_list[i][
                            "Latency_Apnea_BoB(B)_to_VF_Recovery"
                        ] = numpy.subtract(
                            challenge_list[i]["Breathing Only Recovery Start"],
                            challenge_list[i]["Apnea_BoB(B)"]
                        )
    
                        # dual recovery
    
    
                        challenge_list[i][
                            "Latency_Apnea_BoB(B)_to_Accum_Recovery_Met"
                        ] = numpy.subtract(
                            challenge_list[i]["Accum Recovery Met"],
                            challenge_list[i]["Apnea_BoB(B)"]
                        )
                        challenge_list[i][
                            "Latency_Apnea_BoB(B)_to_Accum_Recovery_Start"
                        ] = numpy.subtract(
                            challenge_list[i]["Accum Recovery Start"],
                            challenge_list[i]["Apnea_BoB(B)"]
                        )
                        challenge_list[i][
                            "Latency_Apnea_BoB(B)_to_Consec_Recovery_Met"
                        ] = numpy.subtract(
                            challenge_list[i]["Consec Recovery Met"],
                            challenge_list[i]["Apnea_BoB(B)"]
                        )
                        challenge_list[i][
                            "Latency_Apnea_BoB(B)_to_Consec_Recovery_Start"
                        ] = numpy.subtract(
                            challenge_list[i]["Consec Recovery Start"],
                            challenge_list[i]["Apnea_BoB(B)"]
                        )
    
                        # latency between breathing and hr recovery
                        challenge_list[i][
                            "Latency_HR_Recovery_to_VF_Recovery"
                        ] = numpy.subtract(
                            challenge_list[i]["Heartbeat Only Recovery Start"],
                            challenge_list[i]["Breathing Only Recovery Start"]
                        )
    
                    except Exception:
                        logger.error(
                            f"unable to derive full timestamp set - {i}",
                            exc_info=True,
    
    
                        )
                        
                recovery_dict[cond] = challenge_list
            # %
            try:
                pre_challenge_dfs = {
                    cond:
                    pandas.DataFrame(pre_challenge_dict[cond], index=[0]) for 
                    cond in pre_challenge_dict
                }
                challenge_dfs = {
                    cond:
                    pandas.DataFrame(recovery_dict[cond]) for 
                    cond in recovery_dict
                }

                timestamps_df = pandas.DataFrame(
                    {
                        "ts": harmonized_timestamps.keys(),
                        "text": harmonized_timestamps.values(),
                    }
                )
                

                # add metadata
                meta_list = ['RUID','Line','Genotype','Weight','Sex']
                meta_dict = {}
                for k in meta_list:
                    try:
                        meta_dict[k] = breath_list.iloc[0][k]
                    except:
                        logger.warning(f'unable to locate metadata for field: {k}')
                        meta_dict[k] = 'NA'

                for cond in pre_challenge_dfs:

                    for k in meta_list:
                        pre_challenge_dfs[cond][k] = meta_dict[k]

                    try:
                        pre_challenge_dfs[cond]["Challenges"] = challenge_dfs[cond][
                            "trial_number"
                        ].max()
                    except:
                        pre_challenge_dfs[cond]["Challenges"] = numpy.nan
                    
                for cond in challenge_dfs:
                    for k in meta_list:
                        challenge_dfs[cond][k] = meta_dict[k]
                        
                    # challenge_dfs[cond][
                    #     "challenge_number_prior_to_death"
                    # ] = [
                    #         pre_challenge_dfs[cond]["Challenges"] - i for i in 
                    #         challenge_dfs[cond]["trial_number"]
                    # ]
                    
                    challenge_dfs[cond][
                        "challenge_number_prior_to_death"
                    ] = pandas.Series(
                        
                            pre_challenge_dfs[cond]["Challenges"].values[0] for i in 
                            range(challenge_dfs[cond]["trial_number"].shape[0])
                        
                    ).sub(challenge_dfs[cond]["trial_number"]).values

                for k in meta_list:
                    timestamps_df[k] = meta_dict[k]
                

                try:
                    summary_dfs = {
                        cond:challenge_dfs[cond][
                            [
                                "trial_number",
                                "challenge_number_prior_to_death",
                                "prefill_VF",
                                "prefill_VT",
                                "prefill_HR",
                                "Latency_Apnea_BoB(B)_to_Recovery_Gas",
                                "Latency_Apnea_BoB(B)_to_1st_Gasp",
                                "Latency_Exposure_to_Apnea_BoB(B)",
                                "Duration_Challenge_Start_to_End",
                                "Discrep_Apnea_BoB_B_vs_Arduino",
                                "Latency_Apnea_BoB(B)_to_HR_Recovery",
                                "Latency_Apnea_BoB(B)_to_VF_Recovery",
                                "Latency_Apnea_BoB(B)_to_Accum_Recovery_Met",
                                "Latency_HR_Recovery_to_VF_Recovery",
    
                            ]
                        ] for cond in challenge_dfs
                    }
                except Exception as e:
                    summary_dfs = {'Summary_ERROR':pandas.DataFrame()}
                    logger.error(e)

                for cond in summary_dfs:

                    for k in meta_list:
                        summary_dfs[cond][k] = breath_list.iloc[0][k]
                
                writer=pandas.ExcelWriter(os.path.join(
                    output_path,ruid+'.xlsx'
                    ),engine='xlsxwriter')
                for cond in pre_challenge_dfs:
                    pre_challenge_dfs[cond].to_excel(writer,f'Pre_Challenge_{cond}', index=False)
                for cond in challenge_dfs:
                    challenge_dfs[cond].to_excel(writer,f'Challenge_{cond}', index=False)
                timestamps_df.to_excel(writer,'Timestamps',index=False)
                for cond in summary_dfs:
                    summary_dfs[cond].to_excel(writer,f'Summary_{cond}', index=False)

                writer.close()
                logger.info("Output Saved - {}".format(output_path))

                for cond in pre_challenge_dfs:
                    if not cond in aggregate_baseline:
                        aggregate_baseline[cond] = pandas.DataFrame()
                    aggregate_baseline[cond] = pandas.concat(
                        [aggregate_baseline[cond], pre_challenge_dfs[cond]]
                    )
                for cond in challenge_dfs:
                    if not cond in aggregate_challenge:
                        aggregate_challenge[cond] = pandas.DataFrame()
                    aggregate_challenge[cond] = pandas.concat(
                        [aggregate_challenge[cond], challenge_dfs[cond]]
                    )

                aggregate_timestamps = pandas.concat(
                    [aggregate_timestamps, timestamps_df]
                )

                for cond in summary_dfs:
                    if not cond in aggregate_summary:
                        aggregate_summary[cond] = pandas.DataFrame()
                    aggregate_summary[cond] = pandas.concat(
                        [aggregate_summary[cond], summary_dfs[cond]]
                    )

                logger.info("adding to Aggregate Output")


            except Exception:
                logger.error("!!! Unable to save file !!!", exc_info=True)
        #%
        except Exception:
            logger.error(
        


                f"!!! unable to process {os.path.basename(f)}!!!",
                exc_info=True,
            )

        #%%
    writer = pandas.ExcelWriter(
        os.path.join(output_path, "Aggregate" + ".xlsx"), engine="xlsxwriter"
    )
    for cond in aggregate_baseline:
        aggregate_baseline[cond].to_excel(writer, f'Pre_Challenge_{cond}', index=False)
    for cond in aggregate_challenge:
        aggregate_challenge[cond].to_excel(writer, f'Challenge_{cond}', index=False)
    aggregate_timestamps.to_excel(writer, "Timestamps", index=False)
    for cond in aggregate_summary:
        aggregate_summary[cond].to_excel(writer, f'Summary_{cond}', index=False)
    writer.close()
    logger.info("Aggregate Output Saved - {}".format(output_path))

    # %% output as jsons
    # Aggregate_Baseline.to_json(os.path.join(output_path, "Baseline.json"))
    # Aggregate_Challenge.to_json(os.path.join(output_path, "Challenge.json"))
    # Aggregate_Timestamps.to_json(os.path.join(output_path, "Timestamps.json"))
    # Aggregate_Summary.to_json(os.path.join(output_path, "Summary.json"))



# %% define main


def main():
    # %%
    # collect filepaths
    #  collect as command line arguments
    parser = argparse.ArgumentParser(description="Automated Breath Caller")
    parser.add_argument("-s", action="append", help="Paths to signal files")
    parser.add_argument("-b", help="Path to directory containing beatlists and breathlists.")
    parser.add_argument("-o", help="Path to directory for output files")
    parser.add_argument("-p", help="Path to Analysis Parameters File")

    args, others = parser.parse_known_args()


    #  collect through gui if not set in command line
    #  to settings
    if not args.p:
        settings_path = gui_open_filename(
            {"title": "Select Analysis settings"}
        )
    else:
        settings_path = args.p
    
    if not args.s:
        signal_paths = gui_open_filenames({"title": "Select Signal Files"})
    else:
        signal_paths = args.s
    if not args.b:
        breathbeat_dir = gui_directory({"title": "Select directory containing beatlists and breathlists."})
    else:
        breathbeat_dir = args.b
    #  to output directory
    if not args.o:
        output_path = gui_directory({"title": "Select Output Location"})
    else:
        output_path = args.o

#%%
    run_sham(signal_paths, breathbeat_dir, settings_path, output_path)
#%%


if __name__ == "__main__":
    main()
