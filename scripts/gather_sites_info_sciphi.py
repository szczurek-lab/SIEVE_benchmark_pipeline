import sys
import pandas as pd
import numpy as np
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from scripts.utils import NON_WORDS_PATTERN, which_patterns, vec_pattern, get_matched_files


def get_number_of_sites(file: str) -> tuple[int, int]:
    _numOfSamples: int = 0
    _numOfMuSites: int = 0
    _numOfBackgroundSites: int = 0
    if os.path.isfile(file):
        with open(file, 'r') as fh:
            flag_1: int = 0
            flag_2: int = 0
            flag_3: int = 0
            for line in fh:
                _line = line.strip()
                if flag_1 == 0 and _line == '=numSamples=':
                    flag_1 = 1
                elif flag_1 == 1:
                    _numOfSamples = int(_line)
                    flag_1 = 2
                elif flag_2 == 0 and _line == "=mutations=":
                    flag_2 = 1
                elif flag_2 == 1:
                    if _line == "=background=":
                        flag_2 = 2
                    else:
                        _numOfMuSites += 1
                elif flag_2 == 2 and flag_3 == 0:
                    flag_3 = 1
                elif flag_3 == 1:
                    _comp = _line.split()
                    for i in range(int(len(_comp) / 2)):
                        _numOfBackgroundSites += int(_comp[2 * i + 1])
                    _numOfBackgroundSites /= _numOfSamples

                    if abs(_numOfBackgroundSites - int(_numOfBackgroundSites)) > 1e-8:
                        print(f'Warning: the number of background sites is not an integer: {_numOfBackgroundSites}')
                    else:
                        _numOfBackgroundSites = int(_numOfBackgroundSites)

                    flag_3 = 2
                elif flag_1 == 2 and flag_2 == 2 and flag_3 == 2:
                    break
    return(_numOfMuSites, _numOfBackgroundSites)


search_for_tool_setup_flag = False
if snakemake.params['toolSetup'] is None or type(snakemake.params['toolSetup']) == str:
    tool_setups = [snakemake.params['toolSetup']]
elif type(snakemake.params['toolSetup']) == list:
    search_for_tool_setup_flag = True
    tool_setups = sorted(snakemake.params['toolSetup'])
else:
    raise ValueError(f'toolSetup must be a string or a list of strings ({snakemake.params["toolSetup"]} was given).')

fine_tune_type_flag = False
if snakemake.params['fineTuneType'] is None:
    fine_tune_types = np.repeat(np.nan, len(snakemake.input))
elif len(snakemake.params['fineTuneType']) == 1:
    fine_tune_types = np.repeat(snakemake.params['fineTuneType'], len(snakemake.input))
else:
    fine_tune_types = []
    fine_tune_type_flag = True

data_type_flag = False
if snakemake.params['dataType'] is None:
    data_types = np.repeat(np.nan, len(snakemake.input))
elif len(snakemake.params['dataType']) == 1:
    data_types = np.repeat(snakemake.params['dataType'], len(snakemake.input))
else:
    data_types = []
    data_type_flag = True


tool_setup = []
dataset = []
input_file_names = []
num_mu_sites = []
num_background_sites = []

for __tool_setup in tool_setups:
    for dataset_name in sorted(snakemake.params['datasetNames']):
        if search_for_tool_setup_flag:
            __pat = ''.join(
                [f'(?=.*{i})' for i in vec_pattern([__tool_setup, dataset_name], NON_WORDS_PATTERN)]
            )
        else:
            __pat = vec_pattern(dataset_name, NON_WORDS_PATTERN)

        __idx, __files = get_matched_files(
            snakemake.input,
            __pat
        )
        cnt = np.sum(__idx)

        if cnt == 0:
            print(f"Warning: no matches found for {__pat} in {snakemake.input}.")

        tool_setup.extend(np.repeat(__tool_setup, cnt))
        dataset.extend(np.repeat(dataset_name, cnt))
        input_file_names.extend(__files)

        for file in __files:
            __num_mu_sites, __num_background_sites = get_number_of_sites(file)
            num_mu_sites.append(__num_mu_sites)
            num_background_sites.append(__num_background_sites)

        if fine_tune_type_flag:
            fine_tune_type_idx = which_patterns(
                __files,
                vec_pattern(
                    snakemake.params['fineTuneType'],
                    NON_WORDS_PATTERN
                )
            )
            fine_tune_types.extend([snakemake.params['fineTuneType'][j] for i in fine_tune_type_idx for j in i if j >= 0])
        
        if data_type_flag:
            data_type_idx = which_patterns(
                __files,
                vec_pattern(
                    snakemake.params['dataType'],
                    NON_WORDS_PATTERN
                )
            )
            data_types.extend([snakemake.params['dataType'][j] for i in data_type_idx for j in i if j >= 0])


sites_info_collection = pd.DataFrame(
    {
        'cell_num': snakemake.params['cellNum'],
        'coverage_mean': snakemake.params["covMean"],
        'coverage_variance': snakemake.params["covVariance"],
        'eff_seq_err_rate': snakemake.params['effSeqErrRate'],
        'ado_rate': snakemake.params['adoRate'],
        'deletion_rate': snakemake.params['deletionRate'],
        'insertion_rate': snakemake.params['insertionRate'],
        'dataset': pd.Series(dataset),
        'tool': snakemake.params['tool'],
        'snv_type': snakemake.params['snvType'],
        'tool_setup': pd.Series(tool_setup),
        'fine_tune_type': pd.Series(fine_tune_types),
        'data_type': pd.Series(data_types),
        'num_variant_sites': pd.Series(num_mu_sites),
        'num_background_sites': pd.Series(num_background_sites)
    }
)

sites_info_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
