import sys
import pandas as pd
import numpy as np
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from scripts.utils import NON_WORDS_PATTERN, which_patterns, vec_pattern, get_matched_files


def parse_params(params_file: str) -> tuple[float, float, float, float, float]:
    eff_seq_err_rate = ado = wild_overdis = alt_overdis = zyg = 0.0
    wild_alpha = wild_beta = alt_alpha = alt_beta = -1.0
    if os.path.isfile(params_file):
        with open(params_file, 'r') as fh:
            for line in fh:
                if line.startswith('wildAlpha:'):
                    wild_alpha = float(line.strip().split(
                        '\t')[snakemake.params['estimatesType'] + 1])
                elif line.startswith('wildBeta:'):
                    wild_beta = float(line.strip().split(
                        '\t')[snakemake.params['estimatesType'] + 1])
                elif line.startswith('altAlpha:'):
                    alt_alpha = float(line.strip().split(
                        '\t')[snakemake.params['estimatesType'] + 1])
                elif line.startswith('altBeta:'):
                    alt_beta = float(line.strip().split(
                        '\t')[snakemake.params['estimatesType'] + 1])
                elif line.startswith('mu'):
                    ado = 1 - float(line.strip().split('\t')[snakemake.params['estimatesType'] + 1])
                elif line.startswith('nu'):
                    zyg = float(line.strip().split('\t')[
                                snakemake.params['estimatesType'] + 1])
        eff_seq_err_rate = wild_alpha / (wild_alpha + wild_beta)
        wild_overdis = wild_alpha + wild_beta
        alt_overdis = alt_alpha + alt_beta
        return eff_seq_err_rate, ado, wild_overdis, alt_overdis, zyg


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
eff_seq_err_rates = []
ado_rates = []
wild_overdispersion = []
alternative_overdispersion = []
zygosity_rates = []


for __tool_setup in tool_setups:
    for dataset_name in sorted(snakemake.params['datasetNames']):
        if search_for_tool_setup_flag:
            __pat = ''.join(
                [f'(?=.*{i})' for i in vec_pattern([__tool_setup, dataset_name], NON_WORDS_PATTERN)]
            )
        else:
            __pat = vec_pattern(dataset_name, NON_WORDS_PATTERN)

        inferred_idx, __inferred_params = get_matched_files(
            snakemake.input,
            __pat
        )
        cnt = np.sum(inferred_idx)

        if cnt == 0:
            print(f"Warning: no matches found for {__pat} in {snakemake.input}.")

        tool_setup.extend(np.repeat(__tool_setup, cnt))
        dataset.extend(np.repeat(dataset_name, cnt))

        if fine_tune_type_flag:
            fine_tune_type_idx = which_patterns(
                __inferred_params,
                vec_pattern(
                    snakemake.params['fineTuneType'],
                    NON_WORDS_PATTERN
                )
            )
            fine_tune_types.extend([snakemake.params['fineTuneType'][j] for i in fine_tune_type_idx for j in i if j >= 0])

        for __params in __inferred_params:
            eff_seq_err_rate, ado, wild_overdis, alt_overdis, zyg = parse_params(__params)
            eff_seq_err_rates.append(eff_seq_err_rate)
            ado_rates.append(ado)
            wild_overdispersion.append(wild_overdis)
            alternative_overdispersion.append(alt_overdis)
            zygosity_rates.append(zyg)
        
        if data_type_flag:
            data_type_idx = which_patterns(
                __inferred_params,
                vec_pattern(
                    snakemake.params['dataType'],
                    NON_WORDS_PATTERN
                )
            )
            data_types.extend([snakemake.params['dataType'][j] for i in data_type_idx for j in i if j >= 0])


params_collection = pd.DataFrame(
    {
        'cell_num': snakemake.params['cellNum'],
        'coverage_mean': snakemake.params["covMean"],
        'coverage_variance': snakemake.params["covVariance"],
        'dataset': pd.Series(dataset),
        'tool': snakemake.params['tool'],
        'snv_type': snakemake.params['snvType'],
        'tool_setup': snakemake.params['toolSetup'],
        'fine_tune_type': pd.Series(fine_tune_types),
        'data_type': pd.Series(data_types), 
        'eff_seq_err_rate': pd.Series(eff_seq_err_rates),
        'ado_rate': pd.Series(ado_rates),
        'wild_overdispersion': pd.Series(wild_overdispersion),
        'alternative_overdispersion': pd.Series(alternative_overdispersion),
        'zygosity_rate': pd.Series(zygosity_rates)
    }
)

params_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
