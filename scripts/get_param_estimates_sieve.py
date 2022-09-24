import sys
import pandas as pd
import numpy as np
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from scripts.utils import NON_WORDS_PATTERN, which_patterns, vec_pattern, get_matched_files


def parse_params(
    params_file: str,
    estimate_type: str
) -> tuple[float, float, float, float, float, float, float, float]:
    if os.path.isfile(params_file):
        deletion = insertion = allelic_seq_cov = allelic_seq_var = eff_seq_err_rate = ado = sc_1 = sc_2 = gamma = pop = np.nan
        reach_beginning = reach_ending = False
        with open(params_file, 'r') as fh:
            for line in fh:
                if reach_ending is True:
                    break

                if reach_beginning is True:
                    if not line.strip().startswith('#'):
                        comp = line.strip().split('\t')
                        if comp[1].strip() in estimate_type or estimate_type == 'mode' and comp[1].strip() not in ['mean', 'median']:
                            if comp[0] == 'deletionRate':
                                deletion = float(comp[2].strip())
                            elif comp[0] == 'insertionRate':
                                insertion = float(comp[2].strip())
                            elif comp[0] == 'allelicSeqCov':
                                allelic_seq_cov = float(comp[2].strip())
                            elif comp[0]== 'allelicSeqCovRawVar':
                                allelic_seq_var = float(comp[2].strip())
                            elif comp[0] == 'effSeqErrRate':
                                eff_seq_err_rate = float(comp[2].strip())
                            elif comp[0] == 'adoRate':
                                ado = float(comp[2].strip())
                            elif comp[0] == 'shapeCtrl1':
                                sc_1 = float(comp[2].strip())
                            elif comp[0] == 'shapeCtrl2':
                                sc_2 = float(comp[2].strip())
                            elif comp[0] == 'gammaShape':
                                gamma = float(comp[2].strip())
                            elif comp[0] == 'ePopSize':
                                pop = float(comp[2].strip())
                    else:
                        reach_beginning = False
                        reach_ending = True
                else:
                    if line.startswith('#MCMC samples'):
                        reach_beginning = True

        return deletion, insertion, allelic_seq_cov, allelic_seq_var, eff_seq_err_rate, ado, sc_1, sc_2, gamma, pop


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


candidate_estimate_types = ['mean', 'median', 'mode']

tool_setup = []
dataset = []
tool_setup = []
deletion_rates = []
insertion_rates = []
allelic_seq_covs = []
allelic_seq_vars = []
eff_seq_err_rates = []
ado_rates = []
shape_ctrler_1 = []
shape_ctrler_2 = []
gamma_shape = []
pop_size = []
estimate_types = []


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

        for __params in __inferred_params:
            if snakemake.params['estimatesType'] == 'all':
                for _, _, files in os.walk(os.path.dirname(__params), topdown=True):
                    estimate_type_idx = which_patterns(
                        files,
                        [f'(?=.*{i})(?=.*\.vcf$)' for i in vec_pattern(candidate_estimate_types, NON_WORDS_PATTERN)]
                    )
                    __ret = [candidate_estimate_types[j] for i in estimate_type_idx for j in i if j >= 0]
                    
                    if len(__ret) > 0:
                        estimate_types.extend(__ret)
                        break
            else:
                estimate_types.append(snakemake.params['estimatesType'])

            deletion, insertion, allelic_seq_cov, allelic_seq_var, eff_seq_err_rate, ado, sc_1, sc_2, gamma, pop = parse_params(
                __params, estimate_types[-1])
            deletion_rates.append(deletion)
            insertion_rates.append(insertion)
            allelic_seq_covs.append(allelic_seq_cov)
            allelic_seq_vars.append(allelic_seq_var)
            eff_seq_err_rates.append(eff_seq_err_rate)
            ado_rates.append(ado)
            shape_ctrler_1.append(sc_1)
            shape_ctrler_2.append(sc_2)
            gamma_shape.append(gamma)
            pop_size.append(pop)

        if fine_tune_type_flag:
            fine_tune_type_idx = which_patterns(
                __inferred_params,
                vec_pattern(
                    snakemake.params['fineTuneType'],
                    NON_WORDS_PATTERN
                )
            )
            fine_tune_types.extend([snakemake.params['fineTuneType'][j] for i in fine_tune_type_idx for j in i if j >= 0])
        
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
        'tool_setup': pd.Series(tool_setup),
        'fine_tune_type': pd.Series(fine_tune_types),
        'data_type': pd.Series(data_types),
        'estimates_type': pd.Series(estimate_types),
        'deletion_rate': pd.Series(deletion_rates),
        'insertion_rate': pd.Series(insertion_rates),
        'allelic_seq_cov': pd.Series(allelic_seq_covs),
        'allelic_seq_cov_raw_var': pd.Series(allelic_seq_vars),
        'eff_seq_err_rate': pd.Series(eff_seq_err_rates),
        'ado_rate': pd.Series(ado_rates),
        'wild_overdispersion': pd.Series(shape_ctrler_1),
        'alternative_overdispersion': pd.Series(shape_ctrler_2),
        'gamma_shape': pd.Series(gamma_shape),
        'population_size': pd.Series(pop_size)
    }
)

params_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
