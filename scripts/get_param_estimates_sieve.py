from typing import Tuple
import pandas as pd
from numpy import nan
import os


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


def parse_params(params_file: str, estimate_type: str) -> Tuple[float, float, float, float, float, float, float, float]:
    if os.path.isfile(params_file):
        deletion = insertion = allelic_seq_cov = allelic_seq_var = eff_seq_err_rate = ado = sc_1 = sc_2 = gamma = pop = nan
        reach_beginning = reach_ending = False
        with open(params_file, 'r') as fh:
            for line in fh:
                if reach_ending is True:
                    break

                if reach_beginning is True:
                    if not line.strip().startswith('#'):
                        comp = line.strip().split('\t')
                        if comp[1].strip() in estimate_type:
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


def match_names(file_name: str, template: str, dataset: str) -> bool:
    comp = file_name.strip().split('/')
    found_template: bool = False
    found_dataset: bool = False

    for item in comp:
        if item.strip() == template:
            found_template = True
        elif item.strip() == dataset:
            found_dataset = True
    return found_template and found_dataset


for template in sorted(snakemake.params['configTemplates']):
    for dataset_name in sorted(snakemake.params['datasetNames']):
        tool_setup.append(template)
        dataset.append(dataset_name)
        for params_log in snakemake.input:
            if match_names(params_log, template, dataset_name) == True:
                if snakemake.params['estimatesType'] == 'all':
                    for root, dirs, files in os.walk(os.path.dirname(params_log), topdown=True):
                        found = False
                        for file in files:
                            if '.vcf' in file:
                                found = True
                                start_index = file.index(
                                    dataset_name) + len(dataset_name) + 1
                                estimate_types.append(
                                    file[start_index:file.index('.vcf')])
                                break
                        if found is True:
                            break
                else:
                    estimate_types.append(snakemake.params['estimatesType'])

                deletion, insertion, allelic_seq_cov, allelic_seq_var, eff_seq_err_rate, ado, sc_1, sc_2, gamma, pop = parse_params(
                    params_log, estimate_types[-1])
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

                break

params_collection = pd.DataFrame({'cell_num': snakemake.params['cellNum'],
                                  'coverage_mean': snakemake.params["covMean"],
                                  'coverage_variance': snakemake.params["covVariance"],
                                  'dataset': pd.Series(dataset),
                                  'tool': snakemake.params['tool'],
                                  'snv_type': snakemake.params['snvType'],
                                  'tool_setup': pd.Series(tool_setup), 
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
                                  })

params_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
