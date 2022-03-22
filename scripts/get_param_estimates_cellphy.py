from typing import Tuple
import pandas as pd
from numpy import nan
import os


dataset = []
eff_seq_err_rates = []
ado_rates = []
gamma_shape = []


def parse_params(params_file: str) -> Tuple[float, float, float]:
    if os.path.isfile(params_file):
        eff_seq_err_rate = ado = gamma = nan
        with open(params_file, 'r') as fh:
            for line in fh:
                if 'ERR_P17' in line:
                    index = line.index('ERR_P17')
                    start_index = end_index = index
                    for i in range(index, len(line)):
                        if line[i] == '{':
                            start_index = i + 1
                        elif line[i] == '}':
                            end_index = i
                            break
                    if start_index >= end_index:
                        raise ValueError(
                            'Error! Invalid format of error parameters in this file: ' + params_file)
                    else:
                        params = line[start_index:end_index].strip().split('/')
                        eff_seq_err_rate = float(params[0].strip())
                        ado = float(params[1].strip())

                if 'G4m' in line:
                    index = line.index('G4m')
                    start_index = end_index = index
                    for i in range(index, len(line)):
                        if line[i] == '{':
                            start_index = i + 1
                        elif line[i] == '}':
                            end_index = i
                            break
                    if start_index >= end_index:
                        raise ValueError(
                            'Error! Invalid format of error parameters in this file: ' + params_file)
                    else:
                        gamma = float(line[start_index:end_index].strip())

                return eff_seq_err_rate, ado, gamma


for dataset_name in sorted(snakemake.params['datasetNames']):
    dataset.append(dataset_name)
    for params_log in snakemake.input:
        if dataset_name in params_log:
            eff_seq_err_rate, ado, gamma = parse_params(params_log)
            eff_seq_err_rates.append(eff_seq_err_rate)
            ado_rates.append(ado)
            gamma_shape.append(gamma)

params_collection = pd.DataFrame({'cell_num': snakemake.params['cellNum'],
                                  'coverage_mean': snakemake.params["covMean"],
                                  'coverage_variance': snakemake.params["covVariance"],
                                  'dataset': pd.Series(dataset),
                                  'tool': snakemake.params['tool'],
                                  'snv_type': snakemake.params['snvType'],
                                  'tool_setup': snakemake.params['toolSetup'], 
                                  'eff_seq_err_rate': pd.Series(eff_seq_err_rates),
                                  'ado_rate': pd.Series(ado_rates),
                                  'gamma_shape': pd.Series(gamma_shape)
                                  })

params_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
