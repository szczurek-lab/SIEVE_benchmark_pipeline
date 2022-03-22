from typing import Tuple
import pandas as pd
import os


dataset = []
eff_seq_err_rates = []
ado_rates = []
wild_overdispersion = []
alternative_overdispersion = []
zygosity_rates = []


def parse_params(params_file: str) -> Tuple[float, float, float, float, float]:
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


for dataset_name in sorted(snakemake.params['datasetNames']):
    dataset.append(dataset_name)
    for params_log in snakemake.input:
        if dataset_name in params_log:
            eff_seq_err_rate, ado, wild_overdis, alt_overdis, zyg = parse_params(
                params_log)
            eff_seq_err_rates.append(eff_seq_err_rate)
            ado_rates.append(ado)
            wild_overdispersion.append(wild_overdis)
            alternative_overdispersion.append(alt_overdis)
            zygosity_rates.append(zyg)

params_collection = pd.DataFrame({'cell_num': snakemake.params['cellNum'],
                                  'coverage_mean': snakemake.params["covMean"],
                                  'coverage_variance': snakemake.params["covVariance"],
                                  'dataset': pd.Series(dataset),
                                  'tool': snakemake.params['tool'],
                                  'snv_type': snakemake.params['snvType'],
                                  'tool_setup': snakemake.params['toolSetup'], 
                                  'eff_seq_err_rate': pd.Series(eff_seq_err_rates),
                                  'ado_rate': pd.Series(ado_rates),
                                  'wild_overdispersion': pd.Series(wild_overdispersion),
                                  'alternative_overdispersion': pd.Series(alternative_overdispersion),
                                  'zygosity_rate': pd.Series(zygosity_rates)
                                  })

params_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
