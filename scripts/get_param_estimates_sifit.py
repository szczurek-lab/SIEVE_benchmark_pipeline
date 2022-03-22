from typing import Tuple
import pandas as pd
import os


dataset = []
deletion_rates = []
loh_rates = []
ado_rates = []


def parse_params(params_file: str) -> Tuple[float, float, float]:
    if os.path.isfile(params_file):
        deletion = loh = ado = 0.0
        with open(params_file, 'r') as fh:
            for line in fh:
                if line.startswith('best false negative rate'):
                    ado = float(line.split('=')[1].strip())
                elif line.startswith('best deletion parameter'):
                    deletion = float(line.split('=')[1].strip())
                elif line.startswith('best LOH parameter'):
                    loh = float(line.split('=')[1].strip())
        return deletion, loh, ado


for dataset_name in sorted(snakemake.params['datasetNames']):
    dataset.append(dataset_name)
    for params_log in snakemake.input:
        if dataset_name in params_log:
            deletion, loh, ado = parse_params(params_log)
            deletion_rates.append(deletion)
            loh_rates.append(loh)
            ado_rates.append(ado)

params_collection = pd.DataFrame({'cell_num': snakemake.params['cellNum'],
                                  'coverage_mean': snakemake.params["covMean"],
                                  'coverage_variance': snakemake.params["covVariance"],
                                  'dataset': pd.Series(dataset),
                                  'tool': snakemake.params['tool'],
                                  'snv_type': snakemake.params['snvType'],
                                  'tool_setup': snakemake.params['toolSetup'], 
                                  'deletion_rate': pd.Series(deletion_rates),
                                  'loh_rate': pd.Series(loh_rates),
                                  'ado_rate': pd.Series(ado_rates)
                                  })

params_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
