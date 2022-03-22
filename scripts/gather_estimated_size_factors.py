from typing import List
import pandas as pd
import numpy as np
import os


def _get_specific_file(files: List[str], key1: str, key2: str = '') -> str:
    for file in files:
        if key2 == '':
            if key1 in file:
                return file
        else:
            if key1 in file and key2 in file:
                return file
    raise ValueError('Error! ' + key1 + ' and ' + key2 + ' are not found.')

def _get_content_files(file: str) -> List[str]:
    if os.path.isfile(file):
        results = []
        root = os.path.abspath(os.path.dirname(file))
        with open(file, 'r') as fh:
            for line in fh:
                content = os.path.join(root, line.strip())
                if os.path.isfile(content):
                    results.append(content)
        return results

def _get_cell_names(input: str) -> List[str]:
    if os.path.isfile(input):
        results = []
        with open(input, 'r') as fh:
            for line in fh:
                results.append(line.strip())
        return results

def _get_values(input: str) -> List[float]:
    if os.path.isfile(input):
        results = []
        with open(input, 'r') as fh:
            for line in fh:
                results.append(float(line.strip()))
        return results


tool_setup = []
dataset = []
cell_names = []
true_size_factors = []
size_factors = []
difference = []
difference_in_percentage = []

_true_cell_names = _get_cell_names(snakemake.params['trueCellNames'])

for sieve_run_template in sorted(snakemake.params['sieveRunTemplates']):
        for dataset_name in sorted(snakemake.params['datasetNames']):
            tool_setup.extend([sieve_run_template for _ in range(len(_true_cell_names))])
            dataset.extend([dataset_name for _ in range(len(_true_cell_names))])
            cell_names.extend(_true_cell_names)

            files_list = _get_content_files(_get_specific_file(snakemake.input, sieve_run_template, dataset_name))
            
            _cell_names = _get_cell_names(_get_specific_file(files_list, '.cell_names'))
            if len(_true_cell_names) != len(_cell_names):
                raise ValueError('Error! Unmatched number of cells in ' + snakemake.params['trueCellNames'] + ' and ' + _get_specific_file(files_list, '.cell_names'))
            _cell_index_in_true = []
            for _cell_name in _cell_names:
                _cell_index_in_true.append(_true_cell_names.index(_cell_name))

            _true_size_factors = _get_values(os.path.abspath(snakemake.params['trueSizeFactorsPre'] + '.' + dataset_name))
            true_size_factors.extend(_true_size_factors)
            _size_factors = _get_values(_get_specific_file(files_list, '.size_factors'))
            _ordered_size_factors = [_size_factors[i] for i in _cell_index_in_true]
            size_factors.extend(_ordered_size_factors)
            
            _difference = np.array(_ordered_size_factors, dtype=float) - np.array(_true_size_factors, dtype=float)
            _difference_in_percentage = _difference / np.array(_true_size_factors, dtype=float)
            
            difference.extend(list(_difference))
            difference_in_percentage.extend(list(_difference_in_percentage))

size_factors_collection = pd.DataFrame({
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
    'cell_names': pd.Series(cell_names),
    'true_size_factors': pd.Series(true_size_factors),
    'estimated_size_factors': pd.Series(size_factors),
    'difference': pd.Series(difference),
    'difference_in_percentage': pd.Series(difference_in_percentage)
    })

size_factors_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
