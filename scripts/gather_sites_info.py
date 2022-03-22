from typing import List
import pandas as pd

def get_number_of_sites(files: List[str]) -> pd.DataFrame:
    _template: List[str] = list()
    _dataset: List[int] = list()
    _numOfMuSites: List[int] = list()
    _numOfBackgroundSites: List[int] = list()

    for _file in files:
        _comp: List[str] = _file.strip().split("/")
        _template.append(_comp[-3])
        _dataset.append(int(_comp[-2]))
        with open(_file, 'r') as fh:
            flag_1: int = 0
            flag_2: int = 0
            for line in fh:
                _line = line.strip()
                if flag_1 == 0 and _line == "=numCandidateMutatedSites=":
                    flag_1 = 1
                elif flag_1 == 1:
                    _numOfMuSites.append(int(_line))
                    flag_1 = 2
                elif flag_2 == 0 and _line == "=numBackgroundSites=":
                    flag_2 = 1
                elif flag_2 == 1:
                    _numOfBackgroundSites.append(int(_line))
                    flag_2 = 2
                elif flag_1 == 2 and flag_2 == 2:
                    break
    return(pd.DataFrame({
        'template1': _template,
        'dataset': _dataset,
        'num_candidate_mutated_sites': _numOfMuSites,
        'num_background_sites': _numOfBackgroundSites
    }))


pd.read_csv(
    snakemake.input['treeInfo'],
    sep='\t'
    ).query(
        'tool.str.match(r"(?i)sieve") & tool_setup.str.contains(r"-")'
    ).drop(
        columns=['inferred_tree', 'real_tree']
    ).assign(
        template1=lambda x: x.tool_setup.str.split('-').str[0]
    ).merge(
        get_number_of_sites(snakemake.input['data']),
        on=['template1', 'dataset']
    ).drop(
        columns='template1'
    ).to_csv(
        snakemake.output[0], 
        sep='\t', 
        na_rep='NA', 
        index=False
    )
