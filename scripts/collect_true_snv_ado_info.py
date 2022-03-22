from typing import List, Tuple
import os


def _get_cell_names(input: str) -> List[str]:
    if os.path.isfile(input):
        results = []
        with open(input, 'r') as fh:
            for line in fh:
                results.append(line.strip())
        return results

def _get_loci_info(input: str) -> List[Tuple[str, int]]:
    if os.path.isfile(input):
        results = []
        with open(input, 'r') as fh:
            for line in fh:
                comp = line.strip().split()
                results.append((comp[0].strip(), int(comp[1].strip())))
        return results

def _get_genotypes(input: str) -> List[List[str]]:
    if os.path.isfile(input):
        results = []
        with open(input, 'r') as fh:
            for line in fh:
                results.append(line.strip().split())
        return results

def _collect_genotype_ado(
    input: str,
    cell_names: List[str],
    snv_sites: List[Tuple[str, int]],
    genotypes: List[List[str]]
):
    if os.path.isfile(input):
        with open(input, 'r') as fh:
            for line in fh:
                if not line.strip().startswith('#'):
                    comp: List[str] = line.strip().split()

                    if (comp[0], int(comp[1])) not in snv_sites:
                        continue

                    _snv_index: int = snv_sites.index((comp[0], int(comp[1])))
                    
                    for _cell_index in range(len(cell_names)):
                        _gp: str = genotypes[_snv_index][_cell_index]

                        if _gp == '0/1':
                            if int(comp[2 + _cell_index]) == 0:
                                without_ado_hetero.append((_snv_index, _cell_index))
                            else:
                                with_ado_hetero.append((_snv_index, _cell_index))
                        elif _gp == '1/1' or _gp == '1/2':
                            if int(comp[2 + _cell_index]) == 0:
                                without_ado_homo.append((_snv_index, _cell_index))
                            else:
                                with_ado_homo.append((_snv_index, _cell_index))

def _write_output(
    data: List[Tuple[int, int]],
    output: str
):
    with open(output, 'w') as fh:
        for _item in data:
            fh.write(str(_item[0]))
            fh.write('\t')
            fh.write(str(_item[1]))
            fh.write('\n')


with_ado_hetero: List[Tuple[int, int]] = []
without_ado_hetero: List[Tuple[int, int]] = []
with_ado_homo: List[Tuple[int, int]] = []
without_ado_homo: List[Tuple[int, int]] = []

cell_names: List[str] = _get_cell_names(snakemake.input['cellNames'])
snv_sites: List[Tuple[str, int]] = _get_loci_info(snakemake.input['snvSitesNames'])
genotypes: List[List[str]] = _get_genotypes(snakemake.input['trueGenotypes'])

_collect_genotype_ado(
    snakemake.input['trueAdoStates'],
    cell_names,
    snv_sites,
    genotypes
)

_write_output(with_ado_hetero, snakemake.output['withAdoHetero'])
_write_output(without_ado_hetero, snakemake.output['withoutAdoHetero'])
_write_output(with_ado_homo, snakemake.output['withAdoHomo'])
_write_output(without_ado_homo, snakemake.output['withoutAdoHomo'])
