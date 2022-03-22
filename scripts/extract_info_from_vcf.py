import argparse
import os
from typing import Tuple, Dict, List, Union


def parse_system_args() -> Tuple[argparse.Namespace, Dict[str, str]]:
    sys_args_parser = argparse.ArgumentParser(description='Extract information from vcf file.')
    sys_args_parser.add_argument('input', help='a vcf file or a file list containing a vcf file (for sieve2)', type=str)
    sys_args_parser.add_argument('--mode', help='which type of the input vcf file is? For monovar and sieve, '
                                                'a standard vcf file is input; for sciphi, the input vcf file has no '
                                                'PL filed. For sieve1, --input would be vcf; for sieve2, --input '
                                                'should be a file list containing a vcf file name.',
                                 type=str, choices=['monovar', 'sciphi', 'sieve1', "sieve2"])
    sys_args_parser.add_argument('--prefix', help='prefix of output files; if not provided, prefix of input file will '
                                                  'be used', type=str)
    sys_args_parser.add_argument('--probs', help='a .probs file should be specified when --mode is sciphi', type=str)
    sys_args_parser.add_argument('--sifit',
                                 help='save input data to SIFIT: binary (b), ternary (t), or not saving (n)?',
                                 type=str, choices=['n', 'b', 't'], default='n')
    sys_args = sys_args_parser.parse_args()

    if (sys_args.mode == 'monovar' or sys_args.mode == 'sieve') and sys_args.probs is not None:
        raise ValueError('Error! When in mode \'monovar\' or \'sieve\', no --probs file should be specified.')
    if sys_args.mode == 'sciphi' and sys_args.probs is None:
        raise ValueError('Error! When in mode \'sciphi\', --probs file should be specified.')
    if sys_args.mode == 'sieve' and sys_args.sifit == 'n':
        raise ValueError('Error! When in mode \'sieve\', --sifit should be set to \'b\' or \'t\'.')
    if not os.path.isfile(sys_args.input):
        raise ValueError('Error! File ' + sys_args.input + ' does not exist.')

    if sys_args.prefix is not None:
        prefix = sys_args.prefix
    else:
        prefix = os.path.splitext(sys_args.input)[0]

    output_files = {'cell_names': prefix + '.cell_names',
                    'loci_info': prefix + '.loci_info',
                    'genotypes': prefix + '.genotypes',
                    'ternary': prefix + '.ternary',
                    'genotype_probs': prefix + '.genotype_probs',
                    'sifit_data': prefix + '.sifit_data',
                    'sifit_cell_names': prefix + '.sifit_cell_names'}
    return sys_args, output_files


def is_all_at_least_one(values: List[str]) -> bool:
    if values is None or len(values) == 0:
        return False

    for i in values:
        try:
            value = int(i)
        except ValueError:
            return False

        if value < 1:
            return False

    return True


def is_zero_and_more(values: List[str]) -> bool:
    if values is None or len(values) == 0:
        return False

    results: List[int] = []
    for i in values:
        try:
            value = int(i)
        except ValueError:
            return False

        results.append(value)

    results.sort()

    if results[0] == 0 and results[-1] > 0:
        return True
    else:
        return False


def genotype2ternary(genotype: str) -> int:
    if genotype == '0/0' or genotype == '0|0' or ('.' in genotype and '0' in genotype):
        return 0
    elif is_zero_and_more(genotype.split('|') if '|' in genotype else genotype.split('/')):
        return 1
    elif is_all_at_least_one(genotype.split('|') if '|' in genotype else genotype.split('/')):
        return 2
    elif genotype == '.' or genotype == './.' or genotype == '.|.':
        return 3
    else:
        raise ValueError(
            'Error! Unsupported genotype: ' +
            genotype +
            '. Only support '
            '\'0/0\', \'0|0\', \'./0\', \'.|0\', \'0/.\', \'0|.\', '
            'genotypes with at least one allele mutated, '
            '\'.\', \'./.\', \'.|.\'.'
        )


def get_genotypes_probs(gt00: int, gt01: int, gt11: int) -> Tuple[float, float, float]:
    gt00_f = pow(10, -1 * gt00 / 10)
    gt01_f = pow(10, -1 * gt01 / 10)
    gt11_f = pow(10, -1 * gt11 / 10)

    prob_sum = gt00_f + gt01_f + gt11_f
    return gt00_f / prob_sum, gt01_f / prob_sum, gt11_f / prob_sum


def get_vcf_file(input_file: str) -> str:
    with open(input_file, 'r') as fh:
        for line in fh:
            if line.strip().endswith('.vcf'):
                return os.path.join(os.path.dirname(input_file), line.strip())
    raise ValueError('Error! File ' + input_file + ' does not contain a vcf file.')


def parse_vcf_file(input_file: str, mode: str, probs_file: str) -> Tuple[
    List[str], List[Tuple[str, int]], Dict[str, List[str]], Dict[str, List[int]], Dict[
        str, List[Union[Tuple[float, float, float], Tuple[float, float]]]]]:
    tmp_input_file: str = input_file

    if mode == "sieve2":
        tmp_input_file = get_vcf_file(input_file)

    if os.path.isfile(tmp_input_file):
        cell_names = []
        loci_info = []
        genotypes = {}
        ternary = {}
        genotype_probs = {}

        chrom_index = pos_index = gt_index = pl_index = -1
        cell_info_index = 0

        with open(tmp_input_file, 'r') as fh:
            for line in fh:
                if line.strip().startswith('#'):
                    if line.strip().startswith('#CHROM'):
                        comp = line.strip().split('\t')
                        reach_cells = False
                        for item in comp:
                            if reach_cells is True:
                                cell_name = item.strip().replace('.bam', '')
                                cell_names.append(cell_name)
                                genotypes[cell_name] = []
                                ternary[cell_name] = []
                                genotype_probs[cell_name] = []
                            else:
                                _item = item.strip().upper()
                                if 'CHROM' in _item:
                                    chrom_index = cell_info_index
                                elif _item == 'POS':
                                    pos_index = cell_info_index
                                elif _item == 'FORMAT':
                                    reach_cells = True
                                cell_info_index += 1
                else:
                    comp = line.strip().split('\t')
                    loci_info.append((comp[chrom_index].strip(), int(comp[pos_index].strip())))

                    if mode.lower() == 'monovar' and pl_index == -1:
                        pl_index = 0
                        _comp = comp[cell_info_index - 1].split(':')
                        for item in _comp:
                            if not item.strip().upper() == 'PL':
                                pl_index += 1
                            else:
                                break

                    if gt_index == -1:
                        gt_index = 0
                        _comp = comp[cell_info_index - 1].split(':')
                        for item in _comp:
                            if not item.strip().upper() == 'GT':
                                gt_index += 1
                            else:
                                break

                    for cell_index in range(0, len(cell_names)):
                        cell_info = comp[cell_info_index + cell_index].strip().split(':')
                        genotypes[cell_names[cell_index]].append(cell_info[gt_index].strip())
                        _ternary = genotype2ternary(cell_info[gt_index].strip())
                        ternary[cell_names[cell_index]].append(_ternary)

                        if mode.lower() == 'monovar':
                            if _ternary == 3:
                                genotype_probs[cell_names[cell_index]].append((0, 0, 0))
                            else:
                                _cell_info = cell_info[pl_index].strip().split(',')
                                genotype_probs[cell_names[cell_index]].append(
                                    get_genotypes_probs(int(_cell_info[0].strip()),
                                                        int(_cell_info[1].strip()),
                                                        int(_cell_info[2].strip())))

        if mode.lower() == 'sciphi' and os.path.isfile(probs_file) is True:
            probs_cell_names = []
            cell_index = []
            snv_info = {}
            with open(probs_file, 'r') as fh:
                for line in fh:
                    if line.strip().startswith('chrom'):
                        comp = line.strip().split('\t')
                        for index in range(2, len(comp)):
                            probs_cell_names.append(comp[index].strip().replace('.bam', ''))

                        if len(probs_cell_names) != len(cell_names):
                            raise ValueError('Error! The number of cells in vcf file is different from that in probs '
                                             'file.')

                        for item in cell_names:
                            cell_index.append(probs_cell_names.index(item))
                    else:
                        comp = line.strip().split('\t')
                        key = comp[0].strip() + '_' + comp[1].strip()
                        snv_info[key] = []
                        for index in range(2, len(comp)):
                            _comp = comp[index].strip().split('|')

                            # 0/0, 0/1
                            snv_info[key].append((float(_comp[1].strip()), float(_comp[0].strip())))

            for _chr_index, _pos_index in loci_info:
                for _cell_index in range(0, len(cell_names)):
                    genotype_probs[cell_names[_cell_index]].append(
                        snv_info[_chr_index + '_' + str(_pos_index)][cell_index[_cell_index]])

        return cell_names, loci_info, genotypes, ternary, genotype_probs


def print_cell_names(out_file: str, cell_names: List[str]):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, 'w') as fh:
        for index in range(0, len(cell_names)):
            fh.write(cell_names[index])
            if index < len(cell_names) - 1:
                fh.write('\n')


def print_loci_info(out_file: str, loci_info: List[Tuple[str, int]]):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, 'w') as fh:
        for index in range(0, len(loci_info)):
            fh.write(loci_info[index][0] + '\t' + str(loci_info[index][1]))
            if index < len(loci_info) - 1:
                fh.write('\n')


def print_genotypes(out_file: str, snv_num: int, cell_names: List[str], genotypes: Dict[str, List[str]]):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, 'w') as fh:
        for snv_index in range(0, snv_num):
            for cell_index in range(0, len(cell_names)):
                fh.write(genotypes[cell_names[cell_index]][snv_index])
                if cell_index < len(cell_names) - 1:
                    fh.write('\t')
            if snv_index < snv_num - 1:
                fh.write('\n')


def print_ternary(out_file: str, snv_num: int, cell_names: List[str], ternary: Dict[str, List[int]]):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, 'w') as fh:
        for snv_index in range(0, snv_num):
            for cell_index in range(0, len(cell_names)):
                fh.write(str(ternary[cell_names[cell_index]][snv_index]))
                if cell_index < len(cell_names) - 1:
                    fh.write('\t')
            if snv_index < snv_num - 1:
                fh.write('\n')


def print_genotype_probs(out_file: str, snv_num: int, cell_names: List[str],
                         genotype_probs: Dict[str, List[Union[Tuple[float, float, float], Tuple[float, float]]]]):
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    with open(out_file, 'w') as fh:
        for snv_index in range(0, snv_num):
            for cell_index in range(0, len(cell_names)):
                probs = genotype_probs[cell_names[cell_index]][snv_index]
                for genotype_index in range(0, len(probs)):
                    fh.write(str(probs[genotype_index]))
                    if genotype_index < len(probs) - 1:
                        fh.write(',')
                if cell_index < len(cell_names) - 1:
                    fh.write('\t')
            if snv_index < snv_num - 1:
                fh.write('\n')


def print_for_sifit(out_data_file: str, out_cell_names_file: str, mode: str, snv: List[Tuple[str, int]],
                    cell_names: List[str], ternary: Dict[str, List[int]]):
    if not os.path.exists(os.path.dirname(out_data_file)):
        os.makedirs(os.path.dirname(out_data_file))
    with open(out_data_file, 'w') as fh:
        for snv_index in range(0, len(snv)):
            fh.write(str(snv[snv_index][1]))
            for cell_index in range(0, len(cell_names)):
                fh.write(' ')
                code = ternary[cell_names[cell_index]][snv_index]
                if mode == 'b':
                    if code == 1 or code == 2:
                        fh.write(str(1))
                    else:
                        fh.write(str(code))
                elif mode == 't':
                    fh.write(str(code))
            if snv_index < len(snv) - 1:
                fh.write('\n')

    if not os.path.exists(os.path.dirname(out_cell_names_file)):
        os.makedirs(os.path.dirname(out_cell_names_file))
    with open(out_cell_names_file, 'w') as fh:
        for index in range(0, len(cell_names)):
            fh.write(cell_names[index])
            if index < len(cell_names) - 1:
                fh.write(' ')


def main():
    sys_args, output_files = parse_system_args()
    cell_names, loci_info, genotypes, ternary, genotype_probs = parse_vcf_file(sys_args.input,
                                                                               sys_args.mode,
                                                                               sys_args.probs
                                                                               )
    if sys_args.mode == 'monovar' or sys_args.mode == 'sciphi':
        print_cell_names(output_files['cell_names'], cell_names)
        print_loci_info(output_files['loci_info'], loci_info)
        print_genotypes(output_files['genotypes'], len(loci_info), cell_names, genotypes)
        print_ternary(output_files['ternary'], len(loci_info), cell_names, ternary)
        print_genotype_probs(output_files['genotype_probs'], len(loci_info), cell_names, genotype_probs)

    if sys_args.sifit != 'n':
        print_for_sifit(output_files['sifit_data'],
                        output_files['sifit_cell_names'],
                        sys_args.sifit,
                        loci_info,
                        cell_names, ternary
                        )


if __name__ == '__main__':
    main()
