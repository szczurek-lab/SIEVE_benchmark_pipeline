import argparse
import os
import re

import numpy as np
from typing import Union


DATASET_LABEL_PATTERN = r'/(\d+)/|\.(\d+)$'


OUTPUT_FILE_TEMPLATES: dict[str, dict[str, str]] = {
    'reads': {
        'inc_cnv': '{INC_CNV_LABEL}/{DATASET_LABEL}/read_counts.pileup',
        'exc_cnv': '{EXC_CNV_LABEL}/{DATASET_LABEL}/read_counts.pileup'
    },
    'snv_sites': {
        'inc_cnv': '{INC_CNV_LABEL}/{DATASET_LABEL}/true_snv_sites_names',
        'exc_cnv': '{EXC_CNV_LABEL}/{DATASET_LABEL}/true_snv_sites_names'
    },
    'alpha_gts': {
        'inc_cnv': '{INC_CNV_LABEL}/{DATASET_LABEL}/true_snv_genotypes_alphabetic',
        'exc_cnv': '{EXC_CNV_LABEL}/{DATASET_LABEL}/true_snv_genotypes_alphabetic'
    },
    'num_gts': {
        'inc_cnv': '{INC_CNV_LABEL}/{DATASET_LABEL}/true_snv_genotypes_numeric',
        'exc_cnv': '{EXC_CNV_LABEL}/{DATASET_LABEL}/true_snv_genotypes_numeric'
    },
    'ternary_gts': {
        'inc_cnv': '{INC_CNV_LABEL}/{DATASET_LABEL}/true_snv_ternary_matrix',
        'inc_cnv_sim': '{INC_CNV_LABEL}/{DATASET_LABEL}/true_snv_ternary_matrix.simple',
        'exc_cnv': '{EXC_CNV_LABEL}/{DATASET_LABEL}/true_snv_ternary_matrix',
        'exc_cnv_sim': '{EXC_CNV_LABEL}/{DATASET_LABEL}/true_snv_ternary_matrix.simple'
    },
    'coverage': {
        'inc_cnv': '{INC_CNV_LABEL}/{DATASET_LABEL}/true_snv_coverage',
        'exc_cnv': '{EXC_CNV_LABEL}/{DATASET_LABEL}/true_snv_coverage'
    },
    'evo_events': {
        'inc_cnv': '{INC_CNV_LABEL}/{DATASET_LABEL}/evolutionary_events_summary',
        'exc_cnv': '{EXC_CNV_LABEL}/{DATASET_LABEL}/evolutionary_events_summary'
    }
}


class CellReads:
    __gt_sep__: str = '|'
    __nuc__: str = ['A', 'C', 'G', 'T']
    __qual__: str = '~'

    def __init__(
        self,
        gt: str,
        ref_nuc: str,
        depth: int,
        reads: str,
        update_gt_sep: bool = False,
        has_cnv: bool = False
    ):
        self.ori_gt = gt
        self.gt = gt
        self.__update_gt_flag = True

        if update_gt_sep:
            CellReads.__gt_sep__ = re.findall(r'[^\\.\w]', gt)[0]

        self.alleles = [CellReads.__nuc__.index(i) if i in CellReads.__nuc__ else -1 for i in re.split(r'[/|]', gt)]
        if -1 in self.alleles:
            self.alleles.sort(reverse=True)
        else:
            self.alleles.sort()

        self.ref_nuc = CellReads.__nuc__.index(ref_nuc)
        self.depth = depth
        self.allele_counts = [0 for _ in range(len(CellReads.__nuc__))]
        self.reads_alpha_2_num(reads)

        self.has_cnv = has_cnv

        self.num_gt = None
        self.__convert_gt_2_num_gt()
        self.ternary_gt = None
        self.__convert_gt_2_ternary_gt()

    def __convert_gt_2_num_gt(self):
        if all(np.unique(self.alleles) == -1):
            self.num_gt = '.'
            return

        alt_label = None
        alt_code = 0
        tmp = []
        for i in self.alleles:
            if i == -1:
                tmp.append(100)
            elif i == self.ref_nuc:
                tmp.append(0)
            else:
                if alt_label is None or i != alt_label:
                    alt_label = i
                    alt_code += 1
                    tmp.append(alt_code)
                else:
                    tmp.append(alt_code)
        tmp.sort()
        self.num_gt = CellReads.__gt_sep__.join([str(i) if i != 100 else '.' for i in tmp])

    def __convert_gt_2_ternary_gt(self):
        # -3: missing genotype (.)
        # -2: 1/.
        # -1: 0/.
        # 0: homozygous reference (0/0)
        # 1: heterozygous mutation (0/1)
        # 2: homozygous mutation with homozygous alternative nucletides (1/1)
        # 3: homozygous mutation with heterozygous alternative nucleotides (1/1')
        if all(np.unique(self.alleles) == -1):
            self.ternary_gt = -3
        elif -1 in self.alleles:
            if self.ref_nuc in self.alleles:
                self.ternary_gt = -1
            else:
                self.ternary_gt = -2
        elif len(np.unique(self.alleles)) == 1:
            if self.ref_nuc in self.alleles:
                self.ternary_gt = 0
            else:
                self.ternary_gt = 2
        else:
            if self.ref_nuc in self.alleles:
                self.ternary_gt = 1
            else:
                self.ternary_gt = 2 #3

    def reads_alpha_2_num(self, reads: str):
        if reads.strip():
            for i in reads.strip():
                if i == '*':
                    continue
                elif i == '.':
                    self.allele_counts[self.ref_nuc] += 1
                else:
                    self.allele_counts[CellReads.__nuc__.index(i)] += 1

    @staticmethod
    def get_allele_alpha(num: int):
        if num == -1:
            return '.'
        else:
            return CellReads.__nuc__[num]

    def update_gt(self):
        if not self.__update_gt_flag:
            return
        if all(np.unique(self.alleles) == -1):
            self.gt = '.'
        else:
            self.gt = CellReads.__gt_sep__.join([CellReads.get_allele_alpha(i) for i in self.alleles])
        self.__convert_gt_2_num_gt()
        self.__convert_gt_2_ternary_gt()
        self.__update_gt_flag = False

    def apply_cnv(self, rng, cnv_num: int):
        if cnv_num == 0:
            self.alleles = [-1 for _ in range(len(self.alleles))]
            self.depth = 0
            self.allele_counts = [0 for _ in range(len(CellReads.__nuc__))]
        else:
            while True:
                allele_idx = rng.integers(0, len(self.alleles), 1)[0]
                if self.allele_counts[self.alleles[allele_idx]] > 0:
                    break
            if self.allele_counts[self.alleles[allele_idx]] == self.depth and len(np.unique(self.alleles)) == 1:
                ori_reads = round(self.depth / len(self.alleles))
                left_reads = self.depth - ori_reads
            else:
                ori_reads = self.allele_counts[self.alleles[allele_idx]]
                left_reads = 0
            cnv_reads = ori_reads * (cnv_num - 1)
            self.allele_counts[self.alleles[allele_idx]] = cnv_reads + left_reads
            self.depth = sum(self.allele_counts)

            if cnv_reads == 0:
                self.alleles[allele_idx] = -1

            if -1 in self.alleles:
                self.alleles.sort(reverse=True)
            else:
                self.alleles.sort()

        if cnv_num == 0 or cnv_num == 1:
            self.__update_gt_flag = True

    def is_snv(self):
        for i in self.alleles:
            if i != -1 and i != self.ref_nuc:
                return True
        return False

    def get_gt(self):
        self.update_gt()
        return self.gt

    def get_num_gt(self):
        self.update_gt()
        return self.num_gt

    def get_ternary_gt(self):
        self.update_gt()
        return self.ternary_gt

    def get_simplified_ternary_gt(self):
        self.update_gt()
        if self.ternary_gt == -3:
            return 3
        elif self.ternary_gt == -2 or self.ternary_gt == 2 or self.ternary_gt == 3:
            return 2
        elif self.ternary_gt == -1 or self.ternary_gt == 0:
            return 0
        else:
            return self.ternary_gt

    def __str__(self):
        if self.depth < 0:
            return None
        elif self.depth == 0:
            return '0\t*\t*'
        else:
            return '{}\t{}\t{}'.format(
                self.depth,
                ''.join([('.' if i == self.ref_nuc else CellReads.__nuc__[i]) * self.allele_counts[i] for i in
                         range(len(CellReads.__nuc__))]),
                CellReads.__qual__ * self.depth
            )


class GenomicSite:
    def __init__(self, chr: str, pos: int, ref_nuc: str, is_snv: bool, has_cnv: bool = False, cnv_type: int = -1, cnv: int = 2):
        """
        This function initializes a new `Variant` object with the given parameters

        :param chr: chromosome
        :type chr: str
        :param pos: the position of the site in the genome
        :type pos: int
        :param ref_nuc: the reference nucleotide at this position
        :type ref_nuc: str
        :param is_snv: True if the mutation is a single nucleotide variant (SNV), False if it's not
        :type is_snv: bool
        :param has_cnv: whether the site has a CNV, defaults to False
        :type has_cnv: bool (optional)
        :param cnv_type: -1 = no CNV, 0 = less than 1/3 of cells contains a CNV, 1 = more than 2/3 of cells contains a CNV, defaults to -1
        :type cnv_type: int
        :param cnv: the number of CNV, defaults to 2
        :type cnv: int
        """
        self.chr = chr
        self.pos = pos
        self.ref_nuc = ref_nuc
        self.is_snv = is_snv
        self.has_cnv = has_cnv
        self.cnv_type = cnv_type
        self.cnv = cnv
        if self.has_cnv and self.cnv_type == -1 or not self.has_cnv and self.cnv_type > -1:
            raise ValueError('CNV type and CNV status do not match')
        self.cells = []
        self.num_of_cells_with_cnv = 0

    def update_cnv_snv_properties(self):
        if self.num_of_cells_with_cnv == 0:
            self.has_cnv = False
            self.cnv_type = -1
            self.cnv = 2

        if not self.has_cnv or not self.is_snv:
            return
        for i in self.cells:
            if i.is_snv():
                self.is_snv = True
                return
        self.is_snv = False

    def add_cell(self, cell: CellReads):
        self.cells.append(cell)
        if cell.has_cnv:
            self.num_of_cells_with_cnv += 1

    def __str__(self):
        return '{}\t{}\t{}\t{}'.format(self.chr, self.pos, self.ref_nuc, '\t'.join([str(i) for i in self.cells]))


def input_files_sanity_check(reads_file: str, gt_file: str, cnv_file: str, evo_file: str):
    if not os.path.isfile(reads_file):
        raise FileNotFoundError('{} does not exist'.format(reads_file))
    if not os.path.isfile(gt_file):
        raise FileNotFoundError('{} does not exist'.format(gt_file))
    if not os.path.isfile(cnv_file):
        raise FileNotFoundError('{} does not exist'.format(cnv_file))
    if not os.path.isfile(evo_file):
        raise FileNotFoundError('{} does not exist'.format(evo_file))

    dataset_labels: set[str] = set()
    for i in [reads_file, gt_file, cnv_file, evo_file]:
        dataset_labels.add(*{j for j in re.search(DATASET_LABEL_PATTERN, i).groups() if j is not None})

    if len(dataset_labels) != 1:
        raise ValueError('Input files are not from the same dataset')

    return dataset_labels.pop()


def prepare_output_files(prefix: str, inc_cnv_lab: str, exc_cnv_lab: str, dataset_label: str) -> dict[str, dict[str, str]]:
    ret: dict[str, dict[str, str]] = {}
    for key_1, value_1 in OUTPUT_FILE_TEMPLATES.items():
        ret[key_1] = {}
        for key_2, value_2 in value_1.items():
            ret[key_1][key_2] = ''.join([prefix, value_2.format(
                INC_CNV_LABEL=inc_cnv_lab,
                EXC_CNV_LABEL=exc_cnv_lab,
                DATASET_LABEL=dataset_label
            )])
    return ret


def parse_system_args() -> tuple[argparse.Namespace, dict[str, dict[str, str]]]:
    sys_args_parser = argparse.ArgumentParser(description='Add CNVs to the simulated data.')
    sys_args_parser.add_argument(
        '--reads',
        required=True,
        type=str,
        help='read counts file in mpileup format'
    )
    sys_args_parser.add_argument(
        '--snvs',
        required=True,
        type=str,
        help='file of true SNV sites'
    )
    sys_args_parser.add_argument(
        '--gts',
        required=True,
        type=str,
        help='file of true genotypes'
    )
    sys_args_parser.add_argument(
        '--evo',
        required=True,
        type=str,
        help='file of true evolutionary events'
    )
    sys_args_parser.add_argument(
        '--prob',
        type=float,
        default=1 / 3,
        help='probability of a site containing CNVs'
    )
    sys_args_parser.add_argument(
        '--cnvs',
        nargs=2,
        type=int,
        default=[0, 5],
        help='the number of CNVs to add to the data (min and max)'
    )
    sys_args_parser.add_argument(
        '--seed',
        type=int,
        required=True,
        help='seed for the random number generator'
    )
    sys_args_parser.add_argument(
        '--prefix',
        required=True,
        type=str,
        help='prefix for the output files'
    )
    sys_args_parser.add_argument(
        '--inccnvlabel',
        required=False,
        type=str,
        default='inc_cnv',
        help='label for results including CNVs'
    )
    sys_args_parser.add_argument(
        '--exccnvlabel',
        required=False,
        type=str,
        default='exc_cnv',
        help='label for results excluding CNVs'
    )
    sys_args_parser.add_argument(
        '--keep',
        required=False,
        type=int,
        choices=[1, 2, 3],
        default=3,
        help='keep output types: 1 - only including CNVs, 2 - only excluding CNVs, 3 - both'
    )

    sys_args = sys_args_parser.parse_args()

    if sys_args.prob <= 0 or sys_args.prob > 1:
        sys_args_parser.error('--prob must be between 0 (excluded) and 1 (included)')

    if sys_args.cnvs[0] < 0 or sys_args.cnvs[1] < sys_args.cnvs[0]:
        sys_args_parser.error('--cnvs must be positive and in ascending order')

    if sys_args.cnvs[0] == sys_args.cnvs[1] == 2:
        sys_args_parser.error('The genome is assumed to be diploid, and thus no CNVs will be added.')

    return sys_args, prepare_output_files(
        sys_args.prefix,
        sys_args.inccnvlabel,
        sys_args.exccnvlabel,
        input_files_sanity_check(
            sys_args.reads,
            sys_args.snvs,
            sys_args.gts,
            sys_args.evo
        )
    )


def get_snv_sites(file: str) -> Union[list[tuple[str, int]], None]:
    if os.path.isfile(file):
        ret: list[tuple[str, int]] = []
        with open(file, 'r') as fh:
            for line in fh:
                if line.strip():
                    chrom, pos = line.strip().split()
                    ret.append((chrom, int(pos)))
        return ret
    return None


def get_gt_snv_sites(file: str) -> Union[list[list[str]], None]:
    if os.path.isfile(file):
        ret = []
        with open(file, 'r') as fh:
            for line in fh:
                if line.strip():
                    ret.append(line.strip().split())
        return ret
    return None


def get_cnv_prob(cnv_type: int, rng) -> float:
    if cnv_type == 0:
        return rng.uniform(0, 1 / 3, 1)[0]
    elif cnv_type == 1:
        return rng.uniform(2 / 3, 1, 1)[0]
    else:
        return 0


def perform_cnv(
    reads_file: str,
    snv_sites: list[tuple[str, int]],
    gt_snv_sites: list[list[str]],
    prob: float,
    cnv_range: tuple[int, int],
    seed: int
) -> tuple[Union[list[GenomicSite], None], list[int, int]]:
    cnv_num_range = [2, 2]
    if os.path.isfile(reads_file):
        ret = []
        rng = np.random.default_rng(seed)
        cell_num = len(gt_snv_sites[0])
        infer_gt_sep = True
        with open(reads_file, 'r') as fh:
            for line in fh:
                if line.strip():
                    comb = line.strip().split()
                    snv_idx = snv_sites.index((comb[0], int(comb[1]))) if (comb[0], int(comb[1])) in snv_sites else -1

                    # Whether this site will be affected by CNVs?
                    has_cnv = rng.random(1)[0] <= prob
                    cnv_type = (0 if rng.random(1)[0] <= 0.5 else 1) if has_cnv else -1
                    cnv_prob = get_cnv_prob(cnv_type, rng)
                    if has_cnv:
                        while True:
                            cnv_num = rng.integers(*cnv_range, size=1)[0]
                            if cnv_num != 2:
                                break
                    else:
                        cnv_num = 2

                    if cnv_num < cnv_num_range[0]:
                        cnv_num_range[0] = cnv_num
                    if cnv_num > cnv_num_range[1]:
                        cnv_num_range[1] = cnv_num

                    snv_site = GenomicSite(comb[0], int(comb[1]), comb[2], snv_idx > -1, has_cnv, cnv_type, cnv_num)

                    if len(comb) / 3 - 1 != cell_num:
                        raise ValueError(f'Expected {cell_num} cells, but found {len(comb) / 3 - 1}')

                    for i in range(cell_num):
                        update_gt_sep = infer_gt_sep and snv_idx > -1
                        cell_has_cnv = rng.random(1)[0] <= cnv_prob if has_cnv else False
                        cell_reads = CellReads(
                            gt=gt_snv_sites[snv_idx][i] if snv_idx > -1 else '|'.join(np.repeat(comb[2], 2)),
                            ref_nuc=comb[2],
                            depth=int(comb[3 + 3 * i]),
                            reads=comb[4 + 3 * i],
                            update_gt_sep=update_gt_sep,
                            has_cnv=cell_has_cnv
                        )
                        if update_gt_sep:
                            infer_gt_sep = False

                        if cell_has_cnv:
                            cell_reads.apply_cnv(rng, cnv_num)

                        snv_site.add_cell(cell_reads)

                    snv_site.update_cnv_snv_properties()
                    ret.append(snv_site)
        return ret, cnv_num_range
    return None, cnv_num_range


def get_evo_events(
    file: str,
    snv_sites: list[tuple[str, int]]
) -> Union[tuple[dict[tuple[str, int], list[str]], list[str]], tuple[None, None]]:
    if os.path.isfile(file):
        ret: dict[tuple[str, int], list[str]] = {}
        evo_events = []
        evo_events_comments = []
        with open(file, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line:
                    if line.startswith('#'):
                        evo_events_comments.append(line)
                    else:
                        evo_events.append(line.split())
        if len(snv_sites) != len(evo_events):
            raise ValueError(f'Expected {len(snv_sites)} events, but found {len(evo_events)}')
        for i, j in zip(snv_sites, evo_events):
            ret[i] = j
        return ret, evo_events_comments
    return None, None


def write_to_reads(site_info: list[GenomicSite], output_file: dict[str, str], keep: int = 3) -> None:
    for inc_cnv, (_, file) in zip([True, False], output_file.items()):
        if keep == 3 or keep == 2 and not inc_cnv or keep == 1 and inc_cnv:
            os.makedirs(os.path.dirname(file), exist_ok=True)
            with open(file, 'w') as fh:
                for i in site_info:
                    if inc_cnv or not i.has_cnv:
                        fh.write(f'{str(i)}\n')


def write_to_snv_sites(site_info: list[GenomicSite], output_file: dict[str, str], keep: int = 3) -> None:
    for inc_cnv, (_, file) in zip([True, False], output_file.items()):
        if keep == 3 or keep == 2 and not inc_cnv or keep == 1 and inc_cnv:
            os.makedirs(os.path.dirname(file), exist_ok=True)
            with open(file, 'w') as fh:
                for i in site_info:
                    if i.is_snv and (inc_cnv or not i.has_cnv):
                        fh.write(f'{i.chr} {i.pos}\n')


def write_to_alpha_gts(site_info: list[GenomicSite], output_file: dict[str, str], keep: int = 3) -> None:
    for inc_cnv, (_, file) in zip([True, False], output_file.items()):
        if keep == 3 or keep == 2 and not inc_cnv or keep == 1 and inc_cnv:
            os.makedirs(os.path.dirname(file), exist_ok=True)
            with open(file, 'w') as fh:
                for i in site_info:
                    if i.is_snv and (inc_cnv or not i.has_cnv):
                        for j in range(len(i.cells)):
                            fh.write(i.cells[j].get_gt())
                            if j != len(i.cells) - 1:
                                fh.write(' ')
                        fh.write('\n')


def write_to_num_gts(site_info: list[GenomicSite], output_file: dict[str, str], keep: int = 3) -> None:
    for inc_cnv, (_, file) in zip([True, False], output_file.items()):
        if keep == 3 or keep == 2 and not inc_cnv or keep == 1 and inc_cnv:
            os.makedirs(os.path.dirname(file), exist_ok=True)
            with open(file, 'w') as fh:
                for i in site_info:
                    if i.is_snv and (inc_cnv or not i.has_cnv):
                        for j in range(len(i.cells)):
                            fh.write(i.cells[j].get_num_gt())
                            if j != len(i.cells) - 1:
                                fh.write(' ')
                        fh.write('\n')


def write_to_ternary_gts(site_info: list[GenomicSite], output_file: dict[str, str], keep: int = 3) -> None:
    for inc_cnv, simplified, (_, file) in zip([True, True, False, False], [False, True, False, True], output_file.items()):
        if keep == 3 or keep == 2 and not inc_cnv or keep == 1 and inc_cnv:
            os.makedirs(os.path.dirname(file), exist_ok=True)
            with open(file, 'w') as fh:
                for i in site_info:
                    if i.is_snv and (inc_cnv or not i.has_cnv):
                        for j in range(len(i.cells)):
                            fh.write(f'{str(i.cells[j].get_simplified_ternary_gt()) if simplified else i.cells[j].get_ternary_gt()} ')
                            if j != len(i.cells) - 1:
                                fh.write(' ')
                        fh.write('\n')


def write_to_coverage(site_info: list[GenomicSite], output_file: dict[str, str], keep: int = 3) -> None:
    for inc_cnv, (_, file) in zip([True, False], output_file.items()):
        if keep == 3 or keep == 2 and not inc_cnv or keep == 1 and inc_cnv:
            os.makedirs(os.path.dirname(file), exist_ok=True)
            with open(file, 'w') as fh:
                for i in site_info:
                    if i.is_snv and (inc_cnv or not i.has_cnv):
                        for j in range(len(i.cells)):
                            fh.write(str(i.cells[j].depth))
                            if j != len(i.cells) - 1:
                                fh.write('\t')
                        fh.write('\n')


def write_to_evo_events(
    site_info: list[GenomicSite],
    output_file: dict[str, str],
    evo_events: dict[tuple[str, int], list[str]],
    evo_event_comments: Union[list[str], None],
    cnv_num_range: list[int, int],
    keep: int = 3
) -> None:
    for inc_cnv, (_, file) in zip([True, False], output_file.items()):
        if keep == 3 or keep == 2 and not inc_cnv or keep == 1 and inc_cnv:
            os.makedirs(os.path.dirname(file), exist_ok=True)
            with open(file, 'w') as fh:
                for i in evo_event_comments:
                    if i.startswith('#Cell-wise'):
                        fh.write(f'{i}\n')
                    else:
                        fh.write(
                            f'#Cell-wise number of CNVs = [{cnv_num_range[0] if inc_cnv else 2}, {cnv_num_range[1] if inc_cnv else 2}]\n'
                        )
                        fh.write(f'{i}, CNVs\n')
                for i in site_info:
                    if i.is_snv and (inc_cnv or not i.has_cnv):
                        __evo_events = evo_events[(i.chr, i.pos)]
                        fh.write('{}\n'.format('\t'.join([','.join([__evo_events[j], str(i.cnv) if i.cells[j].has_cnv else str(2)]) for j in range(len(__evo_events))])))


def write_to_output_files(
    site_info: list[GenomicSite],
    output_files: dict[str, dict[str, str]],
    evo_events: dict[tuple[str, int], list[str]],
    evo_event_comments: Union[list[str], None],
    cnv_num_range: list[int, int],
    keep: int = 3
) -> None:
    write_to_reads(site_info, output_files['reads'], keep)
    write_to_snv_sites(site_info, output_files['snv_sites'], keep)
    write_to_alpha_gts(site_info, output_files['alpha_gts'], keep)
    write_to_num_gts(site_info, output_files['num_gts'], keep)
    write_to_ternary_gts(site_info, output_files['ternary_gts'], keep)
    write_to_coverage(site_info, output_files['coverage'], keep)
    write_to_evo_events(site_info, output_files['evo_events'], evo_events, evo_event_comments, cnv_num_range, keep)


def main():
    sys_args, output_files = parse_system_args()
    snv_sites = get_snv_sites(sys_args.snvs)
    gt_snv_sites = get_gt_snv_sites(sys_args.gts)
    sites_with_cnvs, cnv_num_range = perform_cnv(
        sys_args.reads,
        snv_sites,
        gt_snv_sites,
        sys_args.prob,
        sys_args.cnvs,
        sys_args.seed
    )
    evo_events, evo_event_comments = get_evo_events(sys_args.evo, snv_sites)
    write_to_output_files(sites_with_cnvs, output_files, evo_events, evo_event_comments, cnv_num_range, sys_args.keep)


if __name__ == '__main__':
    main()
