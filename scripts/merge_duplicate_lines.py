import argparse
import os
import re
from typing import List, Tuple, Union

CHROMOSOME_PATTERN = re.compile('([a-zA-Z]*)([1-9]|1\\d|2[012]|[XYxy])$')


class LocusInfo:
    def __init__(
            self,
            chr: str,
            pos: int,
            ref_nuc: str = None,
            alt_nuc: str = None
    ):
        matcher = re.match(CHROMOSOME_PATTERN, chr)
        if matcher is None:
            raise ValueError('Error! Unsupported chromosome label.')
        self.chr_full_label = chr
        self.chr_label = matcher.group(2)
        self.pos = pos
        self.ref_nuc = ref_nuc
        self.alt_nuc = alt_nuc

    def __eq__(self, o: Union[None, 'LocusInfo']) -> bool:
        return False if o is None else self.chr_label == o.chr_label and self.pos == o.pos and self.ref_nuc == o.ref_nuc

    def __ne__(self, o: Union[None, 'LocusInfo']) -> bool:
        return True if o is None else self.chr_label != o.chr_label or self.pos != o.pos or self.ref_nuc != o.ref_nuc

    def to_list(self) -> List[str]:
        return [self.chr_full_label, str(self.pos), self.ref_nuc, self.alt_nuc] \
            if self.ref_nuc is not None and self.alt_nuc is not None else \
            [self.chr_full_label, str(self.pos)]


def parse_system_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(description='Merge lines with the same chromosome label and position.')
    sys_args_parser.add_argument('--loci', help='file containing loci information',
                                 metavar='LOCI', type=str, required=True)
    sys_args_parser.add_argument('--ternary', help='file containing ternary genotypes',
                                 metavar='TERNARY', type=str, required=True)
    sys_args_parser.add_argument('--prefix', help='prefix to output file',
                                 metavar='OUTPUT', type=str, required=True)
    sys_args = sys_args_parser.parse_args()
    if not os.path.isfile(sys_args.loci):
        raise ValueError('ERROR! File does not exist: ' + sys_args.loci)
    if not os.path.isfile(sys_args.ternary):
        raise ValueError('ERROR! File does not exist: ' + sys_args.ternary)
    return sys_args


def parse_loci_info(loci_file_name: str) -> List[LocusInfo]:
    results: List[LocusInfo] = []
    with open(loci_file_name, 'r') as fh:
        for line in fh:
            if line is not None:
                comp: List[str] = re.split('\\s+', line.strip())
                if len(comp) == 2:
                    results.append(
                        LocusInfo(
                            comp[0].strip(),
                            int(comp[1])
                        )
                    )
                else:
                    results.append(
                        LocusInfo(
                            comp[0].strip(),
                            int(comp[1]),
                            comp[2].strip(),
                            comp[3].strip()
                        )
                    )
    return results


def merge_genotypes(genotypes1: List[int], genotypes2: List[int], destination: List[int]) -> List[int]:
    for index in range(len(destination)):
        if genotypes1[index] == 2 or genotypes2[index] == 2:
            destination[index] = 2
        else:
            destination[index] = genotypes1[index] + genotypes2[index]


def parse_cell_genotypes(genotype_file_name: str, loci_info: List[LocusInfo]) -> List[Tuple[LocusInfo, List[int]]]:
    results: List[Tuple[LocusInfo, List[int]]] = []

    with open(genotype_file_name, 'r') as fh:
        prev_locus_info = None
        line_index: int = 0

        for line in fh:
            if line is not None:
                comp: List[str] = re.split('\\s+', line.strip())

                if loci_info[line_index] == prev_locus_info:
                    merge_genotypes(results[-1][1], [int(i) for i in comp], results[-1][1])
                else:
                    prev_locus_info = loci_info[line_index]
                    results.append(
                        (
                            loci_info[line_index],
                            [int(i) for i in comp]
                        )
                    )

                line_index += 1

    return results


def prepare_output_files(prefix: str, file_name: str) -> str:
    comp: List[str] = file_name.strip().split('.')
    return '.'.join([prefix, comp[1]])


def write2output(
        merged_genotypes: List[Tuple[LocusInfo, List[int]]],
        loci_info_output_file_name: str,
        ternary_output_file_name: str
) -> None:
    fh_loci_info = open(loci_info_output_file_name, 'w')
    fh_ternary = open(ternary_output_file_name, 'w')

    for index in range(len(merged_genotypes)):
        fh_loci_info.write('\t'.join(merged_genotypes[index][0].to_list()))
        fh_ternary.write('\t'.join([str(i) for i in merged_genotypes[index][1]]))

        if index < len(merged_genotypes) - 1:
            fh_loci_info.write('\n')
            fh_ternary.write('\n')

    fh_loci_info.close()
    fh_ternary.close()


def main() -> None:
    sys_args: argparse.Namespace = parse_system_args()
    loci_info: List[LocusInfo] = parse_loci_info(sys_args.loci)
    merged_genotypes: List[Tuple[LocusInfo, List[int]]] = parse_cell_genotypes(sys_args.ternary, loci_info)
    write2output(
        merged_genotypes,
        prepare_output_files(sys_args.prefix, sys_args.loci),
        prepare_output_files(sys_args.prefix, sys_args.ternary)
    )


if __name__ == '__main__':
    main()
