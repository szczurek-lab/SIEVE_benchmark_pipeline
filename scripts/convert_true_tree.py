import argparse
import os.path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.BaseTree import Tree


def parse_system_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(
        description='Convert branch length of trees from generations to mutations.'
    )
    sys_args_parser.add_argument(
        '-i',
        required=True,
        type=str,
        help='input tree file'
    )
    sys_args_parser.add_argument(
        '-m',
        required=True,
        type=float,
        help='mutation rate')
    sys_args_parser.add_argument(
        '--accu',
        type=str,
        default='25',
        help='accuracy of output branch length')
    sys_args_parser.add_argument(
        '-o',
        required=True,
        type=str,
        help='output tree file'
    )
    args = sys_args_parser.parse_args()
    if not os.path.isfile(args.i):
        raise(ValueError, 'Input file does not exist.')
    os.makedirs(os.path.dirname(args.o), exist_ok=True)
    if os.path.isfile(args.o):
        print('Output file exists; will be overwritten.')
    return args


def convert_branch_length(node: Clade, mu: float) -> None:
    node.branch_length = node.branch_length * mu
    for child in node.clades:
        convert_branch_length(child, mu)


def main() -> None:
    args: argparse.Namespace = parse_system_args()
    tree: Tree = Phylo.read(args.i, 'newick')
    convert_branch_length(tree.root, args.m)
    Phylo.write(tree, args.o, 'newick', format_branch_length='%.' + args.accu + 'f')


if __name__ == '__main__':
    main()
