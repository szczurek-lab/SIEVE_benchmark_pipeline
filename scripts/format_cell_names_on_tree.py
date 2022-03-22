import argparse
import os
import re
from typing import Dict, List, Tuple


CELL_NAME = 'tumcell'
CELL_NAME_PATTERN = '[(,](' + CELL_NAME + '.*?):'


def parse_sys_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(description='Format cell names on a tree according to a reference file.')
    sys_args_parser.add_argument('--tree', help='tree file to be formatted',
                                 metavar='TREE', type=str, required=True)
    sys_args_parser.add_argument('--ref', help='reference file containing standard cell names', 
                                 metavar='REFERENCE CELL NAMES', type=str, required=True)
    sys_args_parser.add_argument('--out', help='output tree file name',
                                 metavar='OUTPUT', type=str, required=True)
    
    sys_args = sys_args_parser.parse_args()
    if not os.path.isfile(sys_args.tree):
        raise ValueError("ERROR! Could not find " + sys_args.tree)
    if not os.path.isfile(sys_args.ref):
        raise ValueError("ERROR! Could not find " + sys_args.ref)
    if not os.path.isdir(os.path.dirname(sys_args.out)):
        os.makedirs(os.path.dirname(sys_args.out), exist_ok=True)

    return sys_args


def build_cell_map(tree_file_name: str, ref_file_name: str) -> Tuple[Dict[str,str], str] :
    tree_content: str = ''
    ref_cell_names: List[str] = []

    with open(tree_file_name, 'r') as fh:
        for line in fh:
            if line is not None:
                tree_content = line.strip()
                break
    
    with open(ref_file_name, 'r') as fh:
        for line in fh:
            if line is not None:
                for item in line.strip().split('\s+'):
                    ref_cell_names.append(item)

    pattern = re.compile(CELL_NAME_PATTERN)
    tree_matched_cells = pattern.findall(tree_content)
    
    if len(tree_matched_cells) != len(ref_cell_names):
        raise ValueError('ERROR! Number of cells in two files don\'t match!')

    copy_ref_cell_names = ref_cell_names.copy()
    copy_tree_matched_cells = tree_matched_cells.copy()

    tree_cell_2_ref = {}

    for ref in ref_cell_names:
        ref_pattern = re.compile(ref)

        for tree_cell in tree_matched_cells:
            if re.match(ref_pattern, tree_cell) is not None:
                tree_cell_2_ref[tree_cell] = ref
                copy_ref_cell_names.remove(ref)
                copy_tree_matched_cells.remove(tree_cell)
    
    if len(copy_tree_matched_cells) != 0 or len(copy_ref_cell_names) != 0:
        raise ValueError('ERROR! Cell names do not match in both files!')

    return tree_cell_2_ref, tree_content


def replace_cell_names(cell_map: Dict[str, str], tree: str) -> str:
    for tree_cell, ref in cell_map.items():
        tree = tree.replace(tree_cell, ref)
    return tree


def print_formatted_tree(tree, output) -> None:
    with open(output, 'w') as fh:
        fh.write(tree)



def main() -> None:
    sys_args: argparse.Namespace = parse_sys_args()
    cell_map, tree = build_cell_map(sys_args.tree, sys_args.ref)
    tree = replace_cell_names(cell_map, tree)
    print_formatted_tree(tree, sys_args.out)


if __name__ == '__main__':
    main()
