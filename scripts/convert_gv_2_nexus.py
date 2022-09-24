import argparse
import os
import sys
from typing import Optional

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from scripts import tree_str


def parse_system_args() -> tuple[argparse.Namespace, dict[str, str]]:
    sys_args_parser = argparse.ArgumentParser(description='Convert gv trees to nexus format.')
    sys_args_parser.add_argument(
        'input',
        nargs='+',
        help='one or more files of gv format'
    )
    sys_args_parser.add_argument(
        '--leafNames',
        required=True,
        type=str,
        help='a list of leaf names used to generate gv files'
    )
    sys_args_parser.add_argument(
        '--output',
        type=str,
        help='output file or directory. If only one input tree is provided, an output file or directory can be specified. If more than one input trees are provided, this option should not be specified as the out put files will be placed with each input tree.'
    )
    sys_args_parser.add_argument(
        '--len',
        type=float,
        default=0.1,
        help='branch length of each node on the output tree.'
    )
    sys_args_parser.add_argument(
        '--startIndex',
        type=int,
        default=-1,
        help='the starting index to leaf nodes in the gv file (must consecutively be leaf nodes from this index). If not provided, the programme will try to parse it from the input file.'
    )
    sys_args = sys_args_parser.parse_args()
    if len(sys_args.input) > 1 and sys_args.output is not None:
        raise ValueError('Error! --out should not be defined as more than one input trees are provided.')
    try:
        return sys_args, prepare_output_tree_files(sys_args.input, sys_args.output)
    except ValueError:
        exit(1)


def get_leaf_names(leaf_names_file: str) -> list[str]:
    leaf_names = []
    with open(leaf_names_file, 'r') as fh:
        for line in fh:
            comp = line.strip().split('\t')
            if len(comp) == 2:
                if comp[1].strip().upper() == 'CT':
                    leaf_names.append(comp[0].strip().split('/')[-1].strip().split('.')[0].strip())
            else:
                raise IOError('Error! The file of leaf names (' + leaf_names_file + ') is not the one used to '
                                                                                    'generate gv file.')
    return leaf_names


def get_root_path(file_name: str) -> str:
    if os.path.isfile(file_name):
        return os.path.dirname(file_name)
    elif os.path.isdir(file_name):
        return file_name


def get_base_name(file_name: str) -> str:
    if os.path.isdir(file_name):
        return os.path.basename(file_name)
    elif os.path.isfile(file_name):
        return '.'.join(os.path.basename(file_name).split('.')[0:-1])


def prepare_output_tree_files(input_tree_files: list[str], output: str = None) -> dict[str, str]:
    input2output = {}
    for input_tree_file in input_tree_files:
        if os.path.isfile(input_tree_file):
            if output is None:
                input2output[input_tree_file] = '/'.join(
                    [
                        get_root_path(input_tree_file),
                        get_base_name(input_tree_file) + '.nexus'
                    ]
                )
            else:
                if os.path.exists(output):
                    if os.path.isfile(output):
                        input2output[input_tree_file] = output
                    elif os.path.isdir(output):
                        input2output[input_tree_file] = '/'.join(
                            [output.rstrip('/'), get_base_name(input_tree_file) + '.nexus'])
                else:
                    input2output[input_tree_file] = output
        else:
            raise ValueError('Error! Input should be a file(s), not a directory: ' + input_tree_file)
    return input2output


def get_index_and_label(line: str) -> tuple[int, int]:
    start = end = -1
    for i in range(0, len(line)):
        if line[i] == '[':
            start = i + 1
        elif line[i] == ']':
            end = i
    if start == -1 or end == -1:
        raise IOError('Error! At least one of \'[\' and \']\' is not in the following string: ' + line)
    elif start == end:
        raise IOError('Error! Nothing between \'[\' and \']\' is detected in the following string: ' + line)
    elif start > end:
        raise IOError('Error! Problematic format detected in the following string: ' + line)
    str_in_brackets_list = line[start:end].split(',')
    for item in str_in_brackets_list:
        str_list = item.strip().split('=')
        if str_list[0].strip() == 'label':
            return int(line[:start - 1]), int(str_list[1].strip().strip('"'))


def get_parent_child(line: str) -> tuple[int, int]:
    str_parent_child = line.split('->')
    if len(str_parent_child) != 2:
        raise IOError(
            'Error! Exact 2 nodes are needed to form a parent-child relationship in the following string: ' + line)
    else:
        return int(str_parent_child[0].strip()), int(str_parent_child[1].strip())


def get_extra_info(line: str) -> Optional[tuple[str, str]]:
    start = end = -1
    for i in range(0, len(line)):
        if line[i] == '[':
            start = i + 1
        elif line[i] == ']':
            end = i
    if start == -1 or end == -1 or start >= end:
        return None
    str_in_brackets_list = line[start:end].split(',')
    for item in str_in_brackets_list:
        str_list = item.strip().split('=')
        if str_list[0].strip() != 'label':
            return str_list[0].strip(), str_list[1].strip()


def read_trees(input_tree_files: list[str],
               branch_length: float,
               start_index: int,
               input_leaf_names: list[str]
               ) -> dict[str, tree_str.Tree]:
    _double_check = False
    if start_index == -1:
        start_index = len(input_leaf_names) - 1
        _double_check = True
    end_index = start_index + len(input_leaf_names) - 1
    trees = {}
    for input_tree_file in input_tree_files:
        tree = tree_str.Tree()
        index2label = {}
        label2node = {}
        with open(input_tree_file, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line.endswith(';'):
                    line = line[:len(line) - 1].strip()
                    if '[' in line and ']' in line:
                        # We are in the definition of nodes
                        try:
                            index, label = get_index_and_label(line)
                            index2label[index] = label
                        except IOError:
                            exit(1)
                        if start_index <= index <= end_index:
                            # A leaf node
                            if _double_check:
                                if get_extra_info(line) is None:
                                    raise ValueError('Error! Without --startIndex specified, the programme cannot be '
                                                     'sure that this node (' + line + ') is a leaf node. Try to '
                                                                                      'specify --startIndex.')
                            node = tree_str.Node(label=label, branch_length=branch_length,
                                                 name=input_leaf_names[label - start_index])
                        else:
                            # An internal node
                            node = tree_str.Node(label=label, branch_length=branch_length)
                        label2node[label] = node
                    elif '->' in line:
                        # We are in the definition of a tree
                        try:
                            parent_index, child_index = get_parent_child(line)
                        except IOError:
                            exit(1)
                        parent_node = label2node[index2label[parent_index]]
                        child_node = label2node[index2label[child_index]]
                        parent_node.add_child(child_node)
                        child_node.set_parent(parent_node)
        # Find the root
        for label, node in label2node.items():
            node.set_attributes()
            if node.is_root():
                if tree.get_tree() is None:
                    tree.set_root(node)
                else:
                    raise IOError('Error! More than one root are detected. Errors in the tree structure.')
            if node.is_leaf():
                tree.add_to_leaf_nodes(node)
            else:
                tree.add_to_internal_nodes(node)
        tree.adjust_nodes_order()
        trees[input_tree_file] = tree
    return trees


def main():
    sys_args, input2output = parse_system_args()
    leaf_names = get_leaf_names(sys_args.leafNames)
    input2tree = read_trees(sys_args.input, sys_args.len, sys_args.startIndex, leaf_names)
    for input_file_name, output_file_name in input2output.items():
        input2tree[input_file_name].write_in_nexus(output_file_name)


if __name__ == '__main__':
    main()
