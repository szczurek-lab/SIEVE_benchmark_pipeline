import _io
import operator
import os
from typing import Optional


class Node:
    def __init__(self, label: int, branch_length: float, name: str = None, parent_node: 'Node' = None,
                 child1: 'Node' = None, child2: 'Node' = None):
        self._label = label
        self._offset_label = -1
        self._branch_length = branch_length
        self._name = name
        self._parent_node = parent_node
        self._children = []
        if child1 is not None:
            self._children.append(child1)
        if child2 is not None:
            self._children.append(child2)
        self._is_leaf = len(self._children) == 0
        self._is_root = parent_node is None

    def get_label(self) -> int:
        return self._label

    def get_parent(self) -> 'Node':
        return self._parent_node

    def get_child(self, child_index: int) -> Optional['Node']:
        if child_index + 1 > len(self._children):
            return None
        return self._children[child_index]

    def get_branch_length(self) -> float:
        return self._branch_length

    def set_node_name(self, name: str) -> None:
        self._name = name

    def get_node_name(self) -> str:
        return self._name

    def set_parent(self, parent: 'Node') -> None:
        self._parent_node = parent
        self._is_root = parent is None

    def add_child(self, child: 'Node') -> None:
        if child is not None:
            for _child in self._children:
                if child is _child:
                    return
                elif child.get_label() == _child.get_label():
                    raise IOError('Duplicate children with the same label (' + str(child.get_label()) + ') detected. '
                                                                                                        'Errors in the '
                                                                                                        'node '
                                                                                                        'definition.')
            self._children.append(child)
            self._is_leaf = len(self._children) == 0

    def set_attributes(self) -> None:
        if self._parent_node is None:
            self._is_root = True
            self._branch_length = 0.0
        else:
            self._is_root = False
        if len(self._children) == 0:
            self._is_leaf = True
        else:
            self._is_leaf = False

    def is_leaf(self) -> bool:
        return self._is_leaf

    def is_root(self) -> bool:
        return self._is_root

    def set_offset_label(self, offset: int) -> None:
        if self.is_leaf():
            self._offset_label = self._label - offset

    def get_offset_label(self) -> int:
        return self._offset_label

    def write_in_nexus(self, fh: _io.TextIOWrapper) -> None:
        if self.is_leaf():
            fh.write('{:d}:{:f}'.format(self._offset_label, self._branch_length))
        else:
            fh.write('(')
            for index in range(len(self._children)):
                self._children[index].write_in_nexus(fh)
                if index < len(self._children) - 1:
                    fh.write(',')
            fh.write('):{0:f}'.format(self._branch_length))


class Tree:
    def __init__(self, root: Node = None):
        self._root = root
        self._internal_nodes = []
        self._leaf_nodes = []

    def get_tree(self) -> Node:
        return self._root

    def set_root(self, root: Node) -> None:
        self._root = root

    def add_to_leaf_nodes(self, node: Node) -> None:
        if node is not None:
            self._leaf_nodes.append(node)

    def add_to_internal_nodes(self, node: Node) -> None:
        if node is not None:
            self._internal_nodes.append(node)

    def adjust_nodes_order(self) -> None:
        self._leaf_nodes = sorted(self._leaf_nodes, key=operator.attrgetter('_label'))
        self._internal_nodes = sorted(self._internal_nodes, key=operator.attrgetter('_label'))
        offset = self._leaf_nodes[0].get_label() - 1
        if offset != 0:
            for node in self._leaf_nodes:
                node.set_offset_label(offset)

    def write_in_nexus(self, output: str) -> None:
        if self._root is None:
            raise ValueError('Error! The tree is empty.')

        if not os.path.exists(os.path.dirname(output)):
            os.makedirs(os.path.dirname(output))

        # Write header
        with open(output, 'w') as fh:
            fh.write('#NEXUS\n\n')

            # Write taxa block
            fh.write('Begin taxa;\n')
            fh.write('\tDimensions ntax=' + str(len(self._leaf_nodes)) + ";\n")
            fh.write('\t\tTaxlabels\n')
            for node in self._leaf_nodes:
                fh.write('\t\t\t' + node.get_node_name() + '\n')
            fh.write('\t\t\t;\n')
            fh.write('End;\n')

            # Write trees block
            fh.write('Begin trees;\n')
            fh.write('\tTranslate\n')
            for node in self._leaf_nodes:
                fh.write('\t\t{0:4d} {1:s},\n'.format(node.get_offset_label(), node.get_node_name()))
            fh.write(';\n')
            fh.write('tree TREE1 = ')
            self._root.write_in_nexus(fh)
            fh.write(';\n')
            fh.write('End;\n')
