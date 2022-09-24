import argparse
import os
import shutil


def parse_sys_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(description='Update outdated files containing commands and options.')
    sys_args_parser.add_argument('--td', help='top directory of command files to be updated',
                                 metavar='TOP DIRECTORY', type=str)
    sys_args_parser.add_argument('--on', help='old command files\' name', metavar='OLD COMMAND FILES NAME', type=str)
    sys_args_parser.add_argument('--os', help='suffix to be appended to old command files\' name',
                                 metavar='SUFFIX TO OLD COMMAND FILES', type=str)
    sys_args_parser.add_argument('--ns', help='suffix to be appened to new command files\' name',
                                 metavar='SUFFIX TO NEW COMMAND FILES', type=str)
    sys_args_parser.add_argument('--stbr', help='strings originally in --on to be replaced',
                                 metavar='STRING TO BE REPLACED', type=str, nargs='*')
    sys_args_parser.add_argument('--osr', help='replacing strings for old command files',
                                 metavar='REPLACING STRING FOR OLD COMMAND FILES', type=str, nargs='*')
    sys_args_parser.add_argument('--nsr', help='replacing strings for new command files',
                                 metavar='REPLACING STRING FOR NEW COMMAND FILES', type=str, nargs='*')
    sys_args_parser.add_argument('--op', help='options to be removed in new command files originally in old command '
                                              'files', metavar='OPTIONS TO BE REMOVED', type=str, nargs='*')
    sys_args = sys_args_parser.parse_args()
    if sys_args.stbr is None or sys_args.osr is None or sys_args.nsr is None:
        if not (sys_args.stbr is None and sys_args.osr is None and sys_args.nsr is None):
            raise ValueError('ERROR! --stbr, --osr, and --nsr should have the same number of arguments.')
    elif len(sys_args.stbr) != len(sys_args.osr) or len(sys_args.stbr) != len(sys_args.nsr):
        raise ValueError('ERROR! The number of strings to be replaced must be the same as the number of replacing '
                         'strings for both old and new command files.')
    return sys_args


def targets_in_line(targets: list[str], line: str) -> list[str]:
    hits: list[str] = []
    results: list[str] = []
    indices: list[int] = []
    for target in targets:
        if ' -' + target + ' ' in line:
            hits.append(' -' + target + ' ')
            indices.append(line.index(hits[-1]))
        if ' --' + target + ' ' in line:
            hits.append(' --' + target + ' ')
            indices.append(line.index(hits[-1]))
    indices = sorted(range(len(indices)), key=lambda i: indices[i])
    for index in indices:
        results.append(hits[index])
    return results


def replace_strings_in_line(line: str, targets: list[str], substitution: list[str]) -> str:
    new_line: str = line
    if targets is not None and len(targets) > 0:
        for index in range(len(targets)):
            new_line = line.replace(targets[index], substitution[index])
    return new_line


def remove_options_in_line(positives: list[str], line: str) -> str:
    if positives is not None and len(positives) > 0:
        new_line: str = ''
        index: int = 0
        for positive in positives:
            positive_index: int = line.index(positive)
            if index < positive_index:
                if new_line != '':
                    new_line += ' '
                new_line += line[index:positive_index]
            index = line.index('-', positive_index + len(positive))
        new_line += ' ' + line[index:]
        return new_line
    else:
        return line


def update_command_files(sys_args: argparse.Namespace):
    for root, dirs, files in os.walk(sys_args.td):
        for file in files:
            if sys_args.on == file:
                shutil.copyfile(os.path.join(root, file), os.path.join(root, file + sys_args.ns))
                shutil.move(os.path.join(root, file), os.path.join(root, file + sys_args.os))

                if sys_args.stbr is not None:
                    _old_content: list[str] = []
                    with open(os.path.join(root, file + sys_args.os), 'r') as fh:
                        for line in fh:
                            _old_content.append(replace_strings_in_line(line, sys_args.stbr, sys_args.osr))

                    with open(os.path.join(root, file + sys_args.os), 'w') as fh:
                        for line in _old_content:
                            fh.write(line)

                if sys_args.stbr is not None or sys_args.op is not None:
                    _new_content: list[str] = []
                    with open(os.path.join(root, file + sys_args.ns), 'r') as fh:
                        for line in fh:
                            hits = targets_in_line(sys_args.op, line)
                            if len(hits) == 0:
                                _new_content.append(replace_strings_in_line(line, sys_args.stbr, sys_args.nsr))
                            else:
                                _new_content.append(replace_strings_in_line(remove_options_in_line(hits, line),
                                                                            sys_args.stbr, sys_args.nsr))

                    with open(os.path.join(root, file + sys_args.ns), 'w') as fh:
                        for line in _new_content:
                            fh.write(line)


def main() -> None:
    sys_args = parse_sys_args()
    update_command_files(sys_args)


if __name__ == "__main__":
    main()
