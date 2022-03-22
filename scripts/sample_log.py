import argparse
import os
from typing import List


def parse_system_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(description='Sample log file from BEAST 2.')
    sys_args_parser.add_argument(
        '--log',
        metavar='LOG',
        type=str,
        required=True,
        help='log file from BEAST 2'
    )
    sys_args_parser.add_argument(
        '--step',
        metavar='STEP',
        type=int,
        required=True,
        help='sample step; one sample is kept for every "--step" samples'
    )
    sys_args_parser.add_argument(
        '--out',
        metavar='OUTPUT',
        type=str,
        required=True,
        help='output file'
    )
    sys_args = sys_args_parser.parse_args()
    if not os.path.isfile(sys_args.log):
        raise ValueError('ERROR! File does not exist: ' + sys_args.log)
    if sys_args.step <= 1:
        raise ValueError('ERROR! Invalid value for "--step".')
    return sys_args


def parse_log(log_file: str, step: int) -> List[str]:
    ret: List[str] = []
    with open(log_file, 'r') as fh:
        reached_log: bool = False
        count: int = step
        for line in fh:
            line = line.strip()
            if reached_log:
                if count == 1:
                    count = step
                    continue
                elif count == step:
                    ret.append(line)
                count -= 1
            else:
                ret.append(line)
                if line.startswith("Sample"):
                    reached_log = True
    return ret


def write_log(contents: List[str], out_file: str) -> None:
    with open(out_file, 'w') as fh:
        for line in contents:
            fh.write(line)
            fh.write('\n')


def main() -> None:
    sys_args: argparse.Namespace = parse_system_args()
    contents: List[str] = parse_log(sys_args.log, sys_args.step)
    write_log(contents, sys_args.out)


if __name__ == '__main__':
    main()
