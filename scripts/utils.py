import os
import numpy as np
import re
from typing import Union

from snakemake import shell


NON_WORDS_PATTERN = '[^a-zA-Z0-9]+'
NONE_PATTERN = r'(?!.*)'


def which_patterns(
    s: str,
    patterns: list[str],
    min_pat_num: int = 0,
    max_pat_num: int = 1
) -> list[list[int]]:
    if max_pat_num <= 0:
        raise ValueError('max_pat_num must be greater than 0.')
    ret: list[list[int]] = []
    for __s in s:
        __ret: list[int] = []
        for id, pat in enumerate(patterns):
            if re.search(re.compile(pat), __s):
                __ret.append(id)
        if len(__ret) > max_pat_num:
            raise ValueError(f'Error! More than {max_pat_num} pattern(s) ({__ret}) is found in {__s}.')
        elif len(__ret) < min_pat_num:
            raise ValueError(f'Error! Less than {min_pat_num} pattern(s) ({__ret}) is found in {__s}.')
        
        if len(__ret) > 0:
            ret.append(__ret)
        else:
            ret.append([-1])

    return ret


vec_search = np.vectorize(lambda s, pat: re.search(pat, s), otypes=[bool])
vec_abspath = np.vectorize(lambda i: os.path.abspath(i), otypes=[str])
vec_pattern = np.vectorize(lambda i, pat: f'{pat}{i}{pat}', otypes=[str])


def get_matched_files(files: list[str], pattern: str) -> tuple[list[bool], list[str]]:
    idx = vec_search(
        files,
        pattern
    )
    matched = vec_abspath(sorted(np.array(files)[idx]))
    return idx, matched


def get_fine_tune_types(
    tool_fine_tune: dict,
    fine_tune_types: dict
) -> list[str]:
    ret: list[str] = []
    for i in get_dict_keys(tool_fine_tune):
        if i in fine_tune_types.keys():
            ret.append(fine_tune_types[i])
        else:
            raise ValueError(f'Error! {i} is not in {fine_tune_types}. Make sure you have listed all fine tune types under "benchmark - fineTuneTypes" in "config.yaml".')
    return ret


def get_data_types(
    used_data_types: list[str],
    available_data_types: dict,
    require_cnv: Union[bool, None] = None
) -> list[str]:
    ret: list[str] = []
    for i in used_data_types:
        if i in available_data_types.keys():
            if require_cnv is None or available_data_types[i]['requireCNV'] == require_cnv:
                ret.append(available_data_types[i]['dir'].rstrip('/'))
        else:
            raise ValueError(f'Error! {i} is not in {available_data_types.keys()}. Do not change any values under "benchmark - simulatedDataTypes"!')
    return ret


def get_cnv_keep_type(key: str, options: dict) -> int:
    ret: list[int] = []
    for i in options.keys():
        if key == options[i]['dir'].rstrip('/'):
            if options[i]['requireCNV']:
                if options[i]['includeCNV']:
                    ret.append(1)
                else:
                    ret.append(2)
    if len(ret) == 0:
        raise ValueError(f'Error! {key} is not in {options.keys()}.')
    elif len(ret) > 1:
        raise ValueError(f'Error! More than one key in {{options.keys()}} is found for {key}.')
    return ret[0]


def get_dict_keys(dict: dict) -> list:
    return list(dict.keys())


def get_dict_values(dict: dict, keys: list[str]):
    tmp = dict
    for key in keys:
        tmp = tmp[key]
    return tmp


def create_sifit_flag_file(file_name: str):
    if os.path.isfile(file_name):
        shell("echo \"Existing flag file detected: %s; skipping creation...\" >&2" % file_name)
        return
    else:
        shell("mkdir -p %s" % os.path.abspath(os.path.dirname(file_name)))
        shell("touch %s" % file_name)
        shell("echo \"Created: %s\" >&2" % file_name)


def generate_simulated_data(base_path,
                            tool_path,
                            config_file_path,
                            log_path,
                            results_path,
                            monovar_with_normal_relative_path,
                            monovar_without_normal_relative_path,
                            sciphi_with_normal_relative_path,
                            sciphi_without_normal_relative_path
                            ):
    """Generate simulated data.

    Args:
        base_path (string): Everything should be moved here.
        tool_path (string): Path to simulator.
        config_file_path (string): Path to simulation configuration file.
        log_path (string): Path to running log.
        results_path (string): Path to generated data defined in the configuration file.
        monovar_with_normal_relative_path (string): Relative to the most recent common parent folder of all simulated files.
        monovar_without_normal_relative_path (string): As above.
        sciphi_with_normal_relative_path (string): As above.
        sciphi_without_normal_relative_path (string): As above.
    """
    _overwrite_flag = False
    if os.path.isdir(results_path):
        # This is a sign of a failure previous run.
        shell("echo \"Existing directory detected: %s\" >&2" % results_path)
        shell("echo \"Removing...\" >&2")
        shell("rm -r %s" % results_path)
        _overwrite_flag = True
    if os.path.isdir(base_path):
        shell("echo \"Existing directory detected: %s\" >&2" % base_path)
        if _overwrite_flag:
            shell("echo \"Overwriting...\" >&2")
            shell("rm -r %s" % base_path)
            _generate_simulated_data(
                base_path, tool_path, config_file_path, log_path, results_path)
        else:
            shell("echo \"Skipping simulated data generation...\" >&2")
            return
    else:
        shell("echo \"Generating simulated data...\" >&2")
        _generate_simulated_data(base_path, tool_path,
                                 config_file_path,
                                 log_path,
                                 results_path
                                 )
    _edit_bam_file_names(base_path, monovar_with_normal_relative_path)
    _edit_bam_file_names(base_path, monovar_without_normal_relative_path)
    _edit_bam_file_names(base_path, sciphi_with_normal_relative_path)
    _edit_bam_file_names(base_path, sciphi_without_normal_relative_path)


def _generate_simulated_data(base_path,
                             tool_path,
                             config_file_path,
                             log_path,
                             results_path
                             ):
    """Generate simulated data. Should only be called by generate_simulated_data internally.

    Args:
        base_path (string): Everything should be moved here.
        tool_path (string): Path to simulator.
        config_file_path (string): Path to simulation configuration file.
        log_path (string): Path to running log.
        results_path (string): Path to generated data defined in the configuration file.
    """
    shell("mkdir -p %s" % base_path)
    shell("%s -F%s &>%s" % (tool_path, config_file_path, log_path))
    shell("cp %s %s" % (config_file_path, base_path))
    shell("mv %s %s" % (results_path, base_path))


def _edit_bam_file_names(base_path, bam_names_file_path):
    if os.path.isfile(bam_names_file_path):
        if not os.path.isabs(base_path):
            base_path = get_abs_working_dir() + base_path

        fh = open(bam_names_file_path, 'r')
        contents = fh.readlines()
        fh.close()

        fh = open(bam_names_file_path, 'w')
        for i in range(0, len(contents)):
            fh.write(base_path + contents[i])
        fh.close()
    else:
        shell("echo \"No such file exists.\" >&2")


def get_dataset_names(path):
    file_full_names = [f for f in os.listdir(
        path) if os.path.isfile(os.path.join(path, f))]
    dataset_names = []
    for f in file_full_names:
        comp = f.split(".")
        dataset_names.append(comp[-1])
    return dataset_names


def get_abs_working_dir():
    return os.getcwd() + "/"


def get_parameter_estimate_type(input):
    if input == "all":
        return ["median", "mean", "mode_gaussian"]
    elif input == "median":
        return ["median"]
    elif input == "mean":
        return ["mean"]
    elif input == "mode":
        return ["mode_gaussian"]
    else:
        shell("echo \"Unrecognized parameter estimater type. Should be one of 'all', 'median', 'mean', and 'mode'.\" >&2")


def get_line_num(file: str) -> int:
    if not os.path.isfile(os.path.abspath(file)):
        return -1

    line_num: int = 0
    with open(file, 'r') as fh:
        for line in fh:
            if not line.strip().startswith(' '):
                line_num += 1
    
    return line_num


def get_item_num(file: str) -> int:
    if not os.path.isfile(os.path.abspath(file)):
        return -1

    item_num: int = 0
    with open(file, 'r') as fh:
        for line in fh:
            if not line.strip().startswith(' '):
                item_num += len(line.strip().split())
    
    return item_num


def get_dir_name(file: str) -> str:
    return os.path.dirname(file) + '/'


def get_specific_file(files: list[str], key1: str, key2: str = '') -> str:
    # shell("echo \"%d\" &>2" % len(files))
    for file in files:
        # shell("echo \"%s\" &>2" % file)
        if key2 == '':
            if key1 in file:
                return file
        else:
            if key1 in file and key2 in file:
                return file
    raise ValueError('Error! ' + key1 + ' and ' + key2 + ' are not found.')


def get_content_files(file: str) -> list[str]:
    if os.path.isfile(file):
        results = []
        root = os.path.abspath(os.path.dirname(file))
        with open(file, 'r') as fh:
            for line in fh:
                content = os.path.join(root, line.strip())
                if os.path.isfile(content):
                    results.append(content)
        return results


def get_thread_num_for_sieve(rawDataFile: str, nrOfSitesPerThread: int) -> int:
    if not os.path.isfile(os.path.abspath(rawDataFile)):
        return -1

    target_line: bool = False
    with open(rawDataFile, 'r') as fh:
        for line in fh:
            if target_line:
                sites_num: int = int(line.strip())
                break
            if line.strip() == '=numCandidateMutatedSites=':
                target_line = True
    
    rawNrOfThreads: float = sites_num / nrOfSitesPerThread
    if int(rawNrOfThreads) == 0:
        # shell("echo \"%f\"" % rawNrOfThreads)
        return 1
    elif 0 <= rawNrOfThreads - int(rawNrOfThreads) < 0.5:
        # shell("echo \"%d\"" % int(rawNrOfThreads))
        return int(rawNrOfThreads)
    elif 0.5 <= rawNrOfThreads - int(rawNrOfThreads) <= 1:
        # shell("echo \"%d\"" % (int(rawNrOfThreads) + 1))
        return int(rawNrOfThreads) + 1


def get_thread_num_for_cellphy(dataFile: str, nrOfSitesPerThread: int) -> int:
    if not os.path.isfile(os.path.abspath(dataFile)):
        return -1
    
    target_line: bool = False
    sites_num: int = 0
    with open(dataFile, 'r') as fh:
        for line in fh:
            if line.strip():
                if target_line:
                    sites_num += 1
                elif line.strip().upper().startswith('#CHROM'):
                    target_line = True
    
    rawNrOfThreads: float = sites_num / nrOfSitesPerThread
    if int(rawNrOfThreads) == 0:
        return 1
    elif 0 <= rawNrOfThreads - int(rawNrOfThreads) < 0.5:
        return int(rawNrOfThreads)
    elif 0.5 <= rawNrOfThreads - int(rawNrOfThreads) <= 1:
        return int(rawNrOfThreads) + 1


def get_sifit_jar_abs_path(
    sifit_jar_name: str,
    lib_paths: list[str] = ['/usr/lib', '/usr/local/lib']
) -> str:
    if os.path.isfile(os.path.abspath(sifit_jar_name)):
        return os.path.abspath(sifit_jar_name)
    
    for p in lib_paths:
        if os.path.isfile(os.path.join(p, sifit_jar_name)):
            return os.path.join(p, sifit_jar_name)
    
    raise ValueError('Error! ' + sifit_jar_name + ' is not found.')

