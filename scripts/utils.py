import os
import os.path
from typing import List

from snakemake import shell

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
                            overwrite_flag,
                            monovar_with_normal_relative_path,
                            monovar_without_normal_relative_path,
                            sciphi_with_normal_relative_path,
                            sciphi_without_normal_relative_path
                            ):
    """Generate simulated data.

    Args:
        base_path (string): Everything should be moved here. If exists, choose to overwrite or not.
        tool_path (string): Path to simulator.
        config_file_path (string): Path to simulation configuration file.
        log_path (string): Path to running log.
        results_path (string): Path to generated data defined in the configuration file.
        overwrite_flag (bool): TRUE or FALSE.
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
        if _overwrite_flag | overwrite_flag:
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


def get_specific_file(files: List[str], key1: str, key2: str = '') -> str:
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


def get_content_files(file: str) -> List[str]:
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
    elif 0 < rawNrOfThreads - int(rawNrOfThreads) <= 0.5:
        # shell("echo \"%d\"" % int(rawNrOfThreads))
        return int(rawNrOfThreads)
    elif 0.5 < rawNrOfThreads - int(rawNrOfThreads) <= 1:
        # shell("echo \"%d\"" % (int(rawNrOfThreads) + 1))
        return int(rawNrOfThreads) + 1

