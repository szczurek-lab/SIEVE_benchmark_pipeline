import sys
import pandas as pd
import numpy as np
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from scripts.utils import NON_WORDS_PATTERN, which_patterns, vec_pattern, get_matched_files


def _get_content_files(file: str) -> list[str]:
    if os.path.isfile(file):
        results = []
        root = os.path.abspath(os.path.dirname(file))
        with open(file, 'r') as fh:
            for line in fh:
                content = os.path.join(root, line.strip())
                if content != os.path.abspath(file) and os.path.isfile(content):
                    results.append(content)
        return results


vec_get_content_files = np.vectorize(_get_content_files, otypes=[list])


def _get_cell_names(input: str) -> list[str]:
    if os.path.isfile(input):
        results = []
        with open(input, 'r') as fh:
            for line in fh:
                results.append(line.strip())
        return results


def _get_loci_info(input: str) -> list[tuple[str, int]]:
    if os.path.isfile(input):
        results = []
        with open(input, 'r') as fh:
            for line in fh:
                comp = line.strip().split()
                results.append((comp[0].strip(), int(comp[1].strip())))
        return results


def _get_ternary(row_num: int, col_num: int, input: str) -> np.ndarray:
    if os.path.isfile(input):
        results = np.empty((row_num, col_num), dtype=np.int32)
        with open(input, 'r') as fh:
            row_index = 0
            for line in fh:
                if row_num <= row_index:
                    raise ValueError('Error! Unmatched size. Expecting ' + str(row_num) + ' rows, detecting ' + str(row_index + 1) + ' rows.')
                comp = line.strip().split()
                if col_num != len(comp):
                    raise ValueError('Error! Unmatched size. Expecting ' + str(col_num) + ' columns, detecting ' + str(len(comp)) + ' columns.')
                for col_index in range(0, len(comp)):
                    results[row_index, col_index] = int(comp[col_index])
                row_index += 1
        return results


def _get_properties(
    true_cell_names: str, 
    true_loci_info: str, 
    true_ternary: str, 
    cell_names: str, 
    loci_info: str, 
    ternary: str
) -> tuple[int, int, int, int, int, int, int, int, int, int, int, int, int, int]:
    _true_cell_names = _get_cell_names(true_cell_names)
    _cell_names = _get_cell_names(cell_names)
    if len(_true_cell_names) != len(_cell_names):
        raise ValueError('Error! Unmatched number of cells in ' + true_cell_names + ' and ' + cell_names)
    _cell_index_in_true = []
    for _cell_name in _cell_names:
        _cell_index_in_true.append(_true_cell_names.index(_cell_name))

    _true_loci_info = _get_loci_info(true_loci_info)
    _loci_info = _get_loci_info(loci_info)
    
    _true_ternary = _get_ternary(len(_true_loci_info), len(_true_cell_names), true_ternary)
    _ternary = _get_ternary(len(_loci_info), len(_cell_names), ternary)
    
    detected_true_loci_info = []
    tp = fp = tn = fn = \
        tp_hetero_mu = fp_hetero_mu = tn_hetero_mu = fn_hetero_mu = \
            tp_homo_mu = fp_homo_mu = tn_homo_mu = fn_homo_mu = \
                t_homo_mu_as_homo_ref = t_homo_mu_as_hetero_mu = \
                    missing_of_homo_ref = missing_of_hetero_mu = missing_of_homo_mu = 0
    
    for _loci_info_index in range(0, len(_loci_info)):
        if _loci_info[_loci_info_index] in _true_loci_info:
            detected_true_loci_info.append(_loci_info[_loci_info_index])
            _true_loci_info_index = _true_loci_info.index(_loci_info[_loci_info_index])
            for _cell_name_index in range(0, len(_cell_names)):
                _true_ternary_value = _true_ternary[_true_loci_info_index, _cell_index_in_true[_cell_name_index]]
                _ternary_value = _ternary[_loci_info_index, _cell_name_index]
                
                if _true_ternary_value == 0:
                    if _ternary_value == 0:
                        tn += 1
                        tn_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 1:
                        fp += 1
                        fp_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 2:
                        fp += 1
                        tn_hetero_mu += 1
                        fp_homo_mu += 1
                    elif _ternary_value == 3:
                        missing_of_homo_ref += 1
                elif _true_ternary_value == 1:
                    if _ternary_value == 0:
                        fn += 1
                        fn_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 1:
                        tp += 1
                        tp_hetero_mu += 1
                        tn_homo_mu += 1
                    elif _ternary_value == 2:
                        tp += 1
                        fn_hetero_mu += 1
                        fp_homo_mu += 1
                    elif _ternary_value == 3:
                        missing_of_hetero_mu += 1
                elif _true_ternary_value == 2:
                    if _ternary_value == 0:
                        fn += 1
                        tn_hetero_mu += 1
                        fn_homo_mu += 1
                        t_homo_mu_as_homo_ref += 1
                    elif _ternary_value == 1:
                        tp += 1
                        fp_hetero_mu += 1
                        fn_homo_mu += 1
                        t_homo_mu_as_hetero_mu += 1
                    elif _ternary_value == 2:
                        tp += 1
                        tn_hetero_mu += 1
                        tp_homo_mu += 1
                    elif _ternary_value == 3:
                        missing_of_homo_mu += 1
                elif _true_ternary_value == 3:
                    if _ternary_value == 0:
                        pass
                    elif _ternary_value == 1:
                        pass
                    elif _ternary_value == 2:
                        pass
                    elif _ternary_value == 3:
                        pass
                elif _true_ternary_value == -3:
                    if _ternary_value == 0:
                        pass
                    elif _ternary_value == 1:
                        pass
                    elif _ternary_value == 2:
                        pass
                    elif _ternary_value == 3:
                        pass
                elif _true_ternary_value == -2:
                    if _ternary_value == 0:
                        fn += 1
                    elif _ternary_value == 1:
                        tp += 1
                    elif _ternary_value == 2:
                        tp += 1
                    elif _ternary_value == 3:
                        pass
                elif _true_ternary_value == -1:
                    if _ternary_value == 0:
                        tn += 1
                    elif _ternary_value == 1:
                        fp += 1
                    elif _ternary_value == 2:
                        fp += 1
                    elif _ternary_value == 3:
                        pass
        else:
            for _cell_name_index in range(0, len(_cell_names)):
                _ternary_value = _ternary[_loci_info_index, _cell_name_index]
                if _ternary_value == 0:
                    tn += 1
                    tn_hetero_mu += 1
                    tn_homo_mu += 1
                elif _ternary_value == 1:
                    fp += 1
                    fp_hetero_mu += 1
                    tn_homo_mu += 1
                elif _ternary_value == 2:
                    fp += 1
                    tn_hetero_mu += 1
                    fp_homo_mu += 1
                elif _ternary_value == 3:
                    missing_of_homo_ref += 1
    
    for _true_loci_info_index in range(0, len(_true_loci_info)):
        if _true_loci_info[_true_loci_info_index] not in detected_true_loci_info:
            for _true_cell_name_index in range(0, len(_true_cell_names)):
                _true_ternary_value = _true_ternary[_true_loci_info_index, _true_cell_name_index]
                if _true_ternary_value == 0:
                    tn += 1
                    tn_hetero_mu += 1
                    tn_homo_mu += 1
                elif _true_ternary_value == 1:
                    fn += 1
                    fn_hetero_mu += 1
                    tn_homo_mu += 1
                elif _true_ternary_value == 2:
                    fn += 1
                    tn_hetero_mu += 1
                    fn_homo_mu += 1
                    t_homo_mu_as_homo_ref += 1
                elif _true_ternary_value == 3:
                    pass
                elif _true_ternary_value == -3:
                    pass
                elif _true_ternary_value == -2:
                    fn += 1
                elif _true_ternary_value == -1:
                    tn += 1
    
    return \
        tp, fp, tn, fn, \
            tp_hetero_mu, fp_hetero_mu, tn_hetero_mu, fn_hetero_mu, \
                tp_homo_mu, fp_homo_mu, tn_homo_mu, fn_homo_mu, t_homo_mu_as_homo_ref, t_homo_mu_as_hetero_mu, \
                    missing_of_homo_ref, missing_of_hetero_mu, missing_of_homo_mu, \
                        len(_true_loci_info) * len(_true_cell_names)


def _get_metrics(
    true_positive: float, 
    false_positive: float, 
    true_negative: float, 
    false_negative: float
) -> tuple[float, float, float, float]:
    if true_positive + false_positive > 0.0:
        _precision = true_positive / (true_positive + false_positive)
    else:
        _precision = float('nan')
    
    if true_positive + false_negative > 0.0:
        _recall = true_positive / (true_positive + false_negative)
    else:
        _recall = float('nan')
    
    if _precision + _recall > 0.0:
        _f1_score = 2 * _precision * _recall / (_precision + _recall)
    else:
        _f1_score = float('nan')

    if true_negative + false_positive > 0.0:
        _fall_out = false_positive / (true_negative + false_positive)
    else:
        _fall_out = float('nan')
    
    return _precision, _recall, _f1_score, _fall_out


search_for_tool_setup_flag = False
if snakemake.params['toolSetup'] is None or type(snakemake.params['toolSetup']) == str:
    tool_setups = [snakemake.params['toolSetup']]
elif type(snakemake.params['toolSetup']) == list:
    search_for_tool_setup_flag = True
    tool_setups = sorted(snakemake.params['toolSetup'])
else:
    raise ValueError(f'toolSetup must be a string or a list of strings ({snakemake.params["toolSetup"]} was given).')

fine_tune_type_flag = False
if snakemake.params['fineTuneType'] is None:
    fine_tune_types = np.repeat(np.nan, len(snakemake.input['allFiles']) if snakemake.params['forSieve'] else len(snakemake.input['lociInfo']))
elif len(snakemake.params['fineTuneType']) == 1:
    fine_tune_types = np.repeat(snakemake.params['fineTuneType'], len(snakemake.input['allFiles']) if snakemake.params['forSieve'] else len(snakemake.input['lociInfo']))
else:
    fine_tune_types = []
    fine_tune_type_flag = True


tool_setup = []
data_type = []
dataset = []
cell_names = []
loci_info = []
genotypes = []
ternary = []
genotype_probs = []
true_loci_info = []
true_genotypes = []
true_ternary = []
TP = []
FP = []
TN = []
FN = []
recall = []
precision = []
fall_out = []
f1_score = []
TP_hetero_mu = []
FP_hetero_mu = []
TN_hetero_mu = []
FN_hetero_mu = []
recall_hetero_mu = []
precision_hetero_mu = []
fall_out_hetero_mu = []
f1_score_hetero_mu = []
TP_homo_mu = []
FP_homo_mu = []
TN_homo_mu = []
FN_homo_mu = []
recall_homo_mu = []
precision_homo_mu = []
fall_out_homo_mu = []
f1_score_homo_mu = []
true_homo_mu = []
true_homo_mu_as_homo_ref = []
true_homo_mu_as_hetero_mu = []
prop_of_missing = []

for __tool_setup in tool_setups:
    for __data_type in sorted(snakemake.params['dataType']):
        for dataset_name in sorted(snakemake.params['datasetNames']):
            if search_for_tool_setup_flag:
                __pat = ''.join(
                    [f'(?=.*{i})' for i in vec_pattern([__tool_setup, __data_type, dataset_name], NON_WORDS_PATTERN)]
                )
            else:
                __pat = ''.join(
                    [f'(?=.*{i})' for i in vec_pattern([__data_type, dataset_name], NON_WORDS_PATTERN)]
                )

            if snakemake.params['forSieve']:
                __files = sorted([j for i in vec_get_content_files(
                    get_matched_files(
                        snakemake.input['allFiles'],
                        __pat
                    )[1]
                ) for j in i])

                inferred_idx, __loci_info = get_matched_files(
                    __files,
                    r'merged\.loci_info$' if snakemake.params['mergedDuplicateLines'] else r'\.loci_info$'
                )
                cnt = np.sum(inferred_idx)

                if cnt == 0:
                    print(f"Warning: no matches found for {'merged.loci_info$' if snakemake.params['mergedDuplicateLines'] else '.loci_info$'} in {__files}.")

                tool_setup.extend(np.repeat(__tool_setup, cnt))
                data_type.extend(np.repeat(__data_type, cnt))
                dataset.extend(np.repeat(dataset_name, cnt))

                loci_info.extend(__loci_info)

                _, __cell_names = get_matched_files(
                    __files,
                    r'\.cell_names$'
                )
                cell_names.extend(__cell_names)
                
                genotypes.extend(
                    get_matched_files(
                        __files,
                        r'\.genotypes$'
                    )[1]
                )

                _, __ternary = get_matched_files(
                    __files,
                    r'merged\.ternary$' if snakemake.params['mergedDuplicateLines'] else r'\.ternary$'
                )
                ternary.extend(__ternary)

                genotype_probs.extend(
                    get_matched_files(
                        __files,
                        r'\.probs$'
                    )[1]
                )
            else:
                inferred_idx, __loci_info = get_matched_files(
                    snakemake.input['lociInfo'],
                    __pat
                )
                cnt = np.sum(inferred_idx)

                if cnt == 0:
                    print(f"Warning: no matches found for {__pat} in {snakemake.input['inferredTrees']}.")

                tool_setup.extend(np.repeat(__tool_setup, cnt))
                data_type.extend(np.repeat(__data_type, cnt))
                dataset.extend(np.repeat(dataset_name, cnt))

                loci_info.extend(__loci_info)

                _, __cell_names = get_matched_files(
                    snakemake.input['cellNames'],
                    __pat
                )
                cell_names.extend(__cell_names)
                
                genotypes.extend(
                    get_matched_files(
                        snakemake.input['genotypes'],
                        __pat
                    )[1]
                )

                _, __ternary = get_matched_files(
                    snakemake.input['ternary'],
                    __pat
                )
                ternary.extend(__ternary)

                genotype_probs.extend(
                    get_matched_files(
                        snakemake.input['genotypeProbs'],
                        __pat
                    )[1]
                )

            if len(__loci_info) != len(__cell_names) != len(__ternary):
                raise ValueError(f'lociInfo, cellNames, and ternary must have the same length ({len(__loci_info)} != {len(__cell_names)} != {len(__ternary)}).')

            __true_pat = ''.join(
                [f'(?=.*{i})' for i in vec_pattern([__data_type, dataset_name], NON_WORDS_PATTERN)]
            )

            true_loci_info_idx, __true_loci_info = get_matched_files(
                snakemake.input['trueSNVSites'],
                __true_pat
            )
            if np.sum(true_loci_info_idx) == 0 or np.sum(true_loci_info_idx) > 1:
                raise ValueError(f"Pattern {__true_pat} should match exactly one object, but {__true_loci_info} found.")
            true_loci_info.extend(
                np.repeat(
                    __true_loci_info,
                    cnt
                )
            )

            true_genotypes_idx, __true_genotypes = get_matched_files(
                snakemake.input['trueSNVGenotypesNU'],
                __true_pat
            )
            if np.sum(true_genotypes_idx) == 0 or np.sum(true_genotypes_idx) > 1:
                raise ValueError(f"Pattern {__true_pat} should match exactly one object, but {__true_genotypes} found.")
            true_genotypes.extend(
                np.repeat(
                    __true_genotypes,
                    cnt
                )
            )

            true_ternary_idx, __true_ternary = get_matched_files(
                snakemake.input['trueSNVGenotypesTer'],
                __true_pat
            )
            if np.sum(true_ternary_idx) == 0 or np.sum(true_ternary_idx) > 1:
                raise ValueError(f"Pattern {__true_pat} should match exactly one object, but {__true_ternary} found.")
            true_ternary.extend(
                np.repeat(
                    __true_ternary,
                    cnt
                )
            )

            for idx in range(len(__loci_info)):
                tp, fp, tn, fn, \
                    tp_hetero_mu, fp_hetero_mu, tn_hetero_mu, fn_hetero_mu, \
                        tp_homo_mu, fp_homo_mu, tn_homo_mu, fn_homo_mu, t_homo_mu_as_homo_ref, t_homo_mu_as_hetero_mu, \
                            missing_of_homo_ref, missing_of_hetero_mu, missing_of_homo_mu, \
                                total_entries = _get_properties(
                                    os.path.abspath(snakemake.input['trueCellNames']),
                                    true_loci_info[-1],
                                    true_ternary[-1],
                                    __cell_names[idx],
                                    __loci_info[idx],
                                    __ternary[idx]
                                )
                
                TP.append(tp)
                FP.append(fp)
                TN.append(tn + missing_of_homo_ref)
                FN.append(fn)
                TP_hetero_mu.append(tp_hetero_mu)
                FP_hetero_mu.append(fp_hetero_mu)
                TN_hetero_mu.append(tn_hetero_mu)
                FN_hetero_mu.append(fn_hetero_mu + missing_of_hetero_mu)
                TP_homo_mu.append(tp_homo_mu)
                FP_homo_mu.append(fp_homo_mu)
                TN_homo_mu.append(tn_homo_mu)
                FN_homo_mu.append(fn_homo_mu + missing_of_homo_mu)
                true_homo_mu.append(tp_homo_mu + fn_homo_mu + missing_of_homo_mu)
                true_homo_mu_as_homo_ref.append(t_homo_mu_as_homo_ref)
                true_homo_mu_as_hetero_mu.append(t_homo_mu_as_hetero_mu)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp), 
                    float(fp), 
                    float(tn + missing_of_homo_ref), 
                    float(fn)
                )
                recall.append(_recall)
                precision.append(_precision)
                f1_score.append(_f1_score)
                fall_out.append(_fall_out)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_hetero_mu), 
                    float(fp_hetero_mu), 
                    float(tn_hetero_mu), 
                    float(fn_hetero_mu + missing_of_hetero_mu)
                )
                recall_hetero_mu.append(_recall)
                precision_hetero_mu.append(_precision)
                f1_score_hetero_mu.append(_f1_score)
                fall_out_hetero_mu.append(_fall_out)

                # if (not snakemake.params['mergedDuplicateLines']) and snakemake.params['tool'].lower() == 'sciphi':
                #     recall_homo_mu.append(np.nan)
                #     precision_homo_mu.append(np.nan)
                #     f1_score_homo_mu.append(np.nan)
                #     fall_out_homo_mu.append(np.nan)
                # else:
                #     _precision, _recall, _f1_score, _fall_out = _get_metrics(float(tp_homo_mu), float(fp_homo_mu), float(tn_homo_mu), float(fn_homo_mu))
                #     recall_homo_mu.append(_recall)
                #     precision_homo_mu.append(_precision)
                #     f1_score_homo_mu.append(_f1_score)
                #     fall_out_homo_mu.append(_fall_out)
                
                _precision, _recall, _f1_score, _fall_out = _get_metrics(
                    float(tp_homo_mu), 
                    float(fp_homo_mu), 
                    float(tn_homo_mu), 
                    float(fn_homo_mu + missing_of_homo_mu)
                )
                recall_homo_mu.append(_recall)
                precision_homo_mu.append(_precision)
                f1_score_homo_mu.append(_f1_score)
                fall_out_homo_mu.append(_fall_out)
            
                prop_of_missing.append((missing_of_homo_ref + missing_of_hetero_mu + missing_of_homo_mu) / total_entries)
            
            if fine_tune_type_flag:
                fine_tune_type_idx = which_patterns(
                    __loci_info,
                    vec_pattern(
                        snakemake.params['fineTuneType'],
                        NON_WORDS_PATTERN
                    )
                )
                fine_tune_types.extend([snakemake.params['fineTuneType'][j] for i in fine_tune_type_idx for j in i if j >= 0])


variant_collection = pd.DataFrame(
    {
        'cell_num': snakemake.params['cellNum'],
        'coverage_mean': snakemake.params["covMean"],
        'coverage_variance': snakemake.params["covVariance"],
        'eff_seq_err_rate': snakemake.params['effSeqErrRate'],
        'ado_rate': snakemake.params['adoRate'],
        'deletion_rate': snakemake.params['deletionRate'],
        'insertion_rate': snakemake.params['insertionRate'],
        'dataset': pd.Series(dataset),
        'tool': snakemake.params['tool'],
        'snv_type': snakemake.params['snvType'],
        'tool_setup': pd.Series(tool_setup),
        'fine_tune_type': pd.Series(fine_tune_types),
        'data_type': pd.Series(data_type),
        'true_cell_names': os.path.abspath(snakemake.input['trueCellNames']),
        'cell_names': pd.Series(cell_names),
        'true_loci_info': pd.Series(true_loci_info),
        'loci_info': pd.Series(loci_info),
        'true_genotypes': pd.Series(true_genotypes),
        'genotypes': pd.Series(genotypes),
        'true_ternary': pd.Series(true_ternary),
        'ternary': pd.Series(ternary),
        'genotype_probs': pd.Series(genotype_probs),
        'true_positive': pd.Series(TP),
        'false_positive': pd.Series(FP),
        'true_negative': pd.Series(TN),
        'false_negative': pd.Series(FN),
        'recall': pd.Series(recall),
        'precision': pd.Series(precision),
        'f1_score': pd.Series(f1_score),
        'fall_out': pd.Series(fall_out),
        'true_positive_hetero_mu': pd.Series(TP_hetero_mu),
        'false_positive_hetero_mu': pd.Series(FP_hetero_mu),
        'true_negative_hetero_mu': pd.Series(TN_hetero_mu),
        'false_negative_hetero_mu': pd.Series(FN_hetero_mu),
        'recall_hetero_mu': pd.Series(recall_hetero_mu),
        'precision_hetero_mu': pd.Series(precision_hetero_mu),
        'f1_score_hetero_mu': pd.Series(f1_score_hetero_mu),
        'fall_out_hetero_mu': pd.Series(fall_out_hetero_mu),
        'true_positive_homo_mu': pd.Series(TP_homo_mu),
        'false_positive_homo_mu': pd.Series(FP_homo_mu),
        'true_negative_homo_mu': pd.Series(TN_homo_mu),
        'false_negative_homo_mu': pd.Series(FN_homo_mu),
        'recall_homo_mu': pd.Series(recall_homo_mu),
        'precision_homo_mu': pd.Series(precision_homo_mu),
        'f1_score_homo_mu': pd.Series(f1_score_homo_mu),
        'fall_out_homo_mu': pd.Series(fall_out_homo_mu),
        'true_homo_mu': pd.Series(true_homo_mu),
        'true_homo_mu_as_homo_ref': pd.Series(true_homo_mu_as_homo_ref),
        'true_homo_mu_as_hetero_mu': pd.Series(true_homo_mu_as_hetero_mu),
        'prop_of_missing': pd.Series(prop_of_missing)
    }
)

variant_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
