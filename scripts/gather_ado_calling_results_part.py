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
                if os.path.isfile(content):
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


def _get_ado_states(row_num: int, col_num: int, input: str) -> np.ndarray:
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
    true_ado_states: str,
    ado_states: str, 
    cell_names: str, 
    loci_info: str
) -> tuple[int, int, int, int, int, int, int, int, int, int, int, int, int, int]:
    _true_cell_names = _get_cell_names(true_cell_names)
    _cell_names = _get_cell_names(cell_names)
    if len(_true_cell_names) != len(_cell_names):
        raise ValueError('Error! Unmatched number of cells in ' + true_cell_names + ' and ' + cell_names)
    _cell_index_in_true = []
    for _cell_name in _cell_names:
        _cell_index_in_true.append(_true_cell_names.index(_cell_name))

    _loci_info = _get_loci_info(loci_info)
    _ado_states = _get_ado_states(len(_loci_info), len(_cell_names), ado_states)
    
    tp = fp = tn = fn = \
        tp_one_ado = fp_one_ado = tn_one_ado = fn_one_ado = \
            tp_two_ados = fp_two_ados = tn_two_ados = fn_two_ados = \
                t_two_ados_as_no_ado = t_two_ados_as_one_ado = 0

    with open(true_ado_states, 'r') as fh:
        for _line in fh:
            if not _line.startswith('#'):
                _comp = _line.strip().split()

                if (_comp[0], int(_comp[1])) in _loci_info:
                    _id = _loci_info.index((_comp[0], int(_comp[1])))

                    for _cell_name_index in range(0, len(_cell_names)):
                        _true_ado_state_value = int(_comp[2 + _cell_index_in_true[_cell_name_index]])
                        _ado_state_value = _ado_states[_id, _cell_name_index]

                        if _true_ado_state_value == 0:
                            if _ado_state_value == 0:
                                tn += 1
                                tn_one_ado += 1
                                tn_two_ados += 1
                            elif _ado_state_value == 1:
                                fp += 1
                                fp_one_ado += 1
                                tn_two_ados += 1
                            elif _ado_state_value == 2:
                                fp += 1
                                tn_one_ado += 1
                                fp_two_ados += 1
                            elif _ado_state_value == 3:
                                pass
                        elif _true_ado_state_value == 1:
                            if _ado_state_value == 0:
                                fn += 1
                                fn_one_ado += 1
                                tn_two_ados += 1
                            elif _ado_state_value == 1:
                                tp += 1
                                tp_one_ado += 1
                                tn_two_ados += 1
                            elif _ado_state_value == 2:
                                tp += 1
                                fn_one_ado += 1
                                fp_two_ados += 1
                            elif _ado_state_value == 3:
                                pass
                        elif _true_ado_state_value == 2:
                            if _ado_state_value == 0:
                                fn += 1
                                tn_one_ado += 1
                                fn_two_ados += 1
                                t_two_ados_as_no_ado += 1
                            elif _ado_state_value == 1:
                                tp += 1
                                fp_one_ado += 1
                                fn_two_ados += 1
                                t_two_ados_as_one_ado += 1
                            elif _ado_state_value == 2:
                                tp += 1
                                tn_one_ado += 1
                                tp_two_ados += 1
                            elif _ado_state_value == 3:
                                pass
                        elif _true_ado_state_value == 3:
                            if _ado_state_value == 0:
                                pass
                            elif _ado_state_value == 1:
                                pass
                            elif _ado_state_value == 2:
                                pass
                            elif _ado_state_value == 3:
                                pass

    return tp, fp, tn, fn, \
        tp_one_ado, fp_one_ado, tn_one_ado, fn_one_ado, \
            tp_two_ados, fp_two_ados, tn_two_ados, fn_two_ados, \
                t_two_ados_as_no_ado, t_two_ados_as_one_ado


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
    fine_tune_types = np.repeat(np.nan, len(snakemake.input['allFiles']))
elif len(snakemake.params['fineTuneType']) == 1:
    fine_tune_types = np.repeat(snakemake.params['fineTuneType'], len(snakemake.input['allFiles']))
else:
    fine_tune_types = []
    fine_tune_type_flag = True


tool_setup = []
data_type = []
dataset = []
cell_names = []
true_ado_states = []
loci_info = []
ado_states = []
TP = []
FP = []
TN = []
FN = []
recall = []
precision = []
fall_out = []
f1_score = []
TP_one_ado = []
FP_one_ado = []
TN_one_ado = []
FN_one_ado = []
recall_one_ado = []
precision_one_ado = []
fall_out_one_ado = []
f1_score_one_ado = []
TP_two_ados = []
FP_two_ados = []
TN_two_ados = []
FN_two_ados = []
recall_two_ados = []
precision_two_ados = []
fall_out_two_ados = []
f1_score_two_ados = []

for __tool_setup in sorted(snakemake.params['toolSetup']):
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
            
            __files = sorted([j for i in vec_get_content_files(
                get_matched_files(
                    snakemake.input['allFiles'],
                    __pat
                )[1]
            ) for j in i])

            inferred_idx, __loci_info = get_matched_files(
                __files,
                r'\.loci_info$'
            )
            cnt = np.sum(inferred_idx)

            if cnt == 0:
                print(f"Warning: no matches found for '.loci_info$' in {__files}.")

            tool_setup.extend(np.repeat(__tool_setup, cnt))
            data_type.extend(np.repeat(__data_type, cnt))
            dataset.extend(np.repeat(dataset_name, cnt))
            true_ado_states.extend(np.repeat(os.path.abspath(snakemake.params['trueAdoStatesPre'] + '.' + dataset_name), cnt))

            loci_info.extend(__loci_info)

            _, __cell_names = get_matched_files(
                __files,
                r'\.cell_names$'
            )
            cell_names.extend(__cell_names)

            _, __ado_states = get_matched_files(
                __files,
                r'\.ado$'
            )
            ado_states.extend(__ado_states)

            for idx in range(len(__loci_info)):
                tp, fp, tn, fn, \
                    tp_one_ado, fp_one_ado, tn_one_ado, fn_one_ado, \
                        tp_two_ados, fp_two_ados, tn_two_ados, fn_two_ados, \
                            t_two_ados_as_no_ado, t_two_ados_as_one_ado = _get_properties(
                                os.path.abspath(snakemake.input['trueCellNames']), 
                                true_ado_states[-1],
                                __ado_states[idx], 
                                __cell_names[idx], 
                                __loci_info[idx]
                            )

                TP.append(tp)
                FP.append(fp)
                TN.append(tn)
                FN.append(fn)

                _precision, _recall, _f1_score, _fall_out = _get_metrics(float(tp), float(fp), float(tn), float(fn))
                recall.append(_recall)
                precision.append(_precision)
                f1_score.append(_f1_score)
                fall_out.append(_fall_out)

                TP_one_ado.append(tp_one_ado)
                FP_one_ado.append(fp_one_ado)
                TN_one_ado.append(tn_one_ado)
                FN_one_ado.append(fn_one_ado)

                _precision_one_ado, _recall_one_ado, _f1_score_one_ado, _fall_out_one_ado = _get_metrics(float(tp_one_ado), float(fp_one_ado), float(tn_one_ado), float(fn_one_ado))
                recall_one_ado.append(_recall_one_ado)
                precision_one_ado.append(_precision_one_ado)
                f1_score_one_ado.append(_f1_score_one_ado)
                fall_out_one_ado.append(_fall_out_one_ado)

                TP_two_ados.append(tp_two_ados)
                FP_two_ados.append(fp_two_ados)
                TN_two_ados.append(tn_two_ados)
                FN_two_ados.append(fn_two_ados)

                _precision_two_ados, _recall_two_ados, _f1_score_two_ados, _fall_out_two_ados = _get_metrics(float(tp_two_ados), float(fp_two_ados), float(tn_two_ados), float(fn_two_ados))
                recall_two_ados.append(_recall_two_ados)
                precision_two_ados.append(_precision_two_ados)
                f1_score_two_ados.append(_f1_score_two_ados)
                fall_out_two_ados.append(_fall_out_two_ados)

            if fine_tune_type_flag:
                fine_tune_type_idx = which_patterns(
                    __loci_info,
                    vec_pattern(
                        snakemake.params['fineTuneType'],
                        NON_WORDS_PATTERN
                    )
                )
                fine_tune_types.extend([snakemake.params['fineTuneType'][j] for i in fine_tune_type_idx for j in i if j >= 0])


ado_collection = pd.DataFrame(
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
        'loci_info': pd.Series(loci_info),
        'true_ado_states': pd.Series(true_ado_states),
        'ado_states': pd.Series(ado_states),
        'true_positive': pd.Series(TP),
        'false_positive': pd.Series(FP),
        'true_negative': pd.Series(TN),
        'false_negative': pd.Series(FN),
        'recall': pd.Series(recall),
        'precision': pd.Series(precision),
        'f1_score': pd.Series(f1_score),
        'fall_out': pd.Series(fall_out),
        'true_positive_one_ado': pd.Series(TP_one_ado),
        'false_positive_one_ado': pd.Series(FP_one_ado),
        'true_negative_one_ado': pd.Series(TN_one_ado),
        'false_negative_one_ado': pd.Series(FN_one_ado),
        'recall_one_ado': pd.Series(recall_one_ado),
        'precision_one_ado': pd.Series(precision_one_ado),
        'f1_score_one_ado': pd.Series(f1_score_one_ado),
        'fall_out_one_ado': pd.Series(fall_out_one_ado),
        'true_positive_two_ados': pd.Series(TP_two_ados),
        'false_positive_two_ados': pd.Series(FP_two_ados),
        'true_negative_two_ados': pd.Series(TN_two_ados),
        'false_negative_two_ados': pd.Series(FN_two_ados),
        'recall_two_ados': pd.Series(recall_two_ados),
        'precision_two_ados': pd.Series(precision_two_ados),
        'f1_score_two_ados': pd.Series(f1_score_two_ados),
        'fall_out_two_ados': pd.Series(fall_out_two_ados)
    }
)

ado_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
