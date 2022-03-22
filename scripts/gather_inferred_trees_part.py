import pandas as pd
import os


dataset = []
inferred_trees_file_names = []
true_trees_file_names = []

if snakemake.params["forSieve"] is False:
    for dataset_name in sorted(snakemake.params['datasetNames']):
        dataset.append(dataset_name)
        true_trees_file_names.append(os.path.abspath(
            snakemake.params['trueTreesDir'] + snakemake.params['trueTreesPre'] + '.' + dataset_name))
        for inferred_tree_file in snakemake.input['inferredTrees']:
            if dataset_name in inferred_tree_file:
                inferred_trees_file_names.append(
                    os.path.abspath(inferred_tree_file))

    tree_collection = pd.DataFrame({'cell_num': snakemake.params['cellNum'],
                                    'coverage_mean': snakemake.params["covMean"],
                                    'coverage_variance': snakemake.params["covVariance"],
                                    'eff_seq_err_rate': snakemake.params['effSeqErrRate'],
                                    'ado_rate': snakemake.params['adoRate'],
                                    'deletion_rate': snakemake.params['deletionRate'],
                                    'insertion_rate': snakemake.params['insertionRate'],
                                    'dataset': pd.Series(dataset),
                                    'tool': snakemake.params['tool'],
                                    'snv_type': snakemake.params['snvType'],
                                    'tool_setup': snakemake.params['toolSetup'], 
                                    'inferred_tree': pd.Series(inferred_trees_file_names),
                                    'real_tree': pd.Series(true_trees_file_names)
                                    })
else:
    tool_setup = []
    for sieve_run_template in sorted(snakemake.params['sieveRunTemplates']):
        for dataset_name in sorted(snakemake.params['datasetNames']):
            tool_setup.append(sieve_run_template)
            dataset.append(dataset_name)
            true_trees_file_names.append(os.path.abspath(
                snakemake.params['trueTreesDir'] + snakemake.params['trueTreesPre'] + '.' + dataset_name))
            for inferred_tree_file in snakemake.input['inferredTrees']:
                if sieve_run_template in inferred_tree_file and dataset_name in inferred_tree_file:
                    inferred_trees_file_names.append(
                        os.path.abspath(inferred_tree_file))
                    break

    tree_collection = pd.DataFrame({'cell_num': snakemake.params['cellNum'],
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
                                    'inferred_tree': pd.Series(inferred_trees_file_names),
                                    'real_tree': pd.Series(true_trees_file_names)
                                    })

tree_collection.to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
