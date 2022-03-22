import pandas as pd


trees_info_per_tool = []

for item in snakemake.input:
    trees_info_per_tool.append(pd.read_csv(
        item, sep='\t', dtype={'dataset': str}))

pd.concat(trees_info_per_tool).to_csv(
    snakemake.output[0], sep='\t', na_rep='NA', index=False)
