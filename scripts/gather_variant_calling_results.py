import pandas as pd


variants_info_per_tool = []

for item in snakemake.input:
    variants_info_per_tool.append(pd.read_csv(item, sep='\t', dtype={'dataset': str}))

pd.concat(variants_info_per_tool).to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
