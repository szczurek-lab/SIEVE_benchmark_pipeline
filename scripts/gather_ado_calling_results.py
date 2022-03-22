import pandas as pd


ado_info_per_tool = []

for item in snakemake.input:
    ado_info_per_tool.append(pd.read_csv(item, sep='\t', dtype={'dataset': str}))

pd.concat(ado_info_per_tool).to_csv(snakemake.output[0], sep='\t', na_rep='NA', index=False)
