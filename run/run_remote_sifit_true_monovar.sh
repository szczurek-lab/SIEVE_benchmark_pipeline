#!/bin/bash

# USAGE: source run/run_remote_sifit_true_monovar.sh


##################################
#           USER SETUP           #
##################################

# conda environment name
ENV_NAME="snake"

# whether to use taskset (1) or not (0)
USETASKSET=1

# if use taskset, set starting and ending cpu indices, both included
STARTCPU=0
ENDCPU=10

# if not use taskset, set the number of cores you want to use
CORENUM=1

# total memory limitation in MB
MEMORY=200000

# other flags of snakemake to attach
ATTACHMENTS=-kp

# where to redirect screen log
LOGFILE="./out"


###################################
#           DO NOT EDIT           #
###################################

conda activate "${ENV_NAME}"

# the start of analysis
SECONDS=0

if [[ ${USETASKSET} -eq 1 ]]; then
  echo "Using ${STARTCPU} to ${ENDCPU} cpus..." >"${LOGFILE}"
  ((CORENUM = ENDCPU - STARTCPU + 1))
  taskset -c ${STARTCPU}-${ENDCPU} snakemake -s remote_sifit_true_monovar.snake --use-conda --resources mem_mb="${MEMORY}" --rerun-trigger mtime --cores "${CORENUM}" "${ATTACHMENTS}" &>"${LOGFILE}"
elif [[ ${USETASKSET} -eq 0 ]]; then
  snakemake -s remote_sifit_true_monovar.snake --use-conda --resources mem_mb="${MEMORY}" --rerun-trigger mtime --cores "${CORENUM}" "${ATTACHMENTS}" &>"${LOGFILE}"
else
  echo "Set USETASKSET to 1 or 0 in order to use taskset or not!"
  exit 1
fi

# the end of analysis
DURATION=$SECONDS

echo "Time elapsed: $((DURATION / 3600)) hours, $(((DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds" >>"${LOGFILE}"

conda deactivate

exit 0
