#!/bin/bash

# USAGE: source run/run.sh path/to/snakefile


##################################
#           USER SETUP           #
##################################

# conda environment name
ENV_NAME="snake"

# whether to use taskset (1) or not (0)
USETASKSET=1

# if use taskset, set starting and ending cpu indices, both included
STARTCPU=0
ENDCPU=7

# if not use taskset, set the number of cores you want to use
CORENUM=1

# total memory limitation in MB
MEMORY=100000

# other flags of snakemake to attach
ATTACHMENTS=-kp

# where to redirect screen log
LOGFILE="./out"


###################################
#           DO NOT EDIT           #
###################################

if [[ -z $1 ]]; then
  echo "> Unset snakemake file; use \"Snakefile\" by default." >>"${LOGFILE}"
  SNAKEFILE=Snakefile
else
  SNAKEFILE=$1
fi

conda activate "${ENV_NAME}"

# the start of analysis
SECONDS=0

if [[ ${USETASKSET} -eq 1 ]]; then
  echo "> Using ${STARTCPU} to ${ENDCPU} cpus..." >>"${LOGFILE}"
  ((CORENUM = ENDCPU - STARTCPU + 1))

  taskset -c ${STARTCPU}-${ENDCPU} nohup snakemake -s "${SNAKEFILE}" --use-conda --resources mem_mb="${MEMORY}" --rerun-triggers mtime --cores "${CORENUM}" "${ATTACHMENTS}" &>"${LOGFILE}" &
elif [[ ${USETASKSET} -eq 0 ]]; then
  nohup snakemake -s "${SNAKEFILE}" --use-conda --resources mem_mb="${MEMORY}" --rerun-triggers mtime --cores "${CORENUM}" "${ATTACHMENTS}" &>"${LOGFILE}" &
else
  echo "> Set USETASKSET to 1 or 0 in order to use taskset or not!"
fi

# the end of analysis
DURATION=$SECONDS

echo "> Time elapsed: $((DURATION / 3600)) hours, $(((DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds" >>"${LOGFILE}"

conda deactivate
