configfile: "config.yaml"
include: "scripts/constants.py"
include: "scripts/utils.py"


#########################################################
#                Generate simulated data                #
#########################################################

generate_simulated_data(
    SIMUDATAROOT,
    config["tools"]["SIEVE_simulator"],
    config["configFiles"]["SIEVE_simulator"],
    SIMUDATAROOT + "simulation.log",
    config["benchmark"]["simulation"]["dataDir"],
    SIMUMONOVARDIR + config["benchmark"]["simulation"]["bamFileNamesWithNormal"],
    SIMUMONOVARDIR + config["benchmark"]["simulation"]["bamFileNamesWithoutNormal"],
    SIMUSCIPHIDIR + config["benchmark"]["simulation"]["bamFileNamesWithNormal"],
    SIMUSCIPHIDIR + config["benchmark"]["simulation"]["bamFileNamesWithoutNormal"]
)


########################################################
#          Include predefined snakemake files          #
########################################################

include: "simulated_data.snake"
include: "sieve_2stage.snake"
include: "sciphi.snake"
include: "monovar.snake"

# Monovar must be included for CellPhy.
include: "cellphy.snake"

# Monovar must be included for SiFit.
# Three options for SiFit is available. Choose one, comment out the other two.
# 1. Run SiFit on the local machine.
include: "sifit_local.snake"
# 2. Results from SiFit are available. Only collect results without running SiFit locally. This mostly happens when SiFit has successfully been run.
# include: "sifit_fixed.snake"
# 3. Run SiFit on a remote server.
# include: "sifit_remote.snake"


#########################################################
#                       All rules                       #
#########################################################

rule all:
    input:
        # Inferred trees
        ANALYSISTREECOMPDIR + "trees_info_updated.tsv",
        
        # Parameter estimates
        ANALYSISPARAMSDIR + "params_info.tsv",

        # Variant calling results
        ANALYSISVARCALLDIR + "variants_info.tsv",
        
        # ADO calling results
        ANALYSISADOCALLDIR + "ado_info.tsv",

        # Allelic information
        ANALYSISALLELICINFO + "allelic_info.rds",

        # Sites information
        ANALYSISSITESINFODIR + "sites_info.tsv",

        # Summary for parallel mutations
        SIMUPARALLELMUTDIR + "summary.tsv"


#######################################################
#                                                     #
#                Gather inferred trees                #
#                                                     #
#######################################################

# Rule: gather all inferred trees from different tools
rule gatherInferredTrees:
    input:
        ANALYSISTREECOMPDIR + "from_sciphi.tsv",
        ANALYSISTREECOMPDIR + "from_sieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISTREECOMPDIR + "from_sieve_stage2_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISTREECOMPDIR + "from_cellphy_" + config["benchmark"]["analysis"]["trueMonovarSNVs"] + "_" + config["benchmark"]["analysis"]["cellphyEPMode"] + ".tsv",
        ANALYSISTREECOMPDIR + "from_sifit_" + config["benchmark"]["analysis"]["trueMonovarSNVs"] + "_" + config["benchmark"]["analysis"]["trueParameters"] + ".tsv"
    output:
        protected(ANALYSISTREECOMPDIR + "trees_info.tsv")
    script:
        "scripts/gather_inferred_trees.py"

# Rule: get tree comparison results
rule getTreeComparisonResults:
    input:
        ANALYSISTREECOMPDIR + "trees_info.tsv"
    output:
        treeInfo=protected(ANALYSISTREECOMPDIR + "trees_info_updated.tsv")
    params:
        log=ANALYSISTREECOMPDIR + "out"
    shell:
        "Rscript --vanilla scripts/compare_trees.R "
        "-i {input} "
        "-o {output.treeInfo} &>{params.log}"


#######################################################
#                                                     #
#               Get parameter estimates               #
#                                                     #
#######################################################

# Rule: get all parameter estimates from different tools
rule getParams:
    input:
        ANALYSISPARAMSDIR + "from_sciphi.tsv",
        ANALYSISPARAMSDIR + "from_sieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISPARAMSDIR + "from_sieve_stage2_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISPARAMSDIR + "from_cellphy_" + config["benchmark"]["analysis"]["trueMonovarSNVs"] + ".tsv",
        ANALYSISPARAMSDIR + "from_sifit_" + config["benchmark"]["analysis"]["trueMonovarSNVs"] + "_" + config["benchmark"]["analysis"]["trueParameters"] + ".tsv"
    output:
        protected(ANALYSISPARAMSDIR + "params_info.tsv")
    script:
        "scripts/get_param_estimates.py"


#######################################################
#                                                     #
#             Get variant calling results             #
#                                                     #
#######################################################

# Rule: gather variant calling results from different tools
rule gatherVarResults:
    input: 
        ANALYSISVARCALLDIR + "from_monovar_" + config["benchmark"]["analysis"]["trueParameters"] + ".tsv",
        ANALYSISVARCALLDIR + "from_sciphi.tsv",
        ANALYSISVARCALLDIR + "from_sieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISVARCALLDIR + "from_sieve_stage2_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv"
    output: 
        protected(ANALYSISVARCALLDIR + "variants_info.tsv")
    script:
        "scripts/gather_variant_calling_results.py"


#######################################################
#                                                     #
#               Get ADO calling results               #
#                                                     #
#######################################################

# Rule: gather ADO calling results from different tools
rule gatherAdoResults:
    input: 
        ANALYSISADOCALLDIR + "from_sieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISADOCALLDIR + "from_sieve_stage2_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv"
    output: 
        protected(ANALYSISADOCALLDIR + "ado_info.tsv")
    script:
        "scripts/gather_ado_calling_results.py"


#######################################################
#                                                     #
#                Get sites information                #
#                                                     #
#######################################################

# Rule: gather sites information from different tools
rule gatherSitesInfo:
    input:
        ANALYSISSITESINFODIR + "from_sciphi.tsv",
        ANALYSISSITESINFODIR + "from_sieve_" + config["benchmark"]["analysis"]["candidateSNVs"] + ".tsv",
        ANALYSISSITESINFODIR + "from_cellphy_" + config["benchmark"]["analysis"]["trueMonovarSNVs"] + ".tsv",
        ANALYSISSITESINFODIR + "from_sifit_" + config["benchmark"]["analysis"]["trueMonovarSNVs"] + "_" + config["benchmark"]["analysis"]["trueParameters"] + ".tsv"
    output:
        protected(ANALYSISSITESINFODIR + "sites_info.tsv")
    script:
        "scripts/gather_sites_info.py"

