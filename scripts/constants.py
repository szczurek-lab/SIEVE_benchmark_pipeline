# This file will be loaded by snakemake files.


configfile: "config.yaml"


##############################################
#         Set up constant variables          #
##############################################

# Constants for simulated data; 'SIMU' for 'SIMULATION'
SIMUDATAROOT = config["benchmark"]["baseDir"] + config["benchmark"]["simulation"]["dir"]
SIMUDATADIR = SIMUDATAROOT + config["benchmark"]["simulation"]["dataDir"]
SIMUTUMORCELLNAMES = SIMUDATADIR + config["benchmark"]["simulation"]["tumorCellNames"]
SIMUSNVSITESDIR = SIMUDATADIR + config["benchmark"]["simulation"]["trueSNVSitesNamesDir"]
SIMUTRUESNVDIR = SIMUDATADIR + config["benchmark"]["simulation"]["trueSNVDir"]
SIMUPARALLELMUTDIR = SIMUDATADIR + config["benchmark"]["simulation"]["parallelMutDir"]
# SIMUTRUESNVADODIR = SIMUDATADIR + config["benchmark"]["simulation"]["trueSNVAdoDir"]
SIMUALLELICSEQINFODIR = SIMUDATADIR + config["benchmark"]["simulation"]["trueAllelicSeqInfoSNVsDir"]
SIMUTRUEADOSTATESDIR = SIMUDATADIR + config["benchmark"]["simulation"]["trueAdoStatesDir"]
SIMUTRUESIZEFACTORSDIR = SIMUDATADIR + config["benchmark"]["simulation"]["trueSizeFactorsDir"]

# Simulated true trees
SIMUTRUETREEDIR = SIMUDATADIR + config["benchmark"]["simulation"]["treesDir"]
SIMUCONVERTEDTRUETREEDIR = SIMUDATADIR + config["benchmark"]["simulation"]["convertedTreesDir"]

# Simulated raw data
SIMURAWDATADIR = SIMUDATADIR + config["benchmark"]["simulation"]["rawDataDir"]

# Simulated aligned reads (MPILEUP)
SIMUMPILEUPDATADIR = SIMUDATADIR + config["benchmark"]["simulation"]["MPileUpDataDir"]
SIMUMPILEUPWITHNORMAL = config["benchmark"]["simulation"]["MPileUpDataWithNormalPre"] + \
    config["benchmark"]["simulation"]["MPileUpDataSuf"]
SIMUMPILEUPWITHOUTNORMAL = config["benchmark"]["simulation"]["MPileUpDataWithoutNormalPre"] + \
    config["benchmark"]["simulation"]["MPileUpDataSuf"]

# BAM file names as input for MONOVAR
SIMUMONOVARDIR = SIMUDATADIR + config["benchmark"]["simulation"]["monovarDir"]

# BAM file names as input for SCIPHI
SIMUSCIPHIDIR = SIMUDATADIR + config["benchmark"]["simulation"]["sciphiDir"]

# Results of SIEVE
SIEVEDIR = config["benchmark"]["baseDir"] + \
    config["benchmark"]["sieve"]["dir"]
SIEVEFROMCANDIDATEDIR = SIEVEDIR + config["benchmark"]["phyloInf"]["fromCandidate"]
SIEVEFROMCANDIDATEPHYINFDIR = SIEVEFROMCANDIDATEDIR + config["benchmark"]["sieve"]["phyloInfDir"]
SIEVEFROMCANDIDATEPHYINFSTAGE2DIR = SIEVEFROMCANDIDATEPHYINFDIR + config["benchmark"]["sieve"]["stage2Dir"]
SIEVEFROMCANDIDATEVARCALDIR = SIEVEFROMCANDIDATEDIR + config["benchmark"]["sieve"]["variantCallingDir"]
SIEVEFROMCANDIDATEVARCALSTAGE2DIR = SIEVEFROMCANDIDATEVARCALDIR + config["benchmark"]["sieve"]["stage2Dir"]

# Results of CELLPHY
CELLPHYDIR = config["benchmark"]["baseDir"] + \
    config["benchmark"]["cellphy"]["dir"]
CELLPHYTRUEMONOVARSNVSDIR = CELLPHYDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]

# Git manager
GITSKELETONLOCALDIR = config["benchmark"]["baseDir"] + \
    config["benchmark"]["git"]["skeletonDir"]
GITREMOTEDIR = config["servers"]["remoteServer"]["rootPath"] + \
    config["benchmark"]["git"]["skeletonDir"]

# Results of SIFIT
# Skeleton
SIFITGITSKELETONLOCALDIR = GITSKELETONLOCALDIR + config["benchmark"]["sifit"]["dir"]
SIFITGITSKELETONTRUEMONOVARLOCALDIR = SIFITGITSKELETONLOCALDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]

SIFITGITREMOTEDIR = GITREMOTEDIR + config["benchmark"]["sifit"]["dir"]
SIFITGITSKELETONTRUEMONOVARREMOTEDIR = SIFITGITREMOTEDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]

# Local
SIFITLOCALDIR = config["benchmark"]["baseDir"] + \
    config["benchmark"]["sifit"]["dir"]

SIFITTRUEMONOVARSNVSLOCALDIR = SIFITLOCALDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]
SIFITLOCALFLAGFILETRUEMONOVARINITSKELETONIN = SIFITTRUEMONOVARSNVSLOCALDIR + config["benchmark"]["sifit"]["flag"]["dir"] + \
    config["benchmark"]["sifit"]["flag"]["initSkeletonIn"]
SIFITLOCALFLAGFILETRUEMONOVARINITSKELETONOUT = SIFITTRUEMONOVARSNVSLOCALDIR + config["benchmark"]["sifit"]["flag"]["dir"] + \
    config["benchmark"]["sifit"]["flag"]["initSkeletonOut"]
SIFITLOCALFLAGFILEPUSHEDTRUEMONOVARDATA = SIFITTRUEMONOVARSNVSLOCALDIR + config["benchmark"]["sifit"]["flag"]["dir"] + \
    config["benchmark"]["sifit"]["flag"]["pushedData"]
SIFITTRUEMONOVARSNVSTRUEPARAMSLOCALDIR = SIFITTRUEMONOVARSNVSLOCALDIR + config["benchmark"]["paramTypes"]["fromTrueParamsDir"]

# Remote, if applicable
SIFITREMOTEDIR = config["servers"]["remoteServer"]["rootPath"] + \
    config["benchmark"]["sifit"]["dir"]
SIFITTRUEMONOVARSNVSREMOTEDIR = SIFITREMOTEDIR + config["benchmark"]["phyloInf"]["fromTrueMonovarDir"]

# Results of MONOVAR
MONOVARDIR = config["benchmark"]["baseDir"] + \
    config["benchmark"]["monovar"]["dir"]
MONOVARTRUEPARAMSDIR = MONOVARDIR + config["benchmark"]["paramTypes"]["fromTrueParamsDir"]

# Results of SCIPHI
SCIPHIDIR = config["benchmark"]["baseDir"] + \
    config["benchmark"]["sciphi"]["dir"]

# Analysis results
ANALYSISDIR = config["benchmark"]["baseDir"] + config["benchmark"]["analysis"]["dir"]

# Analysis results of tree comparison
ANALYSISTREECOMPDIR = ANALYSISDIR + config["benchmark"]["analysis"]["treeComparisonDir"]

# Analysis results of parameter estimates
ANALYSISPARAMSDIR = ANALYSISDIR + config["benchmark"]["analysis"]["paramsEstimatesDir"]

# Analysis results of variant calling
ANALYSISVARCALLDIR = ANALYSISDIR + config["benchmark"]["analysis"]["varCallResultsDir"]

# Analysis results of ado calling
ANALYSISADOCALLDIR = ANALYSISDIR + config["benchmark"]["analysis"]["adoCallResultsDir"]

# Analysis results of allelic information
ANALYSISALLELICINFO = ANALYSISDIR + config["benchmark"]["analysis"]["allelicInfoPerCellDir"]

# Efficiency benchmarking results
#EFFICIENCYDIR = config["benchmark"]["baseDir"] + config["benchmark"]["efficiency"]["dir"]

#EFFICIENCYMONOVARDIR = EFFICIENCYDIR + config["benchmark"]["monovar"]["dir"]

#EFFICIENCYSCIPHIDIR = EFFICIENCYDIR + config["benchmark"]["sciphi"]["dir"]

#EFFICIENCYSIEVEDIR = EFFICIENCYDIR + config["benchmark"]["sieve"]["dir"]

#EFFICIENCYCELLPHYDIR = EFFICIENCYDIR + config["benchmark"]["cellphy"]["dir"]

#EFFICIENCYSIFITDIR = EFFICIENCYDIR + config["benchmark"]["sifit"]["dir"]
