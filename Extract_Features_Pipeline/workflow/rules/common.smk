# modules/rules/common.smk
# functions shared by all modules

#function to get illumina raw reads. returns R1 and R2 reads for paired-end experiments
def get_illumina_raw(wildcards):
    # Get sample index
    samples = config["samples"]["name"]
    sample_idx = samples.index(wildcards.sample)
    
    # Get the paths from config
    r1_paths = config["samples"]["raw_reads_r1"]
    r2_paths = config["samples"]["raw_reads_r2"]
    
    return {
        "r1": r1_paths[sample_idx],
        "r2": r2_paths[sample_idx]
        }

def get_passed_samples():
    # don't run until checkpoint has been done
    checkpoints.filter.get()

    with open("results/00_QC/passed_samples.csv", "r") as f:
        return [line.strip() for line in f if line.strip()]

# Job resources by tier, sized for nodes with 32 CPUs / 128 GB (runtime in minutes):
#   heavy  - bwa mem, Pindel (2 fit on a node)
#   medium - fastp, samtools sort, Picard, MultiQC, mosdepth, SURVIVOR, matrix building
#   light  - small single-threaded steps (indexing, sorting VCFs, copying files)
# Override a tier by adding it under `resources:` in config/config.yaml, e.g.
#   resources:
#     heavy: {runtime: 2880, mem_mb: 32000, cpus: 16}
RESOURCE_TIERS = {
    "heavy":  {"runtime": 1440, "mem_mb": 32000, "cpus": 16},
    "medium": {"runtime": 360,  "mem_mb": 8000,  "cpus": 4},
    "light":  {"runtime": 120,  "mem_mb": 4000,  "cpus": 1},
}

def res(tier, key):
    return config.get("resources", {}).get(tier, RESOURCE_TIERS[tier])[key]
