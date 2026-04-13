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

def get_passed_samples(wildcards):
    with open("passed_file.csv","r") as f:
        samples = [line.strip() for line in f]
    return samples
    
