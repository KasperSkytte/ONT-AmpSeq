import glob
import os

# Helper function to list all fastq files per wildcard (subfolder/sample)
def listFastq(wildcards):
    fastqs = glob.glob(os.path.join(config['input_dir'], wildcards.sample, "*.fastq*"))
    return fastqs

rule concatenate_fastq:
    input:
        listFastq
    output:
        temp(os.path.join(config['tmp_dir'], "samples", "{sample}_concat.fastq"))
    resources:
        mem_mb = 512,
        runtime = "01:00:00"
    threads: 1
    log:
        os.path.join(config["log_dir"], "concatenate_fastq", "{sample}.log")
    shell:
        """
        # Check if input files are compressed
        if [[ $(file -b --mime-type {input}) == "application/gzip" ]]; then
            zcat {input} > {output}
        else
            cat {input} > {output}
        fi
        """