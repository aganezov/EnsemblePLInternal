configfile: "data.yaml"
configfile: "tools.yaml"

import os
import utils

output_dir = config.get(utils.OUTPUT_DIR, ".")
alignment_output_dir = os.path.join(output_dir, utils.ALIGNMENTS)

utils.ensure_samples_correctness(config)
samples_to_reads_paths = utils.get_samples_to_reads_paths(config)
utils.ensure_ref_correctness(config)


ngmlr_config = config.get(utils.TOOLS, {}).get(utils.NGMLR, {})
samtools_config = config.get(utils.TOOLS, {}).get(utils.SAMTOOLS, {})
awk_config = config.get(utils.TOOLS, {}).get(utils.AWK, {})

rule merged_coverage:
    input: os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam")
    output: os.path.join(alignment_output_dir, "{sample}_{tech}.coverage.txt")
    message: "Computing average alignment read depth coverage on {input}"
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.coverage.txt.log")
    params:
        samtools = samtools_config.get(utils.PATH, "samtools"),
        awk = awk_config.get(utils.PATH, "awk")
    shell:
        "{params.samtools} depth -a {input} | {params.awk} \'{{sum += $3}} END {{print \"Average coverage (all) = \",sum/NR}}\' > {output} 2> {log}"

rule merged_sorted:
    input: lambda wildcards: [os.path.join(alignment_output_dir, wildcards.sample + "_" + wildcards.tech + "_" + os.path.basename(read_path) + ".sort.bam") for read_path in samples_to_reads_paths[wildcards.sample]]
    output: protected(os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam"))
    message: f"Combining sorted bam files"
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.sort.bam.log")
    params:
        samtools = samtools_config.get(utils.PATH, "samtools"),
    shell:
         "{params.samtools} merge -o {output} {input} &> {log}"

rule single_sam_to_sort_bam:
    input: os.path.join(alignment_output_dir, "{sample}_{tech}_{read_base}.sam")
    output: temp(os.path.join(alignment_output_dir, "{sample}_{tech}_{read_base}.sort.bam"))
    threads: samtools_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT)
    message: "Transforming an alignment sam file {input} into a sorted bam file {output}"
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{read_base}.sort.bam.log")
    params:
        tmp_dir = lambda wc: os.path.join(config["tools"].get(utils.TMP_DIR, "."), f"samtools_tmp_{wc.sample}"),
        samtools = samtools_config.get(utils.PATH, "samtools"),
    shell:
         "{params.samtools} sort -O bam -o {output} -@ {threads} -T {params.tmp_dir} {input} &> {log}"

rule single_alignment:
    output: temp(os.path.join(alignment_output_dir, "{sample}_{tech}_{read_base}.sam"))
    input: lambda wildcards: [read_path for read_path in samples_to_reads_paths[wildcards.sample] if read_path.endswith(wildcards.read_base)][0]
    threads: ngmlr_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT)
    message: "Aligning reads from {input} with NGMLR to {output}"
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{read_base}.sam.log")
    params:
        ngmlr = ngmlr_config.get(utils.PATH, "ngmlr"),
        tech_config = lambda wildcards: "ont" if wildcards.tech.lower() == "ont" else "pacbio",
        reference = config[utils.REFERENCE],
    shell:
         "{params.ngmlr} -r {params.reference} -q {input} -t {threads} -o {output} --bam-fix -x {params.tech_config} &> {log}"