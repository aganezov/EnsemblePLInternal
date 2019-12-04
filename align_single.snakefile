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

samples_regex = utils.get_samples_regex(samples_to_reads_paths)
read_paths_regex = utils.get_reads_paths_regex(samples_to_reads_paths)
tech_regex = utils.get_tech_regex(config)

rule merged_coverage:
    input: os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam")
    output: os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.coverage.txt")
    message: "Computing average alignment read depth coverage on {input}"
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.coverage.txt.log")
    params:
        samtools = samtools_config.get(utils.PATH, "samtools"),
        awk = awk_config.get(utils.PATH, "awk")
    shell:
        "{params.samtools} depth -a {input} | {params.awk} \'{{sum += $3}} END {{print \"Average coverage (all) = \",sum/NR}}\' > {output} 2> {log}"

rule merge_sorted:
    input: lambda wildcards: [os.path.join(alignment_output_dir, wildcards.sample + "_" + wildcards.tech + "_" + os.path.basename(read_path) + ".sort.bam") for read_path in samples_to_reads_paths[wildcards.sample]]
    output: protected(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech}.sort.bam"))
    message: "Combining sorted bam files. Requested mem {resources.mem_mb}M."
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.sort.bam.log")
    resources:
        mem_mb = lambda wildcards: samtools_config.get(utils.MEM_MB_PER_THREAD, 1000)
    params:
        samtools = samtools_config.get(utils.PATH, "samtools"),
    shell:
         "{params.samtools} merge -o {output} {input} &> {log}"

rule single_sam_to_sort_bam:
    input: os.path.join(alignment_output_dir, "{sample}_{tech}_{read_base}.sam")
    output: temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{read_base," + read_paths_regex + "}.sort.bam"))
    threads: lambda wildcards: min(cluster_config.get("single_sam_to_sort_bam", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), samtools_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Transforming an alignment sam file {input} into a sorted bam file {output}. Requested mem {resources.mem_mb}M on {threads} threads."
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{read_base}.sort.bam.log")
    resources:
        mem_mb = lambda wildcards, threads: samtools_config.get(utils.MEM_MB_PER_THREAD, 1000) *  threads
    params:
        tmp_dir = lambda wc: os.path.join(config["tools"].get(utils.TMP_DIR, "."), f"samtools_tmp_{wc.sample}"),
        samtools = samtools_config.get(utils.PATH, "samtools"),
        mem_mb_per_thread = samtools_config.get(utils.MEM_MB_PER_THREAD, 1000),
    shell:
         "{params.samtools} sort -O bam -o {output} -@ {threads} -m {params.mem_mb_per_thread}M -T {params.tmp_dir} {input} &> {log}"

rule single_alignment:
    output: temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{read_base," + read_paths_regex + "}.sam"))
    input: lambda wildcards: [read_path for read_path in samples_to_reads_paths[wildcards.sample] if read_path.endswith(wildcards.read_base)][0]
    threads: lambda wildcards: min(cluster_config.get("single_alignment", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), ngmlr_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Aligning reads from {input} with NGMLR to {output}. Requested mem {resources.mem_mb}M on {threads} threads. Cluster config {cluster_config}."
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{read_base}.sam.log")
    resources:
        mem_mb = lambda wildcards, threads: ngmlr_config.get(utils.MEM_MB_CORE, 5000) + ngmlr_config.get(utils.MEM_MB_PER_THREAD, 500) * threads
    params:
        ngmlr = ngmlr_config.get(utils.PATH, "ngmlr"),
        tech_config = lambda wildcards: "ont" if wildcards.tech.lower() == "ont" else "pacbio",
        reference = config[utils.REFERENCE],
    shell:
         "{params.ngmlr} -r {params.reference} -q {input} -t {threads} -o {output} --bam-fix -x {params.tech_config} &> {log}"