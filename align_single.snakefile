configfile: "data.yaml"
configfile: "tools.yaml"

import os
import utils
from collections import defaultdict

output_dir = config.get(utils.OUTPUT_DIR, ".")
alignment_output_dir = os.path.join(output_dir, utils.ALIGNMENTS)

utils.ensure_samples_correctness(config)
samples_to_reads_paths = utils.get_samples_to_reads_paths(config)
samples_to_basename_readpaths = defaultdict(dict)
for sample_name, read_paths in samples_to_reads_paths.items():
    for read_path in read_paths:
        samples_to_basename_readpaths[sample_name][os.path.basename(read_path)] = read_path
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
    input:
        fasta=dynamic(os.path.join(alignment_output_dir, "{sample}_{tech}_fasta_{chunk_id}.sort.bam")),
        fastq=dynamic(os.path.join(alignment_output_dir, "{sample}_{tech}_fastq_{chunk_id}.sort.bam")),
    output: protected(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech}.sort.bam"))
    message: "Combining sorted bam files. Requested mem {resources.mem_mb}M."
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.sort.bam.log")
    resources:
        mem_mb = lambda wildcards: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000)
    params:
        samtools = samtools_config.get(utils.PATH, "samtools"),
    shell:
         "{params.samtools} merge -o {output} {input} &> {log}"

rule single_sam_to_sort_bam:
    input: os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}_{chunk_id}.sam")
    output: temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{read_base," + read_paths_regex + "}.sort.bam"))
    threads: lambda wildcards: min(cluster_config.get("single_sam_to_sort_bam", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), samtools_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Transforming an alignment sam file {input} into a sorted bam file {output}. Requested mem {resources.mem_mb}M on {threads} threads. Cluster config "
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{read_base}.sort.bam.log")
    resources:
        mem_mb = lambda wildcards, threads: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000) *  threads
    params:
        tmp_dir = lambda wc: os.path.join(config["tools"].get(utils.TMP_DIR, "."), f"samtools_tmp_{wc.sample}"),
        samtools = samtools_config.get(utils.PATH, "samtools"),
        mem_mb_per_thread = samtools_config.get(utils.MEM_MB_PER_THREAD, 1000),
    shell:
         "{params.samtools} sort -O bam -o {output} -@ {threads} -m {params.mem_mb_per_thread}M -T {params.tmp_dir} {input} &> {log}"

rule single_alignment:
    output: temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format}_{chunk_id}.sam"))
    input: os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}_{chunk_id}.{seq_format}")
    threads: lambda wildcards: min(cluster_config.get("single_alignment", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), ngmlr_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Aligning reads from {input} with NGMLR to {output}. Requested mem {resources.mem_mb}M on {threads} threads. Cluster config "
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{read_base}.sam.log")
    resources:
        mem_mb = lambda wildcards, threads: ngmlr_config.get(utils.MEM_MB_CORE, 5000) + ngmlr_config.get(utils.MEM_MB_PER_THREAD, 500) * threads
    params:
        ngmlr = ngmlr_config.get(utils.PATH, "ngmlr"),
        tech_config = lambda wildcards: "ont" if wildcards.tech.lower() == "ont" else "pacbio",
        reference = config[utils.REFERENCE],
    shell:
         "{params.ngmlr} -r {params.reference} -q {input} -t {threads} -o {output} --bam-fix -x {params.tech_config} &> {log}"

rule ensure_ngmlr_input_extension:
    input: os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format}_{chunk_id}")
    output: temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format}_{chunk_id}.{seq_format}"))
    shell: "mv {input} {output}"

rule split_fastq:
    output:
          temp(dynamic(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fastq_{chunk_id}")))
    input:
        fastq=[],
        fastq_gz=[],
    params:
        prefix=os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format}_"),
    shell:
        "cut <({params.cut} {input.fastq}) & <({params.zcat} {input.fastq_gz}) | split -l {params.fastq_cnt} -a 3 - {params.prefix}"

rule split_fasta:
    output:
          temp(dynamic(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fasta_{chunk_id}")))
    input:
        fasta=[],
        fasta_gz=[],
    params:
        prefix=os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format}_"),
    shell:
        "cut <({params.cut} {input.fasta}) & <({params.zcat} {input.fasta_gz}) | split -l {params.fasta_cnt} -a 3 - {params.prefix}"

localrules: ensure_ngmlr_input_extension
