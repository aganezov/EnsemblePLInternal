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
        samtools=samtools_config.get(utils.PATH, "samtools"),
        awk=awk_config.get(utils.PATH, "awk")
    shell:
        "{params.samtools} depth -a {input} | {params.awk} \'{{sum += $3}} END {{print \"Average coverage (all) = \",sum/NR}}\' > {output} 2> {log}"


def read_extensions_per_sample(sample):
    result = set()
    for read_path in samples_to_reads_paths[sample]:
        if read_path.endswith(("fastq", "fq", "fastq.gz", "fq.gz")):
            result.add("fastq")
        elif read_path.endswith(("fasta", "fa", "fasta.gz", "fa.gz")):
            result.add("fasta")
    return sorted(result)


def aggregated_input_for_bam_merging(wildcards):
    extensions = read_extensions_per_sample(sample=wildcards.sample)
    result = []
    if "fasta" in extensions:
        chekpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
        result.extend(expand(
            os.path.join(alignment_output_dir, f"{wildcards.sample}_{wildcards.tech}_fasta_" + "{chunk_id}.sort.bam"),
            chunk_id=glob_wildcards(os.path.join(chekpoint_output, f"{wildcards.sample}_{wildcards.tech}_fasta_" + "_{chunk_id}")).chunk_id
        ))
    if "fastq" in extensions:
        chekpoint_output = checkpoints.split_fastq.get(**wildcards).output[0]
        result.extend(expand(
            os.path.join(alignment_output_dir, f"{wildcards.sample}_{wildcards.tech}_fastq_" + "{chunk_id}.sort.bam"),
            chunk_id=glob_wildcards(os.path.join(chekpoint_output, f"{wildcards.sample}_{wildcards.tech}_fastq_" + "{chunk_id}")).chunk_id
        ))
    return result

def split_fastx_dirs(wildcards):
    extensions = read_extensions_per_sample(sample=wildcards.sample)
    result = []
    if "fasta" in extensions:
        result.append(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fasta"))
    if "fastq" in extensions:
        result.append(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fastq"))
    return result


rule merge_sorted:
    output: protected(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.sort.bam"))
    input: bams=aggregated_input_for_bam_merging,
           tmp_fastq_dir=split_fastx_dirs,
    message: "Combining sorted bam files. Requested mem {resources.mem_mb}M."
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.sort.bam.log")
    resources:
        mem_mb=lambda wildcards: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000)
    params:
        samtools=samtools_config.get(utils.PATH, "samtools"),
    shell:
        "{params.samtools} merge {output} {input.bams} &> {log}"

rule single_sam_to_sort_bam:
    output:temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format,(fastq|fasta)}_{chunk_id,[a-z]+}.sort.bam"))
    input:
        sam=os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}_{chunk_id}.sam"),
        tmp_dir=lambda wc: os.path.join(config["tools"].get(utils.TMP_DIR, ""), f"samtools_tmp_{wc.sample}_{wc.tech}_{wc.seq_format}_{wc.chunk_id}")
    threads:lambda wildcards: min(cluster_config.get("single_sam_to_sort_bam", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), samtools_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Transforming an alignment sam file {input} into a sorted bam file {output}. Requested mem {resources.mem_mb}M on {threads} threads. Cluster config "
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{seq_format}", "{sample}_{tech}_{seq_format}_{chunk_id}.sort.bam.log")
    resources:
        mem_mb=lambda wildcards, threads: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        tmp_dir=lambda wc: os.path.join(config["tools"].get(utils.TMP_DIR, ""), f"samtools_tmp_{wc.sample}_{wc.tech}_{wc.seq_format}_{wc.chunk_id}" + os.path.sep),
        samtools=samtools_config.get(utils.PATH, "samtools"),
        mem_mb_per_thread=samtools_config.get(utils.MEM_MB_PER_THREAD, 1000),
    shell:
        "{params.samtools} sort -O bam -o {output} -@ {threads} -m {params.mem_mb_per_thread}M -T {params.tmp_dir} {input.sam} &> {log}"

rule single_alignment:
    output:temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format,(fastq|fasta)}" + "_{chunk_id,[a-z]+}.sam"))
    input: os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}_{chunk_id}.{seq_format}")
    threads: lambda wildcards: min(cluster_config.get("single_alignment", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), ngmlr_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Aligning reads from {input} with NGMLR to {output}. Requested mem {resources.mem_mb}M on {threads} threads. Cluster config "
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{seq_format}", "{sample}_{tech}_{seq_format}_{chunk_id}.sam.log")
    resources:
        mem_mb=lambda wildcards, threads: ngmlr_config.get(utils.MEM_MB_CORE, 5000) + ngmlr_config.get(utils.MEM_MB_PER_THREAD, 500) * threads,
    params:
        ngmlr=ngmlr_config.get(utils.PATH, "ngmlr"),
        tech_config=lambda wildcards: "ont" if wildcards.tech.lower() == "ont" else "pacbio",
        reference=config[utils.REFERENCE],
    shell:
        "{params.ngmlr} -r {params.reference} -q {input} -t {threads} -o {output} --bam-fix -x {params.tech_config} &> {log}"

rule ensure_ngmlr_input_extension:
    input: os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}", "{sample}_{tech}_{seq_format}_{chunk_id}")
    output: temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format,(fastq|fasta)}_{chunk_id,[a-z]+}.{seq_format}"))
    shell: "mv {input} {output} && touch {input}"


def get_fastx_files(sample, extension):
    result = []
    for read_path in samples_to_reads_paths[sample]:
        if read_path.endswith(extension):
            result.append(read_path)
    return result


checkpoint split_fastq:
    output:
        temp(directory(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fastq")))
    input:
        fastq=lambda wc: get_fastx_files(sample=wc.sample, extension=("fastq", "fq")),
        fastq_gz=lambda wc: get_fastx_files(sample=wc.sample, extension=("fastq.gz", "fa.gz")),
    params:
        cut_command=lambda wc: "" if len(get_fastx_files(sample=wc.sample, extension=("fastq", "fq"))) == 0 else f"<(cat {' '.join(get_fastx_files(sample=wc.sample, extension=('fastq', 'fq')))})",
        zcat_command=lambda wc: "" if len(get_fastx_files(sample=wc.sample, extension=("fastq.gz", "fq.gz"))) == 0 else f"<(zcat {' '.join(get_fastx_files(sample=wc.sample, extension=('fastq.gz', 'fq.gz')))})",
        prefix=lambda wc: os.path.join(alignment_output_dir, f"{wc.sample}_{wc.tech}_fastq", f"{wc.sample}_{wc.tech}_fastq_"),
        fastq_cnt=lambda wc: config.get(utils.READS_CNT_PER_RUN, 500000) * 4,
    shell:
        "mkdir -p {output} && cat {params.cut_command} {params.zcat_command} | split -l {params.fastq_cnt} -a 3 - {params.prefix}"

checkpoint split_fasta:
    output:
        temp(directory(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fasta")))
    input:
        fasta=lambda wc: get_fastx_files(sample=wc.sample, extension=("fasta", "fa")),
        fasta_gz=lambda wc: get_fastx_files(sample=wc.sample, extension=("fasta.gz", "fa.gz")),
    params:
        cut_command=lambda wc: "" if len(get_fastx_files(sample=wc.sample, extension=("fasta", "fa"))) == 0 else "<(cat {fasta})".format(fasta=" ".join(get_fastx_files(sample=wc.sample, extension=("fasta", "fa")))),
        zcat_command=lambda wc: "" if len(get_fastx_files(sample=wc.sample, extension=("fasta.gz", "fa.gz"))) == 0 else "<(zcat {fasta_gz})".format(fasta_gz=" ".join(get_fastx_files(sample=wc.sample, extension=("fasta.gz", "fa.gz")))),
        prefix=lambda wc: os.path.join(alignment_output_dir, f"{wc.sample}_{wc.tech}_fasta", f"{wc.sample}_{wc.tech}_fasta_"),
        fasta_cnt=lambda wc: config.get(utils.READS_CNT_PER_RUN, 500000) * 2,
    shell:
        "mkdir -p {output} && cat {params.cut_command} {params.zcat_command} | split -l {params.fasta_cnt} -a 3 - {params.prefix}"

rule samtools_tmp_dir:
    output: temp(directory(os.path.join(config["tools"].get(utils.TMP_DIR, ""), "samtools_tmp_{sample}_{tech}_{seq_format}_{chunk_id}")))
    shell: "mkdir -p {output}"

localrules: ensure_ngmlr_input_extension, samtools_tmp_dir
