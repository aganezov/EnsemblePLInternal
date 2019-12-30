configfile: "data.yaml"
configfile: "tools.yaml"

import os
import utils

output_dir = config.get(utils.OUTPUT_DIR, ".")
alignment_output_dir = os.path.join(output_dir, utils.ALIGNMENTS)
svs_output_dir = os.path.join(output_dir, )
raw_svs_output_dir = os.path.join(output_dir, utils.RAW)

utils.ensure_samples_correctness(config)
sample_to_reads_paths = utils.get_samples_to_reads_paths(config)
utils.ensure_ref_correctness(config)

sniffles_sens_suffix = utils.get_sniffles_sens_suffix(config)
samples_regex = utils.get_samples_regex(sample_to_reads_paths)

sniffles_config = config.get(utils.TOOLS, {}).get(utils.SNIFFLES, {})
tech_regex = utils.get_tech_regex(config)

rule sensitive_svs_sniffles:
    input: os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.sort.bam")
    output: protected(os.path.join(raw_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".vcf"))
    threads: lambda wildcards: min(cluster_config.get("sensitive_svs_sniffles", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), sniffles_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Calling sensitive set of SVs on {input} with output stored in {output}. Requested {resources.mem_mb}M of memory on {threads} threads. Cluster config "
    log: os.path.join(raw_svs_output_dir, utils.LOG, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".vcf.log")
    resources:
        mem_mb = lambda wildcards, threads: sniffles_config.get(utils.MEM_MB_CORE, 10000) + sniffles_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        sniffles = sniffles_config.get(utils.PATH, "sniffles"),
        min_length = sniffles_config.get(utils.MIN_LENGTH, 20),
        min_support = sniffles_config.get(utils.MIN_SUPPORT, 2),
        max_num_splits = sniffles_config.get(utils.MAX_NUM_SPLIT_READS, 10),
        max_distance = sniffles_config.get(utils.MAX_DISTANCE, 1000),
        num_reads_report = sniffles_config.get(utils.NUM_READS_REPORT, -1),
        min_seq_size = sniffles_config.get(utils.MIN_SEQ_SIZE, 1000),
        ccs_flag = lambda wc: "--ccs_reads" if wc.tech.lower() in ["pacbioccs", "pbccs"] else "",
    shell:
        "{params.sniffles} -m {input} -v {output} --threads {threads} --min_support {params.min_support} --max_distance {params.max_distance} --max_num_splits {params.max_num_splits} "
        "--min_length {params.min_length} --num_reads_report {params.num_reads_report} --min_seq_size {params.min_seq_size} {params.ccs_flag} &> {log}"

include: "align_single.snakefile"