configfile: "data.yaml"
configfile: "tools.yaml"

import os
import utils

output_dir = config.get(utils.OUTPUT_DIR, "")
alignment_output_dir = os.path.join(output_dir, utils.ALIGNMENTS)
svs_output_dir = os.path.join(output_dir, utils.SVS)
raw_svs_output_dir = os.path.join(svs_output_dir, utils.RAW)

utils.ensure_samples_correctness(config)
sample_to_reads_paths = utils.get_samples_to_reads_paths(config)
utils.ensure_ref_correctness(config)

sniffles_sens_suffix = utils.get_sniffles_sens_suffix(config)
samples_regex = utils.get_samples_regex(sample_to_reads_paths)

sniffles_config = config.get(utils.TOOLS, {}).get(utils.SNIFFLES, {})
jasmine_config=config.get(utils.TOOLS, {}).get(utils.JASMINE, {})
iris_config=config.get(utils.TOOLS, {}).get(utils.IRIS, {})
tech_regex = utils.get_tech_regex(config)
java_config=config.get(utils.TOOLS, {}).get(utils.JAVA, {})

rule get_raw_specific:
    output: protected(os.path.join(raw_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".specific.vcf"))
    input: os.path.join(raw_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".specific_marked.vcf")
    run:
       shell('grep "#" {input[0]} > {output[0]}')
       shell('grep "IN_SPECIFIC=1" {input[0]} >> {output[0]}')

rule mark_specific_in_raw:
    output: temp(os.path.join(raw_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".specific_marked.vcf"))
    input: os.path.join(raw_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".vcf")
    threads: lambda wc: min(cluster_config.get("sensitive_ins_to_dup_conversion", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), jasmine_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(raw_svs_output_dir, utils.LOG, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".specific_marked.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 4000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        ins_to_dup_script_name=jasmine_config.get(utils.SCRIPT_NAME, "InsertionsToDuplications")
    shell:
        "{params.java} -cp {params.java_src} {params.ins_to_dup_script_name} {input} {output} &> {log}"

rule sensitive_svs_sniffles:
    input: os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam")
    output: protected(os.path.join(raw_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".vcf"))
    threads: lambda wildcards: min(cluster_config.get("sensitive_svs_sniffles", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), sniffles_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
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
        min_seq_size = sniffles_config.get(utils.MIN_SEQ_SIZE, 1000)
    shell:
        "{params.sniffles} -m {input} -v {output} --threads {threads} --min_support {params.min_support} --max_distance {params.max_distance} --max_num_splits {params.max_num_splits} --min_length {params.min_length} --num_reads_report {params.num_reads_report} --min_seq_size {params.min_seq_size} &> {log}"

include: "align_single.snakefile"