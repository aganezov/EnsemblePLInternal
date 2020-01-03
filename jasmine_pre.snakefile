configfile: "data.yaml"
configfile: "tools.yaml"

import os
import utils

output_dir = config.get(utils.OUTPUT_DIR, "")
alignment_output_dir = os.path.join(output_dir, utils.ALIGNMENTS)
svs_output_dir = os.path.join(output_dir, utils.SVS)
raw_svs_output_dir = os.path.join(svs_output_dir, utils.RAW)
refined_svs_output_dir = os.path.join(svs_output_dir, utils.REFINED)
ins_to_dup_output_dir = os.path.join(refined_svs_output_dir, utils.INS_TO_DUP)
iris_refined_output_dir = os.path.join(refined_svs_output_dir, utils.IRIS_REFINED)
specific_marked_output_dir = os.path.join(refined_svs_output_dir, utils.SPECIFIC_MARKED)

sniffles_sens_suffix = utils.get_sniffles_sens_suffix(config)
tech_regex = utils.get_tech_regex(config)
samples_regex = utils.get_samples_regex(sample_to_reads_paths)

samtools_config = config.get(utils.TOOLS, {}).get(utils.SAMTOOLS, {})
java_config=config.get(utils.TOOLS, {}).get(utils.JAVA, {})
jasmine_config=config.get(utils.TOOLS, {}).get(utils.JASMINE, {})
iris_config=config.get(utils.TOOLS, {}).get(utils.IRIS, {})
minimap2_config=config.get(utils.TOOLS, {}).get(utils.MINIMAP2, {})
racon_config=config.get(utils.TOOLS, {}).get(utils.RACON, {})


def get_min_support(coverage_file, min_support_fixed_cnt, min_support_fraction):
    coverage = 100
    with open(coverage_file, "rt") as source:
        for line in source:
            coverage = int(float(line.strip().split("=")[1].strip()))
            break
    return min(int(min_support_fixed_cnt), int(coverage * min_support_fraction))



rule specific_or_sv_types:
    input: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.specific.vcf")
    run:
        shell(f'grep "#" {input[0]} > {output[0]}')
        shell(f'grep "IN_SPECIFIC=1" {input[0]} >> {output[0]}')

rule spec_marked_sensitive_or_sv_types_ins_to_dup:
    input: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.vcf")
    threads: lambda wc: min(cluster_config.get("spec_marked_sensitive_or_sv_types_ins_to_dup", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), jasmine_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(refined_svs_output_dir, utils.LOG, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 4000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        ins_to_dup_script_name=jasmine_config.get(utils.SCRIPT_NAME, "InsertionsToDuplications")
    shell:
        "{params.java} -cp {params.java_src} {params.ins_to_dup_script_name} {input} {output} &> {log}"

rule specific_new_sv_types:
    input: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.specific.vcf")
    run:
        shell(f'grep "#" {input[0]} > {output[0]}')
        shell(f'grep "IN_SPECIFIC=1" {input[0]} >> {output[0]}')

rule spec_marked_sensitive_new_sv_type_final_location:
    input: vcf=os.path.join(specific_marked_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined_markedSpec.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.vcf")
    params:
        min_support_fixed=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_READS_FIXED, 10),
        min_support_fraction=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_READS_FRACTION, 0.25),
        min_length=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_LEN, 30),
    run:
        shell("mv {input} {output}")

rule spec_marked_sensitive_new_sv_types:
    input: vcf=os.path.join(iris_refined_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined.vcf"),
           coverage=os.path.join(alignment_output_dir, "{sample}_{tech}.coverage.txt"),
           vcf_file_list=os.path.join(iris_refined_output_dir, "{sample}_{tech}_sniffles", "vcf_list_dupToIns_irisRefined.txt")
    output: vcf=os.path.join(specific_marked_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined_markedSpec.vcf"),
            vcf_file_list=os.path.join(specific_marked_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "vcf_list_dupToIns_irisRefined_markedSpec.txt")
    threads: lambda wc: min(cluster_config.get("sensitive_ins_to_dup_conversion", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), jasmine_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(specific_marked_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined_markedSpec.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 4000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        output_dir=lambda wc: os.path.join(specific_marked_output_dir, wc.sample + "_" + wc.tech + "_sniffles"),
        min_support_fixed=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_READS_FIXED, 10),
        min_support_fraction=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_READS_FRACTION, 0.25),
        min_length=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_LEN, 30),
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
    run:
        min_support=get_min_support(input.coverage, params.min_support_fixed, params.min_support_fraction)
        shell("{params.java} -cp {params.java_src} Main file_list={input.vcf_file_list} --preprocess_only --mark_specific out_dir={params.output_dir} spec_reads=" + str(min_support) + "spec_len={params.min_length} out_file=test.vcf &> {log}")

rule refined_sensitive_new_sv_types:
    input: vcf=os.path.join(ins_to_dup_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns.vcf"),
           vcf_file_list=os.path.join(ins_to_dup_output_dir, "{sample}_{tech}_sniffles", "vcf_list_dupToIns.txt"),
           bam=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam"),
           bam_bai=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam.bai"),
           bam_file_list=os.path.join(ins_to_dup_output_dir, "{sample}_{tech}_sniffles", "bams.txt"),
    output: vcf=temp(os.path.join(iris_refined_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined.vcf")),
            vcf_file_list=temp(os.path.join(iris_refined_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "vcf_list_dupToIns_irisRefined.txt"))
    threads: lambda wc: min(cluster_config.get("refined_sensitive_new_sv_types", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), iris_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(iris_refined_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: iris_config.get(utils.MEM_MB_CORE, 4000) + iris_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        output_dir=lambda wc: os.path.join(iris_refined_output_dir, wc.sample + "_" + wc.tech + "_sniffles"),
        iris_output_dir=lambda wc: os.path.join(iris_refined_output_dir, wc.sample + "_" + wc.tech + "_sniffles", "iris"),
        samtools=samtools_config.get(utils.PATH, "samtools"),
        ref_genome=config[utils.REFERENCE],
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        minimap2=minimap2_config.get(utils.PATH, "minimap2"),
        racon=racon_config.get(utils.PATH, "racon"),
    shell:
        "{params.java} -cp {params.java_src} Main file_list={input.vcf_file_list} --run_iris --preprocess_only genome_file={params.ref_genome} bam_list={input.bam_file_list} "
        "--iris_args=minimap_path={params.minimap2},racon_path={params.racon},samtools_path={params.samtools},threads={threads},out_dir={params.iris_output_dir} out_dir={params.output_dir} out_file=test.vcf &> {log}"


rule create_bam_file_list:
    output: temp(os.path.join(ins_to_dup_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "bams.txt"))
    input: os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam")
    run:
        with open(output[0], "wt") as dest:
            print(input[0], file=dest)


rule sensitive_ins_to_dup_conversion:
    input: vcf=os.path.join(raw_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".vcf"),
           vcf_file_list=os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles", "vcf_list.txt")
    output: vcf=temp(os.path.join(ins_to_dup_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns.vcf")),
            vcf_file_list=temp(os.path.join(ins_to_dup_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "vcf_list_dupToIns.txt"))
    threads: lambda wc: min(cluster_config.get("sensitive_ins_to_dup_conversion", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), jasmine_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(ins_to_dup_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 4000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        output_dir=lambda wc: os.path.join(ins_to_dup_output_dir, wc.sample + "_" + wc.tech + "_sniffles"),
        ref_genome=config[utils.REFERENCE],
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        samtools=samtools_config.get(utils.PATH, "samtools"),
        max_dup_length=jasmine_config.get(utils.INS_TO_DUP, {}).get(utils.MAX_DUP_LENGTH, 10000),
    shell:
         "{params.java} -cp {params.java_src} Main file_list={input.vcf_file_list} --dup_to_ins genome_file={params.ref_genome} "
         "--preprocess_only out_dir={params.output_dir} threads={threads} samtools_path={params.samtools} max_dup_length={params.max_dup_length} out_file=test.vcf &> {log}"

rule create_first_vcf_file_list:
    output: temp(os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "vcf_list.txt"))
    input: os.path.join(raw_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".vcf")
    run:
        with open(output[0], "wt") as dest:
            print(input[0], file=dest)

localrules: create_first_vcf_file_list, create_bam_file_list

include: "call_svs_sniffles_single.snakefile"
