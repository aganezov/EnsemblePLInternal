import os

files = []
out_dir = config["out_dir"]
suffix = config.get("name", "regions")
input_by_base = {}
for path in config["bams"]:
    basename = os.path.basename(path)
    base = os.path.splitext(basename)[0]
    input_by_base[base] = path
    files.append(os.path.join(out_dir, f"{base}.{suffix}.bam"))
    files.append(os.path.join(out_dir, f"{base}.{suffix}.bam.bai"))


rule all:
    input: files


rule index_bam:
    output: os.path.join(out_dir, "{base}.{suffix," + suffix + "}.bam.bai")
    input: os.path.join(out_dir, "{base}.{suffix}.bam")
    shell:
        "samtools index {input}"


rule create_cut_bam:
    output: os.path.join(out_dir, "{base}.{suffix," + suffix + "}.bam")
    input: bam=lambda wc: input_by_base[wc.base],
           bed=config["regions"]
    shell:
        "bedtools intersect -abam {input.bam} -b {input.bed} > {output}"