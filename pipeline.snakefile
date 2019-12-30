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
utils.ensure_enabled_sv_tools(config)


# during development this thing guarantees that only the latest supported part of pipeline produces results
overall_expected_files = []
for (sample, tech) in sample_to_reads_paths.keys():
    for sv_tool in config[utils.SV_TOOLS_ENABLED]:
        if sv_tool == "sniffles":
            suffix = utils.get_sniffles_sens_suffix(config) + "."
        else:
            suffix = ""
        overall_expected_files.append(os.path.join(alignment_output_dir, f"{sample}_{tech}.coverage.txt"))
        overall_expected_files.append(os.path.join(raw_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}vcf"))

rule main:
    input: overall_expected_files

include: "call_svs_sniffles_single.snakefile"
include: "align_single.snakefile"