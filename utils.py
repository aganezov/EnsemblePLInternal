from collections import defaultdict

DEFAULT_THREAD_CNT = 24

OUTPUT_DIR = "output_dir"
ALIGNMENTS = "alignments"
READS_PATHS = "reads_paths"
SAMPLES = "samples"
LOG = "log"
TOOLS = "tools"
THREADS = "threads"
SAMTOOLS = "samtools"

AWK = "awk"
NGMLR = "ngmlr"
PATH = "path"
TMP_DIR = "tmp_dir"
TECH = "tech"
REFERENCE = "ref"


RAW = "raw"

SNIFFLES = "sniffles"
MIN_SUPPORT = "min_support"
MIN_LENGTH = "min_length"
MAX_NUM_SPLIT_READS = "max_num_splits"
MAX_DISTANCE = "max_distance"
NUM_READS_REPORT = "num_reads_report"
MIN_SEQ_SIZE = "min_seq_size"

SV_TOOLS_ENABLED = "sv_tools_enabled"


def ensure_samples_correctness(config):
    if SAMPLES not in config or not isinstance(config[SAMPLES], dict) or len(config[SAMPLES]) < 1:
        raise ValueError("Configuration data.yaml file is missing information about samples or the setup is not dictionary-like")


def get_samples_to_reads_paths(config):
    sample_to_reads_paths = defaultdict(list)
    for sample_name, sample_data in config["samples"].items():
        if READS_PATHS not in sample_data or not isinstance(sample_data[READS_PATHS], list) or len(sample_data[READS_PATHS]) < 1:
            raise ValueError(
                f"Error when parsing reads paths for sample {sample_name} sample. Make sure the entries are formatted as a list of strings under the {READS_PATHS} key")
        for read_path in sample_data[READS_PATHS]:
            if not read_path.endswith(("fastq", "fq", "fastq.gz", "fq.gz")):
                raise ValueError(f"Unsupported input format for read path {read_path}. Only 'fastq', 'fq', 'fastq.gz', and 'fq.gz' are supported")
            sample_to_reads_paths[sample_name].append(read_path)
        if TECH not in sample_data or sample_data[TECH].lower() not in ["ont", "pb", "pacbio"]:
            raise ValueError(
                f"incorrect or missing tech {sample_data[TECH]} specified for sample {sample_name} in data.yaml. Only ONT or PB are supported, and tech specification is required")
    return sample_to_reads_paths


def ensure_ref_correctness(config):
    if REFERENCE not in config:
        raise ValueError(f"No reference fasta file specified under 'ref' key in data.yaml. Reference is required.")


def get_sniffles_sens_suffix(config):
    min_support = config.get(TOOLS, {}).get(SNIFFLES, {}).get(MIN_SUPPORT, 2)
    min_length = config.get(TOOLS, {}).get(SNIFFLES, {}).get(MIN_LENGTH, 20)
    return f"s{min_support}l{min_length}"


SUPPORTED_SV_TOOLS = {"sniffles"}


def ensure_enabled_sv_tools(config):
    for tool in config[SV_TOOLS_ENABLED]:
        if tool.lower() not in SUPPORTED_SV_TOOLS:
            raise ValueError(f"Attempt to enable unsupported SV inference tool {tool}. Only {','.join(SUPPORTED_SV_TOOLS)} are supported")
