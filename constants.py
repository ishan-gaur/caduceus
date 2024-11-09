from pathlib import Path
from Bio import SeqIO

OUTPUT_DIR = Path(__file__).parent / "output"
if not OUTPUT_DIR.exists():
    OUTPUT_DIR.mkdir()

SAMPLES = ["LIB", "S1", "S2", "S3", "TEST", "LIB_TEST", "LIB_TEST_HEAD", "LIB_TEST_TAIL"]
SELECTION_SAMPLES = ["S1", "S2", "S3"]
REFERENCE = "LIB"
SAMPLE_LET = {"LIB":"A", "S1":"B", "S2":"C", "S3":"D", "TEST":"D", "LIB_TEST":"A", "LIB_TEST_HEAD":"A", "LIB_TEST_TAIL":"A"}
SAMPLE_NUM = {"LIB":"1", "S1":"2", "S2":"3", "S3":"4", "TEST":"T", "LIB_TEST":"L", "LIB_TEST_HEAD":"LH", "LIB_TEST_TAIL":"LT"}

def get_sample_files(sample):
    data_dir = Path("/home/ishan/mg_analysis/data/microglia_promoter_screen/")
    r1 = data_dir / f"2024_DSSV001{SAMPLE_LET[sample]}_S{SAMPLE_NUM[sample]}_L008_R1_001.fastq"
    r2 = data_dir / f"2024_DSSV001{SAMPLE_LET[sample]}_S{SAMPLE_NUM[sample]}_L008_R2_001.fastq"
    return r1, r2

TFBS_FW = [
    "CGGA",
    "GTAT",
    "TTAT",
    "GAAG",
    "CGTC",
    "TAGC",
    "GTGC",
    "CGAC",
    "CACG",
]
TFBS_BW = [
    "CGGG",
    "AAGC",
    "AGCA",
    "GGAA",
    "CGGC",
    "TGGC",
    "AATT",
    "TAGG",
    "AATA",
]

TFBS_LABELS = [
    "ZFX",
    "SP1/SP2",
    "SP3",
    "KLF5",
    "FOSL1",
    "FOSL2",
    "JUNB",
    "JUND",
    "BACH2",
]

TFBS_GENES = [
    "ZFX_A",
    "SP1/SP2_A",
    "SP3_A",
    "KLF5_A",
    "FOSL1_A",
    "FOSL2_A",
    "JUNB_A",
    "JUND_A",
    "BACH2_A",
    "ZFX_B",
    "SP1/SP2_B",
    "SP3_B",
    "KLF5_B",
    "FOSL1_B",
    "FOSL2_B",
    "JUNB_B",
    "JUND_B",
    "BACH2_B"
]

TFBS = TFBS_FW + TFBS_BW
TFBS = {tfbs: i for (i, tfbs) in enumerate(TFBS)}

assert len(TFBS_FW) == len(TFBS_BW)
assert len(TFBS_FW) == len(TFBS_LABELS)

# amplicon_template_file = Path(__file__).parent / 'elips-amplicon-example.fasta'
# amplicon_template = list(SeqIO.parse(amplicon_template_file, 'fasta'))[0].seq

# p5/7 variable regions define the sample
# SAMPLE_P5 = {
#     "LIB": "TTCGCAGT",
#     "S1": "CGAGACTA",
#     "S2": "ACAGCTCA",
#     "S3": "GGCATACT",
#     "TEST": "GGCATACT",
# }

# SAMPLE_P7 = {
#     "LIB": "CGGATTGA",
#     "S1": "ATGTAGCG",
#     "S2": "AGTGGATC",
#     "S3": "TCGTGGAT",
#     "TEST": "TCGTGGAT",
# }

# BC_P5_HEAD = "AATGATACGGCGACCACCGAGATCTACAC"
# p5_bc_len = 10
BC_TRUSEQ1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
BC_TRUSEQ1 = BC_TRUSEQ1[-6:] # because we actually don't get long enough reads to get the whole thing
# umi1_len = 12 # actually 12-16
BC_10X = "TTGCTAGGACCGGCCTTAAAGC"
GFP_TAG = "TTACTTGTACAGCTCGTCCATGCC"
# umi2_len = 4 # actually 0-4
BC_TRUSEQ2 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
BC_TRUSEQ2 = BC_TRUSEQ2[:6] # because we actually don't get long enough reads to get the whole thing
# p7_bc_len = 10
# BC_P7_TAIL = "ATCTCGTATGCCGTCTTCTGCTTG"