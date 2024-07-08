import os
import pytest
from subprocess import getstatusoutput, getoutput
import sample_parsing as sp
import glob
import shutil
import pandas as pd


prg = './sample_parsing.py'
sample_id_data = "testdata/sample-id.txt"

# -----------------------------------------------------
##################
# FIXTURES
##################
@pytest.fixture
def create_BC():
    BC02 = sp.BarcodeAttributes("BC02", "3531017", "ENV", "VP1_DENOVO", "testdata/sample-id.txt")
    return BC02

@pytest.fixture
def BC_dir(create_BC):
    return create_BC.dir_path

@pytest.fixture
def tmp_fq(create_BC):
    """make BC directory & concat. fastq"""
    return create_BC.concat_raw_fastq

@pytest.fixture
def readsize(tmp_fq, create_BC):
    fastq = glob.glob(create_BC.dir_path + "/*.fastq.gz")
    return sp.read_lengths(fastq)

@pytest.fixture
def sumstats(create_BC, readsize):
    return sp.summary_statistics(readsize, create_BC)

# -----------------------------------------------------
##################
# TESTS
##################

def test_exists():
    """exists"""
    assert os.path.isfile(prg)

# -----------------------------------------------------
def test_usage():
    """usage"""
    for flag in ['-h', '--help']:
        rv, out = getstatusoutput(f'{prg} {flag}')
        assert rv == 0
        assert out.lower().startswith('usage')

# -----------------------------------------------------
def test_runname():
    """get runname"""
    assert sp.get_runname(sample_id_data) == "20220921_X4_NBA_EXP02_IFIK_test96"


# -----------------------------------------------------
def test_dirpath(create_BC):
    """get BC directory path"""
    wd = os.getcwd()
    assert create_BC.dir_path == os.path.join(wd, "BC02")

# -----------------------------------------------------
def test_rawfastq(create_BC): ## change maybe!
    """get raw fastq files per barcode"""
    assert create_BC.get_raw_fastq == glob.glob("/storage/tmp/20220921_X4_NBA_EXP02_IFIK_test96.tmp/fastq_pass/barcode02/*.gz")

# -----------------------------------------------------
def test_concatfastq(tmp_fq, BC_dir):
    """test concatination of raw fastq files"""
    try:
        assert os.path.exists(BC_dir + "/BC02.fastq.gz")
    except Exception:
        assert False
    finally:
         shutil.rmtree(BC_dir)

# -----------------------------------------------------
def test_ampliconlength(create_BC):
    """test amplicon size by gene name"""
    assert create_BC.amplicon_length == 400

# -----------------------------------------------------
def test_readlengths(tmp_fq, create_BC, BC_dir):
    """get read lengths from fastq"""
    try:
        fastq = glob.glob(create_BC.dir_path + "/*.fastq.gz")
        sizes = sp.read_lengths(fastq)
        assert len(sizes) == 27593
    except Exception:
        assert False
    finally:
         shutil.rmtree(BC_dir)

# -----------------------------------------------------
def test_histogram(readsize, create_BC, BC_dir):
    """get read length distribution histogram"""
    try:
        sp.read_lengths_hist(readsize, create_BC)
        assert os.path.exists(BC_dir + "/" + create_BC.BC + ".png")
    except Exception:
        assert False
    finally:
        shutil.rmtree(BC_dir)

# -----------------------------------------------------
def test_summarystats(readsize, create_BC, BC_dir):
    """get summary statistics as table"""
    try:
        sp.summary_statistics(readsize, create_BC)
        assert os.path.exists(BC_dir + "/" + create_BC.BC + "_stats.txt")
    except Exception:
        assert False
    finally:
        shutil.rmtree(BC_dir)

def test_concat_summarystats(sumstats, BC_dir):
    try:
        wd=os.getcwd()
        sp.concat_summarystats()
        assert os.path.exists(wd + "/summary_stats.txt")
    except:
        assert False
    finally:
        shutil.rmtree(BC_dir)
        os.remove(wd + "/summary_stats.txt")



