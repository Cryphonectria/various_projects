import argparse
import tempfile
import re
import os
from subprocess import getstatusoutput, getoutput
import glob
import pandas as pd
import pytest
import shutil
import check_contaminations as cc

prg = './check_contaminations.py'

@pytest.fixture
def create_dummy_files():
    # Create dummy subdirectories and 'list.mapped.taxo' files in the current directory
    tmpdirs = []

    for i in range(3):
        tmpdir = tempfile.mkdtemp(dir=os.getcwd())
        tmpdirs.append(tmpdir)
        tmp_file = tempfile.NamedTemporaryFile(suffix='list.mapped.taxo', dir=tmpdir, delete=False)
        tmp_file.close()

    yield tmpdirs

    # Clean up: Remove the temporary directories and their contents after the test
    for tmpdir in tmpdirs:
        shutil.rmtree(tmpdir)



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

def test_find_taxo_files(create_dummy_files):
    tmpdir = create_dummy_files

    # Call the function with the temporary directory
    taxo_files = cc.find_taxo_files(os.getcwd())
    print(tmpdir)
    print(taxo_files)
    # Check if the correct number of files are found
    assert len(taxo_files) == 3

    # Check if all files have the expected path
    expected_paths = [
        os.path.join(tmpdir, f"subdir_{i}/list.mapped.taxo") for i in range(3)
    ]
    assert sorted(taxo_files) == sorted(expected_paths)
    