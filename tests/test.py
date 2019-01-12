import os
import subprocess

import pytest

TEST_DATA_DIR = 'tests/data/'
EXPECTED_DIR = os.path.join(TEST_DATA_DIR, 'expected')

DNA_NORMAL_VAF_FP = os.path.join(TEST_DATA_DIR, 'dna_a.vaf')
DNA_TUMOR_VAF_FP = os.path.join(TEST_DATA_DIR, 'dna_t.vaf')
RNA_NORMAL_VAF_FP = os.path.join(TEST_DATA_DIR, 'rna_a.vaf')
RNA_TUMOR_VAF_FP = os.path.join(TEST_DATA_DIR, 'rna_t.vaf')

def test_simple():
    tool_args = ['python', 'phyllite/phyllite.py',
            '--dna-normal-vaf', DNA_NORMAL_VAF_FP,
            '--dna-tumor-vaf', DNA_TUMOR_VAF_FP,
            '--rna-normal-vaf', RNA_NORMAL_VAF_FP,
            '--rna-tumor-vaf', DNA_TUMOR_VAF_FP]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert open('output.tsv').read() == open(os.path.join(EXPECTED_DIR, 'simple.tsv')).read()

def test_thresholds():
    tool_args = ['python', 'phyllite/phyllite.py',
            '--dna-normal-vaf', DNA_NORMAL_VAF_FP,
            '--dna-tumor-vaf', DNA_TUMOR_VAF_FP,
            '--rna-normal-vaf', RNA_NORMAL_VAF_FP,
            '--rna-tumor-vaf', DNA_TUMOR_VAF_FP,
            '--min-depth', '1',
            '--max-dna-minor-vaf', '0.05',
            '--min-rna-minor-vaf', '0.09',
            '--table-output', 'output.tsv',
            '--json-output', 'output.json']
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert open('output.tsv').read() == open(os.path.join(EXPECTED_DIR, 'threshold.tsv')).read()
