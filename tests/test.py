import os
import subprocess

import pytest

TEST_DATA_DIR = 'tests/data/'
EXPECTED_DIR = os.path.join(TEST_DATA_DIR, 'expected')

DNA_NORMAL_VAF_FP = os.path.join(TEST_DATA_DIR, 'dna_a.vaf')
DNA_TUMOR_VAF_FP = os.path.join(TEST_DATA_DIR, 'dna_t.vaf')
RNA_NORMAL_VAF_FP = os.path.join(TEST_DATA_DIR, 'rna_a.vaf')
RNA_TUMOR_VAF_FP = os.path.join(TEST_DATA_DIR, 'rna_t.vaf')

DNA_NORMAL_COUNT_VAF_FP = os.path.join(TEST_DATA_DIR, 'counts/dna_a.vaf')
DNA_TUMOR_COUNT_VAF_FP = os.path.join(TEST_DATA_DIR, 'counts/dna_t.vaf')
RNA_NORMAL_COUNT_VAF_FP = os.path.join(TEST_DATA_DIR, 'counts/rna_a.vaf')
RNA_TUMOR_COUNT_VAF_FP = os.path.join(TEST_DATA_DIR, 'counts/rna_t.vaf')


def test_simple():
    tool_args = ['python', 'phyllite/phyllite.py',
            '--dna-normal-vaf', DNA_NORMAL_VAF_FP,
            '--dna-tumor-vaf', DNA_TUMOR_VAF_FP,
            '--rna-normal-vaf', RNA_NORMAL_VAF_FP,
            '--rna-tumor-vaf', DNA_TUMOR_VAF_FP,
            '--max-dna-minor-count', '10',
            '--min-rna-minor-count', '1']
    
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
            '--max-dna-minor-count', '10',
            '--min-rna-minor-count', '1',
            '--table-output', 'output.tsv',
            '--json-output', 'output.json']
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert open('output.tsv').read() == open(os.path.join(EXPECTED_DIR, 'threshold.tsv')).read()

def test_count_thresholds():
    tool_args = ['python', 'phyllite/phyllite.py',
            '--dna-normal-vaf', DNA_NORMAL_COUNT_VAF_FP,
            '--dna-tumor-vaf', DNA_TUMOR_COUNT_VAF_FP,
            '--rna-normal-vaf', RNA_NORMAL_COUNT_VAF_FP,
            '--rna-tumor-vaf', DNA_TUMOR_COUNT_VAF_FP,
            '--min-depth', '10',
            '--max-dna-minor-vaf', '0.1',
            '--min-rna-minor-vaf', '0.1',
            '--max-dna-minor-count', '1',
            '--min-rna-minor-count', '3',
            '--table-output', 'output.tsv',
            '--json-output', 'output.json']
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert open('output.tsv').read() == open(os.path.join(EXPECTED_DIR, 'counts.tsv')).read()
