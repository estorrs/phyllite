import os
import subprocess

import pytest

TEST_DATA_DIR = 'tests/data/'

DNA_TUMOR_COUNT_VAF_FP = os.path.join(TEST_DATA_DIR, 'counts/dna_t.vaf')
RNA_TUMOR_COUNT_VAF_FP = os.path.join(TEST_DATA_DIR, 'counts/rna_t.vaf')


def test_simple():
    tool_args = ['python', 'phyllite/phyllite.py',
            '--dna-vaf', DNA_TUMOR_COUNT_VAF_FP,
            '--rna-vaf', RNA_TUMOR_COUNT_VAF_FP,
            '--min-depth', '10',
            '--max-dna-minor-vaf', '0.1',
            '--min-rna-minor-vaf', '0.1',
            '--max-dna-minor-count', '1',
            '--min-rna-minor-count', '3',
            '--table-output', 'output.tsv',
            '--json-output', 'output.json']
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    f = open('output.tsv')
    output = f.read()
    assert '.9' in output and '.5' in output and '\tG\t' in output
