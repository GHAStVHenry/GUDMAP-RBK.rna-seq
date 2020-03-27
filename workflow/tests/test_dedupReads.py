#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
	'/../../'


@pytest.mark.dedupData
def test_dedupData():
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.deduped.bam'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.deduped.bam.bai'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.deduped.chr1.bam'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.deduped.chr1.bam.bai'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.deduped.chr2.bam'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.deduped.chr2.bam.bai'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.deduped.chr8.bam'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.deduped.chr8.bam.bai'))