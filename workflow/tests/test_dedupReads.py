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
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.deduped.chrNC_000001.11.bam'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.deduped.chrNC_000001.11.bam.bai'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.deduped.chrNC_000002.12.bam'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.deduped.chrNC_000002.12.bam.bai'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.deduped.chrNC_000019.10.bam'))
    assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.deduped.chrNC_000019.10.bam.bai'))