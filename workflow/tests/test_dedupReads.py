#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.dedupData
def test_dedupData():
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.deduped.bam'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.deduped.bam.bai'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.deduped.chr8.bam'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.deduped.chr8.bam.bai'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.deduped.chr4.bam'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.deduped.chr4.bam.bai'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.deduped.chrY.bam'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.deduped.chrY.bam.bai'))
