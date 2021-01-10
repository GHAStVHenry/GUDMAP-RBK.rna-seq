#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os
import json

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.consistencySE
def test_consistencySE():
    assert os.path.exists(os.path.join(
        test_output_path, 'SE_multiqc_data.json'))

    with open(os.path.join(
        test_output_path, 'SE_multiqc_data.json')) as f:
        assigned_reads_json = json.load(f)
    assigned_reads = assigned_reads_json['report_general_stats_data'][4]['16-1ZX4_sorted']['Assigned']
    baseline = 7746121
    baseline_hi = baseline+(baseline*0.01)
    baseline_lo = baseline-(baseline*0.01)
    self.assertTrue ( baseline_lo <= assigned_reads <=  baseline_hi )


@pytest.mark.consistencyPE
def test_consistencyPE():
    assert os.path.exists(os.path.join(
        test_output_path, 'PE_multiqc_data.json'))

    with open(os.path.join(
        test_output_path, 'PE_multiqc_data.json')) as f:
        assigned_reads_json = json.load(f)
    assigned_reads = assigned_reads_json['report_general_stats_data'][4]['Q-Y5JA_sorted']['Assigned']
    baseline = 2596053
    baseline_hi = baseline+(baseline*0.01)
    baseline_lo = baseline-(baseline*0.01)
    self.assertTrue ( baseline_lo <= assigned_reads <=  baseline_hi )
