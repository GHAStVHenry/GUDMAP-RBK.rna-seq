#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.alignData
def test_alignData_se():
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.unal.gz'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.bam'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.se.sorted.bam.bai'))


@pytest.mark.alignData
def test_alignData_pe():
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.pe.unal.gz'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.pe.sorted.bam'))
    assert os.path.exists(os.path.join(
        data_output_path, '16-WJRA_1M.pe.sorted.bam.bai'))
