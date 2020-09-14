#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.trimData
def test_trimData_se():
    assert os.path.exists(os.path.join(
        test_output_path, 'Q-Y5F6_1M.se_trimmed.fq.gz'))


@pytest.mark.trimData
def test_trimData_pe():
    assert os.path.exists(os.path.join(
        test_output_path, 'Q-Y5F6_1M.pe_R1_val_1.fq.gz'))
    assert os.path.exists(os.path.join(
        test_output_path, 'Q-Y5F6_1M.pe_R2_val_2.fq.gz'))
