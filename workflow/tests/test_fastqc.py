#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.fastqc
def test_fastqc():
    assert os.path.exists(os.path.join(
        test_output_path, '16-WJRA_1M.R1_fastqc.zip'))
