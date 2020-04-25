#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
                '/../../'

@pytest.mark.downsampleData
def test_downsampleData():
    assert os.path.exists(os.path.join(test_output_path, 'sampled.1.fq'))