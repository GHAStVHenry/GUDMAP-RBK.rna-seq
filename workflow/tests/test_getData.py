#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
                '/../../Replicate_16-1ZX4/'

@pytest.mark.getData
def test_getData():
    assert os.path.exists(os.path.join(test_output_path, 'bagit.txt'))
    assert os.path.exists(os.path.join(test_output_path, 'data/assets/Study/16-1ZWP/Experiment/16-1ZWR/Replicate/16-1ZX4/3_1_single.R1.fastq.gz'))