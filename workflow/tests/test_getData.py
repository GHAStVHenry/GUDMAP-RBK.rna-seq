#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
                '/../../'

@pytest.mark.getData
def test_getData():
    assert os.path.exists(os.path.join(test_output_path, 'Replicate_16-1ZX4/bagit.txt'))
    assert os.path.exists(os.path.join(test_output_path, 'Study_Q-Y4H0/data/assets/Study/Q-Y4H0/Experiment/Q-Y4BY/Replicate/Q-Y5F8/hMARIS_SIX2+_RiboDep#1.gene.rpkm.txt'))