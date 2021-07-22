#!/usr/bin/env python3
#test_getData.py
#*
#* --------------------------------------------------------------------------
#* Licensed under MIT (https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/blob/14a1c222e53f59391d96a2a2e1fd4995474c0d15/LICENSE)
#* --------------------------------------------------------------------------
#*

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.getData
def test_getData():
    assert os.path.exists(os.path.join(
        test_output_path, 'Q-Y5F6_inputBag/bagit.txt'))
    assert os.path.exists(os.path.join(
        test_output_path, 'Q-Y5F6_inputBag/data/assets/Study/Q-Y4GY/Experiment/Q-Y4DP/Replicate/Q-Y5F6/mMARIS_Six2-#3.gene.rpkm.txt'))
