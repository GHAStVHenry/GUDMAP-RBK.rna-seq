#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.getData
def test_getData():
    assert os.path.exists(os.path.join(
        test_output_path, '16-DNBG_inputBag/bagit.txt'))
    assert os.path.exists(os.path.join(
        test_output_path, '16-DNBG_inputBag/data/assets/Study/16-DMQG/Experiment/16-DMW0Replicate/16-DNBG/22967_gene_counts.txt'))
