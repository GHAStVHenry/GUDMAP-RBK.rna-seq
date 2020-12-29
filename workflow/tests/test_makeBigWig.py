#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.makeBigWig
def test_makeBigWig():
    assert os.path.exists(os.path.join(data_output_path, '16-DNBG_1M.se.bw'))
