#!/usr/bin/env python3
#test_getBag.py
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


@pytest.mark.getBag
def test_getBag():
    assert os.path.exists(os.path.join(
        test_output_path, 'Q-Y5F6_inputBag.zip'))
