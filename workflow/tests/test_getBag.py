#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.getBag
def test_getBag():
    assert os.path.exists(os.path.join(
        test_output_path, '16-WJRA_inputBag.zip'))
