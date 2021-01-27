#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os
import json

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'

@pytest.mark.completionMultiqc
def test_multiqcExist(filename):
    assert os.path.exists(os.path.join(
        test_output_path, filename))