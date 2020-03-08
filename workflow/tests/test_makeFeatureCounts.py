#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
	'/../../'


@pytest.mark.makeFeatureCounts
def test_makeFeatureCounts():
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.featureCounts'))
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.countTable.csv'))
