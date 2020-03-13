#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
	'/../../'


@pytest.mark.dedupData
def test_dedupData():
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.deduped.bam'))
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.deduped.bam.bai'))
