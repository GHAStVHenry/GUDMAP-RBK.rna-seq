#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
	'/../../'
logs_output_path = os.path.dirname(os.path.abspath(__file__)) + \
	'/../../'


@pytest.mark.dedupData
def test_dedupData_se():
	assert os.path.exists(os.path.join(data_output_path, '16-1ZX4.deduped.bm'))
	assert os.path.exists(os.path.join(data_output_path, '16-1ZX4.deduped.Metrics.txt'))

@pytest.mark.dedupData
def test_dedupData_pe():
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.deduped.bm'))
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.deduped.Metrics.txt'))


@pytest.mark.dedupLogs
def test_dedupLogs_se():
	assert os.path.exists(os.path.join(logs_output_path, '16-1ZX4.dedup.err'))
	assert os.path.exists(os.path.join(logs_output_path, '16-1ZX4.dedup.out'))

@pytest.mark.dedupLogs
def test_dedupLogs_pe():
	assert os.path.exists(os.path.join(logs_output_path, 'Q-Y5JA.dedup.err'))
	assert os.path.exists(os.path.join(logs_output_path, 'Q-Y5JA.dedup.out'))
