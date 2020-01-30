#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
	'/../../'


@pytest.mark.alignData
def test_alignData_se():
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.unal.gz'))
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.bam'))
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.se.sorted.bai'))


@pytest.mark.alignData
def test_alignData_pe():
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.pe.unal.gz'))
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.pe.sorted.bam'))
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA_1M.pe.sorted.bai'))


@pytest.mark.alignLogs
def test_alignLogs_se():
	assert os.path.exists(os.path.join(data_output_path, '16-1ZX4.align.err'))
	assert '34497376 reads; of these:' in open(os.path.join(data_output_path, '16-1ZX4.align.err')).readlines()[0]
	assert os.path.exists(os.path.join(data_output_path, '16-1ZX4.align.out'))


@pytest.mark.alignLogs
def test_alignLogs_pe():
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.align.err'))
	assert utils.count_lines(os.path.join(data_output_path, 'Q-Y5JA.align.err')) == 7
	assert '15824858 reads; of these:' in open(os.path.join(data_output_path, 'Q-Y5JA.align.err')).readlines()[0]
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.align.out'))
	assert utils.count_lines(os.path.join(data_output_path, 'Q-Y5JA.align.out')) == 0
