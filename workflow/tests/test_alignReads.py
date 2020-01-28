#!/usr/bin/env python3

import pytest
import pandas as pd
import os
import utils

data_output_path = os.path.dirname(os.path.abspath(__file__)) + \
	'/../../'
logs_output_path = os.path.dirname(os.path.abspath(__file__)) + \
	'/../../'

@pytest.mark.alignData
def test_alignData_se():
	assert os.path.exists(os.path.join(data_output_path, '16-1ZX4.unal.gz'))
	assert utils.count_lines(os.path.join(data_output_path, '16-1ZX4.unal.gz')) == 3070528
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.bam'))
	assert utils.count_lines(os.path.join(data_output_path, 'Q-Y5JA.sorted.bam')) == 5805611
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.bai'))
	assert utils.count_lines(os.path.join(data_output_path, 'Q-Y5JA.sorted.bai')) == 12824

@pytest.mark.alignData
def test_alignData_pe():
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.unal.gz'))
	assert utils.count_lines(os.path.join(data_output_path, 'Q-Y5JA.unal.gz')) == 0
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.bam'))
	assert utils.count_lines(os.path.join(data_output_path, 'Q-Y5JA.sorted.bam')) == 5805611
	assert os.path.exists(os.path.join(data_output_path, 'Q-Y5JA.sorted.bai'))
	assert utils.count_lines(os.path.join(data_output_path, 'Q-Y5JA.sorted.bai')) == 12824


@pytest.mark.alignLogs
def test_alignLogs_se():
	assert os.path.exists(os.path.join(logs_output_path, '16-1ZX4.align.err'))
	assert utils.count_lines(os.path.join(logs_output_path, '16-1ZX4.align.err')) == 7
	assert '34497376 reads; of these:' in open(os.path.join(logs_output_path, '16-1ZX4.align.err')).readlines()[0]
	assert os.path.exists(os.path.join(logs_output_path, '16-1ZX4.align.out'))
	assert utils.count_lines(os.path.join(logs_output_path, '16-1ZX4.align.out')) == 0

@pytest.mark.alignLogs
def test_alignLogs_pe():
	assert os.path.exists(os.path.join(logs_output_path, 'Q-Y5JA.align.err'))
	assert utils.count_lines(os.path.join(logs_output_path, 'Q-Y5JA.align.err')) == 7
	assert '15824858 reads; of these:' in open(os.path.join(logs_output_path, 'Q-Y5JA.align.err')).readlines()[0]
	assert os.path.exists(os.path.join(logs_output_path, 'Q-Y5JA.align.out'))
	assert utils.count_lines(os.path.join(logs_output_path, 'Q-Y5JA.align.out')) == 0
