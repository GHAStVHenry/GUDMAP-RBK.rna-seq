#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
                '/../../'
test_sampleSE = '16-1ZX4'
test_assignededSE = '7742416'
test_samplePE = 'Q-Y5JA'
test_assignedPE = '2599149'

@pytest.mark.consistencySE
def test_consistencySE():
    assert os.path.exists(os.path.join(test_output_path, test_sampleSE, '_multiqc_data.json'))
    assert readAssigned("assignedSE.txt",test_assignedSE)

@pytest.mark.consistencyPE
def test_consistencyPE():
    assert os.path.exists(os.path.join(test_output_path, test_samplePE, '_multiqc_data.json'))
    assert readAssigned("assignedPE.txt",test_assignedPE)

def readAssgned(assignedFile,assigned):
    data = False
    file = open(assignedFile, "r")
    line = file.readline()
    if line.strip() == assigned:
        data = True

    return data
