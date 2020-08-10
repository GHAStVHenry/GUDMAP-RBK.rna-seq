#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
                '/../../'

@pytest.mark.consistencySE
def test_consistencySE():
    assert os.path.exists(os.path.join(test_output_path, 'SE_multiqc_data.json'))
    assert readAssigned("assignedSE.txt","assignedExpectSE.txt")

@pytest.mark.consistencyPE
def test_consistencyPE():
    assert os.path.exists(os.path.join(test_output_path, 'PE_multiqc_data.json'))
    assert readAssigned("assignedPE.txt","assignedExpectPE.txt")

def readAssigned(fileAssigned,fileExpectAssigned):
    data = False
    assigned = open(fileAssigned, "r")
    expect = open(fileExpectAssigned, "r")
    lineAssigned = assigned.readline()
    lineExpect = expect.readline()
    if lineAssigned.strip() < (lineExpect.strip()+(lineExpect.strip()*0.00001)) and lineAssigned.strip() > (lineExpect.strip()-(lineExpect.strip()*0.00001)):
        data = True

    return data
