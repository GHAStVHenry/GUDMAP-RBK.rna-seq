#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
                '/../../'

@pytest.mark.parseMetadata
def test_parseMetadata():
    assert os.path.exists(os.path.join(test_output_path, 'design.csv'))
    assert readLine(os.path.join(test_output_path, 'design.csv'))

def readLine(fileName):
    data = False
    file = open(fileName, "r")
    if file.readlines() == "uk,se,no,no,Homo sapiens,Experiment_RID,Study_RID,Replicate_RID":
        data = True

    return data
