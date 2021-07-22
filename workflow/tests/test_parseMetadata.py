#!/usr/bin/env python3
#test_parseMetadata.py
#*
#* --------------------------------------------------------------------------
#* Licensed under MIT (https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/blob/14a1c222e53f59391d96a2a2e1fd4995474c0d15/LICENSE)
#* --------------------------------------------------------------------------
#*

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
    line = file.readline()
    if line.strip() == "uk,uk,se,unstranded,f,Homo sapiens,75,Experiment_RID,Study_RID,Replicate_RID":
        data = True

    return data
