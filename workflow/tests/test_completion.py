#!/usr/bin/env python3

import pytest
import pandas as pd
from io import StringIO
import os
import json

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
    '/../../'


@pytest.mark.completionIntegration_se
def test_multiqc():
    assert os.path.exists(os.path.join(
        test_output_path, 'SE_multiqc_data.json'))


@pytest.mark.completionIntegration_pe
def test_multiqc():
    assert os.path.exists(os.path.join(
        test_output_path, 'PE_multiqc_data.json'))


@pytest.mark.completionOverride_inputBag
def test_multiqc():
    assert os.path.exists(os.path.join(
        test_output_path, 'inputBagOverride_PE_multiqc_data.json'))


@pytest.mark.completionOverride_fastq
def test_multiqc():
    assert os.path.exists(os.path.join(
        test_output_path, 'fastqOverride_PE_multiqc_data.json'))


@pytest.mark.completionOverride_species
def test_multiqc():
    assert os.path.exists(os.path.join(
        test_output_path, 'speciesOverride_PE_multiqc_data.json'))


@pytest.mark.completionOverride_stranded
def test_multiqc():
    assert os.path.exists(os.path.join(
        test_output_path, 'strandedOverride_PE_multiqc_data.json'))


@pytest.mark.completionOverride_spike
def test_multiqc():
    assert os.path.exists(os.path.join(
        test_output_path, 'spikeOverride_PE_multiqc_data.json'))