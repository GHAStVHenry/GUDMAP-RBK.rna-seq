#!/usr/bin/env python3
import pytest


def pytest_addoption(parser):
    parser.addoption("--filename", action="store")

@pytest.fixture(scope='session')
def filename(request):
    filename_value = request.config.option.filename
    if filename_value is None:
        pytest.skip()
    return filename_value