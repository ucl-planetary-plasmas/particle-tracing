# -*- coding: utf-8 -*-
"""
Created on 2021-10-23 

@author: I Kit Cheng

Title: test_load_data.py
"""
import pytest
import numpy as np
from ..load_data import get_keys


def test_get_keys():
    x = np.array(
        [("Rex", 9, 81.0), ("Fido", 3, 27.0)],
        dtype=[("name", "U10"), ("age", "i4"), ("weight", "f4")],
    )
    assert get_keys(x) == ("name", "age", "weight")
