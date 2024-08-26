#import pytest
import numpy as np 


def test_norm(): 
    d = np.loadtxt('./mode_normalisation/status.txt')    
    assert d == 0 
