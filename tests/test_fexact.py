import numpy as np
print(np.__version__)
from fexact import fexact

def test_fexact_small():

    k=np.array([[5,1,5],[1,1,5]], dtype=np.intc)
    res=fexact(k)
    assert np.abs(0.4004524886877815 -res)<10e-5
