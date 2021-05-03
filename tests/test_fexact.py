import numpy as np
print(np.__version__)
from fexact import fexact

def test_fexact_small():

    k=np.array([[5,1,5],[1,1,5]], dtype=np.int64)
    res=fexact(k)
    assert np.abs(0.4004524886877815 -res)<10e-5

def test_fexact_small_ws():

    k=np.array([[5,1,5],[1,1,5]], dtype=np.int64)
    res=fexact(k,workspace=10000000)
    assert np.abs(0.4004524886877815 -res)<10e-5

def test_fexact_small_ws_kw():

    k=np.array([[5,1,5],[1,1,5]], dtype=np.int64)
    res=fexact(k,workspace=10000000)
    assert np.abs(0.4004524886877815 -res)<10e-5


def test_fexact_hybrid():
    k=np.array([[10,6,7],[5,7,9]], dtype=np.int64)
    res=fexact(k,workspace=10000000,hybrid=True)
    assert np.abs(0.3962 -res)<10e-3

    res=fexact(k,workspace=10000000,hybrid=False)
    assert np.abs(0.4158 -res)<10e-3


import gc
def test_memleak():
    for i in range(1):
        k=np.array([[10,6,7]*10,[5,7,9]*10], dtype=np.int64)
        res=fexact(k,workspace=10000000,hybrid=False)
        print(res)
        gc.collect()

