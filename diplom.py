from pprint import pprint
from scipy import misc, integrate, special
from scipy.fftpack import diff, sc_diff
import sympy as sp
import numpy as np
import matplotlib as plt
import pyquaternion as pq

def _Lambda(t):
    return pq.Quaternion(    
        np.array(
            [
                t,
                t,
                t,
                t
            ]
        )    
    )

def _w(t):
    return pq.Quaternion(    
        np.array([
            w0,
            t,
            t,
            t
            ]))

def _theta(t):
    return integrate.quad(lambda t: b[2]*_w(t), 0, t)

def _M(t):
    return np.array(
        [
            0,
            t,
            t,
            t
        ]
    )
def u(t):
    return 0

def _J(a):
    return integrate.quad(lambda t: (1+(I[1]*u(t))**2), 0, T)

def _beta(t, tau):
    return np.exp(pq.Quaternion(
        0, 
        b[2]*integrate.quad(
            integrate.quad(
                _M(tau)[1]/I[0]+_w(0)[1],
                0, 
                tau
                ), 
            0, 
            t
            ),
        0,
        0
    )
)

def _B(t):
    return pq.Quaternion(np.cos(_theta(t/2)), np.sin(_theta(t/2)))

def _dev_w():
    return pq.Quaternion(
        m[0],
        b[1]*m[1]-b[0]*_w(t)[1]*_w(t)[2],
        b[1]*m[2]-b[0]*_w(t)[1]*_w(t)[3]
    )

eps = 0.001
n = 10
I_1 = 1
I_2 = 1
T = 1
t = np.linspace(0, T, n)
I = np.array(
    [
        I_1, 
        I_2
    ]
)
c = I[0]/I[1]
b = np.array(
    [
        c-1,
        c,
        (c-1)/c
    ]
)
m = np.array(
    [
        _M(t)[1]/I[0],
        _M(t)[2]/I[0],
        _M(t)[3]/I[0]
    ]
)
w0 = pq.Quaternion(0,0,0,0)
wT = _B(_theta(T)).inverse*pq.Quaternion(0,0,0,0)*_B(_theta(T))

