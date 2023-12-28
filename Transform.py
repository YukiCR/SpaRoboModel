"""
providing basic transformation needed
"""
import numpy as np

def augment(R:np.ndarray=np.eye(3), p:np.ndarray=np.zeros(3)) -> np.ndarray: 
    """
    get the augumented matrix T, given R or p or both of them\n
    aug(R,p) equals T(R)@T(p)
    """
    return np.vstack(( np.hstack((R, p[:,np.newaxis])),  np.array([0,0,0,1]) ))

def augP(p:np.ndarray) -> np.ndarray:
    """
    augment a P array to 4 dimension
    """
    return np.hstack( (p, np.array([1])) )

def Rx(x):
    pass

def Ry(y):
    pass

def Rz(theta):
    """
    rotate along the z axis for angle theta, \n
    theta in radians
    """
    return np.array([[np.cos(theta),   -np.sin(theta),   0],
                     [np.sin(theta),    np.cos(theta),   0],
                     [       0     ,         0       ,   1]])
    