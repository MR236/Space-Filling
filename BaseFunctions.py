import numpy as np

def quaternary(n):
    k = n
    nums = []
    i = 1
    if k == 0:
        return [0]
    while k != 0:
        q = int(np.floor(k*(4**i)))
        k -= q/(4**i)
        nums.append(q)
        i += 1
    return nums


def sgn(q):
    if q == 0:
        return 0
    else:
        return 1


def Hilbert_Curve(q):
    x = 0
    y = 0
    e0j = 0
    e3j = 0
    dj = 0
    for j in range(1, len(q) + 1):
        x += ((1/2)**j)*((-1)**e0j)*(sgn(q[j - 1]))*((1-dj)*q[j - 1] - 1)
        y += ((1/2)**j)*((-1)**e0j)*(sgn(q[j - 1]))*(1-dj*q[j-1])
        if q[j - 1] == 0:
            e0j = (e0j+1)%2
            dj = (e0j+e3j)%2
        elif q[j - 1] == 3:
            e3j = (e3j+1)%2
            dj = (e0j+e3j)%2
    return [x, y]


def Hilbert_Polygon(n):
    nodal_points = []
    for i in range(0, 4**n):
        nodal_points.append(Hilbert_Curve(quaternary(i/4**n)))
    nodal_points.append([1, 0])
    return nodal_points


def Sierpinski_Curve(q):
    nj = 0
    dj = 0
    x = 0
    y = 0
    for j in range(1, len(q) + 1):
        x += (((-1)**nj)/(2**j))*sgn(q[j-1])*((1-dj)*(1+(-1)**dj)+0.5*(dj-2)*(1-(-1)**dj)*(1+(-1)**q[j-1]))
        y += (((-1)**nj)/(2**j))*sgn(q[j-1])*((2-dj)*(1-(-1)**dj)+0.5*(1-dj)*(1+(-1)**dj)*(1+(-1)**q[j-1]))
        if q[j-1] == 2:
            nj = (nj+1)%2
        if q[j-1] == 1 or q[j-1] == 2:
            dj = (dj+1)%4
    return [x,y]

def Sierpinski_Polygon(n):
    nodal_points = []
    for i in range(0, 2**n):
        nodal_points.append(Sierpinski_Curve(quaternary(i/2**n)))
    nodal_points.append([2, 0])
    return nodal_points