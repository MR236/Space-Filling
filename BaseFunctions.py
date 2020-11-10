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


def Polya_Dissection(n, angle):
    TriangleN = [[np.array((0,0)), np.array((2*np.cos(angle)**2, 2*np.sin(angle)*np.cos(angle))), np.array((2,0))]]
    for i in range(n):
        NewTriangles = []
        for x in TriangleN:
            a = np.linalg.norm(x[1]-x[0])
            c = np.linalg.norm(x[2]-x[0])
            r = a**2/c
            bisect_p = np.add(np.multiply(r/c, x[2]-x[0]), x[0])
            NewTriangles.append([x[0], bisect_p, x[1]])
            NewTriangles.append([x[1], bisect_p, x[2]])
        TriangleN = NewTriangles
    Triangles = []
    for i in TriangleN:
        points = []
        points.append([i[0][0], i[0][1]])
        points.append([i[1][0], i[1][1]])
        points.append([i[2][0], i[2][1]])
        points.append([i[0][0], i[0][1]])
        Triangles.append(points)
    return Triangles


def Polya_Curve(q, angle):
    e0j = 0
    e3j = 0
    nj = 0
    dj = 0
    x = 2*np.cos(angle)**2
    y = 2*np.sin(angle)*np.cos(angle)
    s = [np.array([0,0]), np.array([x*2, 0]), np.array([x*2, y*2]), np.array([x*2, 0])]
    S1 = np.array([[0.0, -1.0],
                   [1.0, 0.0]])
    image = np.array([0.0,0.0])
    for j in range(1, len(q) + 1):
        vec = np.matmul(np.linalg.matrix_power(S1, dj%4), s[q[j - 1]])
        coef = (1/2**j)*(y**dj)*(x**e0j)*((2-x)**e3j)*((-1)**(nj%2))
        image += np.multiply(coef, vec)
        if q[j-1] == 0:
            e0j += 1
        elif q[j-1] == 1:
            dj += 1
        elif q[j-1] == 2:
            dj += 1
            nj += 1
        else:
            e3j += 1
    return [image[0], image[1]]


def Polya_Polygon(n, angle):
    nodal_points = []
    for i in range(0, 2 ** n):
        nodal_points.append(Polya_Curve(quaternary(i / 2 ** n), angle))
    nodal_points.append([2, 0])
    return nodal_points

