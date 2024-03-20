from itertools import combinations
from random import random
from time import time

import matplotlib.pyplot as plt
from numpy import array, abs, sign, sqrt
from numpy.linalg import norm, det
TOL = 1e-8

# Symmetric Polynomials
p1 = lambda a, b, c: a + b + c
p2 = lambda a, b, c: a*b + b*c + c*a
p3 = lambda a, b, c: a*b*c

class Circle:
    def __init__(self, x, y, k):
        self.x = x
        self.y = y
        self.k = k
        self.r = 1/abs(k)

    def draw(self, ax):
        circle = plt.Circle((self.x, self.y), self.r, fill=False)
        ax.add_patch(circle)

    def is_tangent(self, other):
        return abs((self.x - other.x)**2 + (self.y - other.y)**2 - (sign(self.k) * self.r + sign(other.k) * other.r)**2) < TOL
        
class Gasket:
    def __init__(self, depth):
        outer_exists = False
        while not outer_exists:
            A, B, C = self.rand_pts()
            self.circles = self.get_soddy(A, B, C)

            for circle in self.get_circles(*self.circles):
                if circle.k < 0:
                    outer_exists = True
                    break 
        self.iterate(depth)

    def draw(self):
        c0 = next(filter(lambda c: c.k < 0, self.circles))

        fig, ax = plt.subplots()
        ax.set_xlim(c0.x - c0.r, c0.x + c0.r)
        ax.set_ylim(c0.y - c0.r, c0.y + c0.r)
        ax.axis('off')

        for circle in self.circles:
            circle.draw(ax)
        return fig

    def iterate(self, depth, circles=None, idx=0):
        if circles is None:
            circles = self.circles[:]

        if idx < depth:
            new_circles = self.get_circles(*circles)
            self.circles.extend(new_circles)

            for circle_pair in combinations(circles, 2):
                for new_circle in new_circles:
                    self.iterate(depth, [*circle_pair, new_circle], idx+1)
    @staticmethod
    def rand_pts():
        colinear = True
        while colinear:
            A = array([random(), random()])
            B = array([random(), random()])
            C = array([random(), random()])
            colinear = det(array([
                [A[0], B[0], C[0]],
                [A[1], B[1], C[1]],
                [   1,    1,    1]
            ])) == 0
        return A, B, C
    
    @staticmethod
    def get_soddy(A, B, C):
        a = norm(B - C)
        b = norm(C - A)
        c = norm(A - B)

        p = (a + b + c)/2

        return [
            Circle(*A, 1/(p - a)),
            Circle(*B, 1/(p - b)),
            Circle(*C, 1/(p - c))
        ]

    @staticmethod
    def get_circles(c1, c2, c3):
        rhsPos = lambda a, b, c: p1(a, b, c) + 2 * sqrt(p2(a, b, c))
        rhsNeg = lambda a, b, c: p1(a, b, c) - 2 * sqrt(p2(a, b, c))

        # Descarte (1643)
        kPos = rhsPos(c1.k, c2.k, c3.k)
        kNeg = rhsNeg(c1.k, c2.k, c3.k)

        # Wilks et al. (2002)
        zPos = rhsPos(c1.k*(c1.x + c1.y*1j), c2.k*(c2.x + c2.y*1j), c3.k*(c3.x + c3.y*1j))
        zNeg = rhsNeg(c1.k*(c1.x + c1.y*1j), c2.k*(c2.x + c2.y*1j), c3.k*(c3.x + c3.y*1j))

        new_circles = [
            Circle(zPos.real/kPos, zPos.imag/kPos, kPos),
            Circle(zPos.real/kNeg, zPos.imag/kNeg, kNeg),
            Circle(zNeg.real/kPos, zNeg.imag/kPos, kPos),
            Circle(zNeg.real/kNeg, zNeg.imag/kNeg, kNeg)
        ]

        for circle in [c1, c2, c2]:
                for new_circle in new_circles:
                    if not circle.is_tangent(new_circle):
                        new_circles.remove(new_circle)
        return new_circles
    
if __name__ == '__main__':
    start = time()

    depth = 6
    gasket = Gasket(depth)
    if len(gasket.circles) < 3**depth + 2:
        print('Warning: Not enough circles generated. Try decreasing tolerance.')
    fig = gasket.draw()
    fig.savefig('gasket.png')

    end = time()
    print(f'Generated {len(gasket.circles)} circles in {end - start:.2f} seconds')