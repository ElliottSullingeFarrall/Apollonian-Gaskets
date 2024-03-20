from itertools import combinations
from random import uniform

import matplotlib.pyplot as plt
from mpmath import mp
TOL = 1e-6

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

    def draw(self, ax, label=None):
        circle = plt.Circle((self.x, self.y), self.r, fill=False)
        ax.annotate(label, xy=(self.x, self.y), fontsize=6, ha="center")
        ax.add_patch(circle)

    def is_tangent(self, other):
        return mp.almosteq((self.x - other.x)**2 + (self.y - other.y)**2, (mp.sign(self.k) * self.r + mp.sign(other.k) * other.r)**2, rel_eps=TOL)
        

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
        ax.set_xlim(float(c0.x - c0.r), float(c0.x + c0.r))
        ax.set_ylim(float(c0.y - c0.r), float(c0.y + c0.r))
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
            A = mp.matrix([uniform(-1, +1), uniform(-1, +1)])
            B = mp.matrix([uniform(-1, +1), uniform(-1, +1)])
            C = mp.matrix([uniform(-1, +1), uniform(-1, +1)])
            colinear = mp.det(mp.matrix([
                [A[0], B[0], C[0]],
                [A[1], B[1], C[1]],
                [   1,    1,    1]
            ])) == 0
        return A, B, C
    
    @staticmethod
    def get_soddy(A, B, C):
        a = mp.norm(B - C)
        b = mp.norm(C - A)
        c = mp.norm(A - B)

        p = (a + b + c)/2

        return [
            Circle(*A, 1/(p - a)),
            Circle(*B, 1/(p - b)),
            Circle(*C, 1/(p - c))
        ]

    @staticmethod
    def get_circles(c1, c2, c3):
        rhsPos = lambda a, b, c: p1(a, b, c) + 2 * mp.sqrt(p2(a, b, c))
        rhsNeg = lambda a, b, c: p1(a, b, c) - 2 * mp.sqrt(p2(a, b, c))

        # Descarte (1643)
        kPos = rhsPos(c1.k, c2.k, c3.k).real
        kNeg = rhsNeg(c1.k, c2.k, c3.k).real

        # Wilks et al. (2002)
        zPos = rhsPos(c1.k*mp.mpc(c1.x, c1.y), c2.k*mp.mpc(c2.x, c2.y), c3.k*mp.mpc(c3.x, c3.y))
        zNeg = rhsNeg(c1.k*mp.mpc(c1.x, c1.y), c2.k*mp.mpc(c2.x, c2.y), c3.k*mp.mpc(c3.x, c3.y))

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
    gasket = Gasket(5)
    fig = gasket.draw()
    fig.savefig('gasket.png')