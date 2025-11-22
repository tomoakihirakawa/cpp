import math

L = 1.
w = 2.*math.pi
k = 4.*math.pi
c1 = 0.2
c2 = 0.2


def yLH(x, t):
    return (c1 * x / L + c2 * (x / L)**2) * math.sin(k * (x / L) - w * t)


print("y = ", yLH(0.1, 0.1))
