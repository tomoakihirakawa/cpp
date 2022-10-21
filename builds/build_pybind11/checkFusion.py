
from math import *
from fundamental import *
from time import sleep, time_ns

red = "\033[31m"
blue = "\033[34m"
default = "\033[39m"

theta = -54/180*pi
G0 = (0, 0., -1.)
M0 = (cos(theta), 0., sin(theta))
fusion = Fusion(G0, M0)
# -------------------------------------------------------- #
s = time_ns()
axis = (1., 0., 0.)
angle = 70./180.*pi
q = Quaternion(axis, angle)

print("q.Rs(M0) = ", q.Rs(M0))
print("q.Rs(G0) = ", q.Rs(G0))

ans = fusion.solveForQuaternion(q.Rs(G0), q.Rs(M0), [0, 0, 0], 0.)
print(Norm(Subtract(q(), ans())))
# Rsでいいのか？
print("time elapsed", (time_ns()-s)*10**-9)

# -------------------------------------------------------- #
s = time_ns()
for i in range(100):
    sleep(1.)
    axis = (1., 0., 0.)
    angle = 10.*i/180.*pi
    q = Quaternion(axis, angle)
    t = time_ns()
    intp = fusion.interp((t-s)*10**-9)
    ans = fusion.solveForQuaternion(
        q.Rs(G0), q.Rs(M0), [0, 0, 0], (t-s)*10**-9)
    print(ans(), ",", intp[0]())
    # print("time elapsed", (time_ns()-s)*10**-9)
