
import math

def f(x):
    return math.cos(x)

def ans(a,b):
    return math.sin(b) - math.sin(a)

def trapezoidal_rule(f, a, b, N):
    h = (b-a)/N
    return (h/2) * (f(a) + 2*sum(f(a + h*i) for i in range(1, N)) + f(b))

print("周期関数の場合，0~math.piまでの台計則による数値積分はものすごく精度が高いので，Nが小さくてとてもいい結果が出ます．")
a = 0
b = math.pi
for n in range(2, 20):
    result = trapezoidal_rule(f, a, b, n) - ans(a,b)
    print(f'n={n} 差 : {result}')

print("周期関数の場合，0~4までの台計則による数値積分")
print("nが大きくなると精度がよくります")
a = 0
b = 4.

A = []

for n in range(2, 20):
    result = trapezoidal_rule(f, a, b, n) - ans(a,b)
    A.append(result)
    print(f'n={n} 差 : {result}')


import matplotlib.pyplot as plt

plt.plot(A)
plt.xlabel('n')
plt.ylabel('Difference')
plt.title('Difference between Trapezoidal Rule and Actual Integral')
plt.show()


