# 39	Default foreground color
# echo -e "\e[39mDefault"
# Default Default
# 30	Black
# echo -e "Default \e[30mBlack"
# Default Black
# 31	Red
# echo -e "Default \e[31mRed"
# Default Red
# 32	Green
# echo -e "Default \e[32mGreen"
# Default Green
# 33	Yellow
# echo -e "Default \e[33mYellow"
# Default Yellow
# 34	Blue
# echo -e "Default \e[34mBlue"
# Default Blue
# 35	Magenta
# echo -e "Default \e[35mMagenta"
# Default Magenta
# 36	Cyan
# echo -e "Default \e[36mCyan"
# Default Cyan
# 37	Light gray
# echo -e "Default \e[37mLight gray"
# Default Light gray
# 90	Dark gray
# echo -e "Default \e[90mDark gray"
# Default Dark gray
# 91	Light red
# echo -e "Default \e[91mLight red"
# Default Light red
# 92	Light green
# echo -e "Default \e[92mLight green"
# Default Light green
# 93	Light yellow
# echo -e "Default \e[93mLight yellow"
# Default Light yellow
# 94	Light blue
# echo -e "Default \e[94mLight blue"
# Default Light blue
# 95	Light magenta
# echo -e "Default \e[95mLight magenta"
# Default Light magenta
# 96	Light cyan
# echo -e "Default \e[96mLight cyan"
# Default Light cyan
# 97	White
# echo -e "Default \e[97mWhite"
from math import *
from fundamental import *
from time import sleep

red = "\033[31m"
blue = "\033[34m"
default = "\033[39m"

'''
回転軸と回転角度を与えた場合に，正しく回転を計算できるかチェック
'''

print("------------- チェック -------------")
print("y,z平面のpi/4回転")
axis = (1, 0., 0.)
angle = pi/4.
q = Quaternion(axis, angle)
gravity = (0, 0, -sqrt(2))


'''
        A z
        |
        |
y <-----+
        |  gravity = (0, 0, -sqrt(2))
        |
        V
'''
print("   座標を回転", "(0,  1, -1)", "<->", "結果", red, q.Rs(gravity), default)
print("ベクトルを回転", "(0, -1, -1)", "<->", "結果", red, q.Rv(gravity), default)

print("-------------------------------")
'''
        A z
        |
        |
        +------>x
        |  gravity = (0, 0, -sqrt(2))
        |
        V
'''
print("x,z平面のpi/4回転")
axis = (0, 1, 0)
angle = pi/4.
q = Quaternion(axis, angle)
print("   座標を回転", "(-1, 0, -1)", "<->", "結果", red, q.Rs(gravity), default)
print("ベクトルを回転", "( 1, 0, -1)", "<->", "結果", red, q.Rv(gravity), default)

print("-------------------------------")
'''
        A y
        |
        |
        +------>x
X gravity = (0, 0, -sqrt(2))
'''
print("x,y平面のpi/4回転")
axis = (0, 0, 1)
angle = pi/4.
q = Quaternion(axis, angle)
print("   座標を回転", "(0, 0, -sqrt(2))", "<->", "結果", red, q.Rs(gravity), default)
print("ベクトルを回転", "(0, 0, -sqrt(2))", "<->", "結果", red, q.Rv(gravity), default)

print("-------------------------------")

'''
各軸の角速度をクォターニオンに変換
'''
print("各軸の角速度をクォターニオンに変換")

q = Quaternion(axis, 0.)
rad = pi/1000.
for i in range(1000):
    sleep(.1)
    w = (rad, 0., 0.)  # [rad/s]
    q = Add(q, q.d_dt(w))
    estimate_rotation = (i+1)*rad
    print("roll = (予測)", estimate_rotation/pi * 180., "<->",
          red, "（結果）", q.YPR()[2]/pi*180., default)
    p = Quaternion((1, 0, 0), estimate_rotation)
    print("   座標を回転", p.Rs(gravity), "<->",
          red, "結果", q.Rs(gravity), default)
    print("ベクトルを回転", p.Rv(gravity), "<->",
          red, "結果", q.Rv(gravity), default)
