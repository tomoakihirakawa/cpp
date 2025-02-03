# import numpy as np

# # 固定パラメータ
# E = 205e9  # ヤング率 (Pa)
# b = 0.035  # 幅 (m)
# h = 0.002  # 高さ (m)
# rho = 7.8e3  # 密度 (kg/m^3)
# A = b * h  # 断面積 (m^2)
# I = b * h**3 / 12  # 断面二次モーメント (m^4)
# L = 0.9  # 梁の長さ (m)

# # 固有振動数を計算
# def calc_natural_frequencies(num_modes):
#     """
#     固有振動数（角周波数）を計算する関数

#     Parameters:
#         num_modes (int): 計算するモード数

#     Returns:
#         omega (numpy.ndarray): 固有振動数（角周波数）
#     """
#     # 固有振動モードの波数 k
#     k = np.array([(i + 1) * np.pi / L for i in range(num_modes)])  # 固定端と自由端の境界条件に基づく
#     # 固有振動数（角周波数）を計算
#     omega = np.sqrt((E * I * k**4) / (rho * A))
#     return omega

# # モード形状
# def mode_shape(x, n):
#     """
#     固有モード形状を計算する関数

#     Parameters:
#         x (float): x座標
#         n (int): モード番号

#     Returns:
#         (float): モード形状
#     """
#     return np.sin((n + 1) * np.pi * x / L)

# # 振幅の計算
# def calc_amplitudes(num_modes, initial_displacement_func):
#     """
#     初期条件に基づく各モードの振幅を計算する関数

#     Parameters:
#         num_modes (int): 計算するモード数
#         initial_displacement_func (function): 初期変位を返す関数

#     Returns:
#         amplitudes (numpy.ndarray): 各モードの振幅
#     """
#     x_values = np.linspace(0, L, 1000)  # 空間分割
#     amplitudes = []
#     for n in range(num_modes):
#         phi_n = mode_shape(x_values, n)
#         u0 = initial_displacement_func(x_values)
#         # 振幅はモード形状と初期変位の内積で決まる
#         A_n = np.trapz(rho * A * phi_n * u0, x_values) / np.trapz(rho * A * phi_n**2, x_values)
#         amplitudes.append(A_n)
#     return np.array(amplitudes)

# # 初期変位の定義
# def initial_displacement(x):
#     return 0.000001 * x / L  # 例: 正弦波の初期変位
#     # return 0.01 * np.sin(np.pi * x / L)  # 例: 正弦波の初期変位

# # 計算モード数を指定
# num_modes = 20
# freqs = calc_natural_frequencies(num_modes)
# amplitudes = calc_amplitudes(num_modes, initial_displacement)

# # 結果をフォーマット
# freqs_rad_s = ", ".join([f"{val:.6g}" for val in freqs])
# freqs_hz = ", ".join([f"{val:.6g}" for val in freqs / (2 * np.pi)])
# amplitudes_str = ", ".join([f"{val:.6g}" for val in amplitudes])

# # 結果を表示
# print(f"固有振動数（角周波数 rad/s）: {{{freqs_rad_s}}}")
# print(f"固有振動数（Hz）: {{{freqs_hz}}}")
# print(f"モード振幅: {{{amplitudes_str}}}")


import numpy as np
from scipy.optimize import root_scalar

# 固定パラメータ
E = 205e9  # ヤング率 (Pa)
b = 0.035  # 幅 (m)
h = 0.002  # 高さ (m)
rho = 7.8e3  # 密度 (kg/m^3)
A = b * h  # 断面積 (m^2)
I = b * h**3 / 12  # 断面二次モーメント (m^4)
L = 0.9  # 梁の長さ (m)

# 固有振動数を計算
def characteristic_equation(beta):
    return np.cos(beta) * np.cosh(beta) + 1

def solve_beta_roots(num_modes):
    roots = []
    for n in range(1, num_modes + 1):
        guess = (2 * n - 1) * np.pi / 2
        sol = root_scalar(characteristic_equation, bracket=[guess - np.pi / 4, guess + np.pi / 4])
        if sol.converged:
            roots.append(sol.root)
    return np.array(roots)

def calc_natural_frequencies_continuous(num_modes):
    beta_n = solve_beta_roots(num_modes)
    omega_n = (beta_n**2) * np.sqrt(E * I / (rho * A))
    return omega_n, beta_n

def mode_shape_continuous(x, beta_n, L):
    return np.cosh(beta_n * x / L) - np.cos(beta_n * x / L) - \
           (np.sinh(beta_n) - np.sin(beta_n)) / (np.cosh(beta_n) + np.cos(beta_n)) * \
           (np.sinh(beta_n * x / L) - np.sin(beta_n * x / L))

# 計算モード数を指定
num_modes = 25
freqs_rad_s, beta_n = calc_natural_frequencies_continuous(num_modes)
freqs_hz = freqs_rad_s / (2 * np.pi)

# モード形状を計算
x_vals = np.linspace(0, L, 100)
mode_shapes = [mode_shape_continuous(x_vals, beta, L) for beta in beta_n]

# Mathematica向けにフォーマット
freqs_hz_list = "{" + ", ".join([f"{val:.6f}" for val in freqs_hz]) + "}"
freqs_rad_s_list = "{" + ", ".join([f"{val:.6f}" for val in freqs_rad_s]) + "}"
beta_n_list = "{" + ", ".join([f"{val:.6f}" for val in beta_n]) + "}"
mode_shapes_list = "{" + ", ".join(["{" + ", ".join([f"{val:.6f}" for val in shape]) + "}" for shape in mode_shapes]) + "}"
x_vals_list = "{" + ", ".join([f"{val:.6f}" for val in x_vals]) + "}"

# 結果を表示
print("(* Mathematica 用の出力 *)")
print(f"自然周波数 (Hz): freqHz = {freqs_hz_list};")
print(f"自然周波数 (rad/s): freqRadS = {freqs_rad_s_list};")
print(f"β の値: betaVals = {beta_n_list};")
print(f"x 座標: xVals = {x_vals_list};")
print(f"モード形状: modeShapes = {mode_shapes_list};")