import numpy as np
import matplotlib.pyplot as plt
import pywt

# サンプリング周波数と信号の長さ
fs = 1000
N = 10000
t = np.linspace(0, 1, N, endpoint=False)

# 複数の周波数成分を持つ信号を生成
np.random.seed(0)  # 再現性のため
frequencies = [50, 80, 120, 200]
original_signal = sum(np.sin(2 * np.pi * f * t) for f in frequencies)

# ガウシアンノイズを追加
noisy_signal = original_signal + np.random.normal(size=t.shape)*np.random.normal(size=t.shape)

# ウェーブレット変換
wavelet = 'db4'
coeffs = pywt.wavedec(noisy_signal, wavelet, level=10)

# 高周波成分を除去
threshold = 0.8 * np.std(coeffs[-1])
for i in range(1, len(coeffs)):
    coeffs[i] = pywt.threshold(coeffs[i], threshold)

# 逆ウェーブレット変換
reconstructed_signal = pywt.waverec(coeffs, wavelet)

# 結果をプロット
plt.figure(figsize=(15, 10))

# 元の信号（ノイズなし）
plt.subplot(4, 1, 1)
plt.plot(t, original_signal, label='Original Signal (No Noise)')
plt.title('Original Signal (No Noise)')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

# ノイズの乗った信号のプロット
plt.subplot(4, 1, 2)
plt.plot(t, noisy_signal, label='Original Signal with Noise')
plt.title('Original Signal with Noise')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

# ウェーブレット変換を使って復元した信号
plt.subplot(4, 1, 3)
plt.plot(t, reconstructed_signal, label='Recovered Signal with Wavelet Transform')
plt.title('Recovered Signal with Wavelet Transform')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

# 元の信号と復元信号の比較
plt.subplot(4, 1, 4)
plt.plot(t, original_signal, label='Original Signal (No Noise)', alpha=0.5)
plt.plot(t, reconstructed_signal, label='Recovered Signal with Wavelet Transform', alpha=0.7)
plt.title('Comparison of Original and Recovered Signals')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
