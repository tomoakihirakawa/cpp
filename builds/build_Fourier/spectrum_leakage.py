import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
from scipy.signal import get_window

# サンプリング周波数と信号の長さ
fs = 1000
N = 10000
t = np.linspace(0, 1, N, endpoint=False)

# 複数の周波数成分を持つ信号を生成
np.random.seed(0)  # 再現性のため
frequencies = [50, 80, 120, 200]
original_signal = sum(np.sin(2 * np.pi * f * t) for f in frequencies)

# ガウシアンノイズを追加
noisy_signal = original_signal + np.random.normal(size=t.shape)

# 窓関数の選択
window = get_window('hann', N)

# 窓関数を適用してフーリエ変換
windowed_signal = noisy_signal * window
fft_noisy_signal = fft(noisy_signal)
fft_windowed_signal = fft(windowed_signal)

# 高周波成分を除去
cutoff_frequency = 210  # 適切なカットオフ周波数を選択
fft_noisy_signal_filtered = np.copy(fft_noisy_signal)
fft_windowed_signal_filtered = np.copy(fft_windowed_signal)

fft_noisy_signal_filtered[cutoff_frequency:] = 0
fft_noisy_signal_filtered[-cutoff_frequency:] = 0
fft_windowed_signal_filtered[cutoff_frequency:] = 0
fft_windowed_signal_filtered[-cutoff_frequency:] = 0

# 逆フーリエ変換
recovered_signal_no_window = ifft(fft_noisy_signal_filtered).real
recovered_signal_with_window = ifft(fft_windowed_signal_filtered).real

# 結果をプロット
freqs = np.fft.fftfreq(N, 1/fs)

plt.figure(figsize=(15, 12))

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

# 窓関数を使わずに復元した信号
plt.subplot(4, 1, 3)
plt.plot(t, recovered_signal_no_window, label='Recovered Signal without Window')
plt.title('Recovered Signal without Window')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

# 窓関数を使って復元した信号
plt.subplot(4, 1, 4)
plt.plot(t, recovered_signal_with_window, label='Recovered Signal with Window')
plt.title('Recovered Signal with Window')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
