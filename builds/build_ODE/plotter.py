import csv
import matplotlib.pyplot as plt

# CSVファイルのパス
# node points are stored in result.csv
file_path = 'result.csv'

# データを読み込む
x_data = []
y_data = []

time_x_data = []
time_y_data = []

with open(file_path, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)  # ヘッダーを読み飛ばす
    i = 0
    for row in reader:
        for i in range(7):
            x_data.append(float(row[i]))
            y_data.append(i)

# プロット
plt.figure(figsize=(10, 6))
plt.plot(x_data, y_data, marker='o', linestyle='-', label=f'{header[0]} vs {header[1]}')

# グラフの装飾
plt.title('CSV Data Plot', fontsize=16)
plt.xlabel(header[0], fontsize=12)
plt.ylabel(header[1], fontsize=12)
plt.legend()
plt.grid(True)

# 表示
plt.show()