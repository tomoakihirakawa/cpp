# Contents


---
実行方法：

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

オイラー法，Leap-Frog法，Runge-Kutta法
を用いて，弾性体の動きをシミュレーション．

`const double stiffness = 1000;`の場合
![sample.gif](sample.gif)

`const double stiffness = 10000;`の場合
![sample_2.gif](sample_2.gif)

[./main.cpp#L11](./main.cpp#L11)

---
