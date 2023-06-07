# Contents

- [🐋pybind11の例](#🐋pybind11の例)


---
# 🐋pybind11の例 

この例は，c++のNewton法を利用して作ったLight Hill Robotをpythonで使うためのもの.
以下をターミナルで実行して`make`すると，Macだと`LightHillRobot_pybind.cpython-311-darwin.so`が作られる.

```
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./
make
```

`run_robot.py`にあるように`import`して利用できる．


[./LightHillRobot_pybind.cpp#L1](./LightHillRobot_pybind.cpp#L1)


---
出力結果

![sample.gif](sample.gif)


[./run_robot.py#L1](./run_robot.py#L1)


---
