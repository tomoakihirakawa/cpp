# Contents

- [🐋pybind11の使い方](#🐋pybind11の使い方)
    - [⛵️Light Hill Robot](#⛵️Light-Hill-Robot)
        - [🪸コンパイル方法](#🪸コンパイル方法)


---
# 🐋pybind11の使い方 

## ⛵️Light Hill Robot 

### 🪸コンパイル方法 

この例は，c++のNewton法を利用して作ったLight Hill Robotをpythonで使うためのもの.
以下をターミナルで実行して`make`すると，Macだと`LightHillRobot_pybind.cpython-311-darwin.so`が作られる.

```
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./ -DINPUT=LightHillRobot.cpp -DOUTPUT=shared_file_name_that_will_be_generated
make
```

[このように](../../builds/build_pybind11/runLightHillRobot.py#L13)`import`して利用できる．


[./LightHillRobot.cpp#L1](./LightHillRobot.cpp#L1)


⚠️ `cmake`の`-DOUTPUT`オプションで指定した名前と同じ`shared_file_name_that_will_be_generated`を指定する．

```
PYBIND11_MODULE(shared_file_name_that_will_be_generated, m) {
py::class_<class_name_declared_in_cpp>(m, "class_name_read_from_python")
...
```


[./LightHillRobot.cpp#L98](./LightHillRobot.cpp#L98)


---
出力結果

![sample.gif](sample.gif)


[./runLightHillRobot.py#L1](./runLightHillRobot.py#L1)


---
