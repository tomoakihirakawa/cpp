# Contents

- [🐋pybind11の使い方](#🐋pybind11の使い方)
    - [⛵️Lighthill Robot](#⛵️Lighthill-Robot)
        - [🪸コンパイル方法](#🪸コンパイル方法)


---
# 🐋pybind11の使い方 

## ⛵️Lighthill Robot 

ラズパイでサーボモーターを動かすには，pythonを使うのが簡単．
ただ，数値計算においては，pythonの速度が遅いため実用的でなくなる場合があり，それを考慮しながらやっていくことは面倒．
そこで，pybind11を使って，pythonからでも読み込める共有ライブラリをc++を元に作る．

### 🪸コンパイル方法 

この例は，c++のNewton法を利用して作った[Lighthill Robot](../../include/rootFinding.hpp#L214)をpythonで使うためのもの.

以下をターミナルで実行して`make`すると，Macだと`LighthillRobot_pybind.cpython-311-darwin.so`が作られる.

```
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./ -DINPUT=LighthillRobot.cpp -DOUTPUT=shared_file_name_that_will_be_generated
make
```

[このように](../../builds/build_pybind11/runLightHillRobot.py#L15)`import`して利用できる．


[./LighthillRobot.cpp#L1](./LighthillRobot.cpp#L1)


⚠️ `cmake`の`-DOUTPUT`オプションで指定した名前と同じ`shared_file_name_that_will_be_generated`を指定する．

```
PYBIND11_MODULE(shared_file_name_that_will_be_generated, m) {
py::class_<class_name_declared_in_cpp>(m, "class_name_read_from_python")
...
```


[./LighthillRobot.cpp#L36](./LighthillRobot.cpp#L36)


---
出力結果

|水族館の魚|ロボット|
|:---:|:---:|
| <img src="sample_aquarium.gif"  width="80%" height="80%"> | ![sample.gif](sample.gif) |


[./runLightHillRobot.py#L1](./runLightHillRobot.py#L1)


---
