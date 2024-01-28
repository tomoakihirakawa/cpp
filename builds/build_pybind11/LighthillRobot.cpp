/*DOC_EXTRACT 0_1_0_how_to_make_shared_library

## pybind11で共有ライブラリを作成

この例は，c++のNewton法を利用して作った\ref{newton:LighthillRobot}{Lighthill Robot}をpythonで使うためのもの.

このディレクトリにCMakelists.txtを用意しているので，
それを使って，以下のようにターミナル上で実行・`make`すると，
Macだと`LighthillRobot_pybind.cpython-311-darwin.so`が作られる.

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./ -DINPUT=LighthillRobot.cpp -DOUTPUT=shared_file_name_that_will_be_generated
make
```

具体的には，以下のようにコンパイルする．

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./ -DINPUT=LighthillRobot.cpp -DOUTPUT=LighthillRobot -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-13
make
```

### ラズパイでのコンパイル

c++で-std=c++17を使うためには，gcc-9.1.0以上が必要．
[https://gist.github.com/sol-prog/95e4e7e3674ac819179acf33172de8a9#file-commands-sh](https://gist.github.com/sol-prog/95e4e7e3674ac819179acf33172de8a9#file-commands-sh)
ここを参考にして，まずはgcc-9.1.0をインストールする．
git cloneをする際にプロキシを通す必要がある場合は，以下のようにする．

```sh
# Commands used in the video https://youtu.be/-bCG87jBDqA :

sudo apt update && sudo apt upgrade -y

git clone https://bitbucket.org/sol_prog/raspberry-pi-gcc-binary.git
cd raspberry-pi-gcc-binary
tar -xjvf gcc-9.1.0-armhf-raspbian.tar.bz2
sudo mv gcc-9.1.0 /opt
cd ..
rm -rf raspberry-pi-gcc-binary

cd ~
echo 'export PATH=/opt/gcc-9.1.0/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/opt/gcc-9.1.0/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
. ~/.bashrc
sudo ln -s /usr/include/arm-linux-gnueabihf/sys /usr/include/sys
sudo ln -s /usr/include/arm-linux-gnueabihf/bits /usr/include/bits
sudo ln -s /usr/include/arm-linux-gnueabihf/gnu /usr/include/gnu
sudo ln -s /usr/include/arm-linux-gnueabihf/asm /usr/include/asm
sudo ln -s /usr/lib/arm-linux-gnueabihf/crti.o /usr/lib/crti.o
sudo ln -s /usr/lib/arm-linux-gnueabihf/crt1.o /usr/lib/crt1.o
sudo ln -s /usr/lib/arm-linux-gnueabihf/crtn.o /usr/lib/crtn.o

g++-9.1 -std=c++17 -Wall -pedantic test_fs.cpp -o test_fs
./test_fs
```

```sh
git config --global http.proxy http://書き換え:8080
git config --global https.proxy http://書き換え:8080
```



*/
#define NOMINMAX
#define _CRT_SECURE_NO_WARNINGS

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include "../../include/minMaxOfFunctions.hpp"
#include "../../include/rootFinding.hpp"

namespace py = pybind11;

/*DOC_EXTRACT 0_0_0_how_to_write_for_pybind11

# pybind11の使い方

ラズパイでサーボモーターを動かすには，pythonを使うのが簡単．
ただ，数値計算においては，pythonの速度が遅いため実用的でなくなる場合があり，それを考慮しながらやっていくことは面倒．
そこで，pybind11を使って，pythonからでも読み込める共有ライブラリをc++を元に作った．

## pybind11の書き方

WARNING `cmake`の`-DOUTPUT`オプションで指定した名前と同じ`shared_file_name_that_will_be_generated`を指定する．

```cpp
PYBIND11_MODULE(LighthillRobot, m) {
   py::class_<LighthillRobot>(m, "LighthillRobot")
       .def(py::init<double, double, double, double, double, int>())
       .def_readwrite("c1", &LighthillRobot::c1)
       .def_readwrite("c2", &LighthillRobot::c2)
       .def("yLH", &LighthillRobot::yLH)
       .def("X_RB", &LighthillRobot::X_RB)
       .def("f", &LighthillRobot::f)
       .def("ddx_yLH", &LighthillRobot::ddx_yLH)
       .def("ddq_f", &LighthillRobot::ddq_f)
       .def("getAngles", &LighthillRobot::getAngles)
       .def("anglesToX", &LighthillRobot::anglesToX);
}
```

*/

PYBIND11_MODULE(LighthillRobot, m) {
   py::class_<LighthillRobot>(m, "LighthillRobot")
       .def(py::init<double, double, double, double, double, int>())
       .def_readwrite("c1", &LighthillRobot::c1)
       .def_readwrite("c2", &LighthillRobot::c2)
       .def("yLH", &LighthillRobot::yLH)
       .def("X_RB", &LighthillRobot::X_RB)
       .def("f", &LighthillRobot::f)
       .def("ddx_yLH", &LighthillRobot::ddx_yLH)
       .def("ddq_f", &LighthillRobot::ddq_f)
       .def("getAngles", &LighthillRobot::getAngles)
       .def("anglesToX", &LighthillRobot::anglesToX);
}
