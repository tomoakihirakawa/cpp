/*DOC_EXTRACT pybind11

# pybind11の使い方

## Lighthill Robot

ラズパイでサーボモーターを動かすには，pythonを使うのが簡単．
ただ，数値計算においては，pythonの速度が遅いため実用的でなくなる場合があり，それを考慮しながらやっていくことは面倒．
そこで，pybind11を使って，pythonからでも読み込める共有ライブラリをc++を元に作る．

### コンパイル方法

この例は，c++のNewton法を利用して作った\ref{newton:LighthillRobot}{Lighthill Robot}をpythonで使うためのもの.

以下をターミナルで実行して`make`すると，Macだと`LighthillRobot_pybind.cpython-311-darwin.so`が作られる.

```
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ./ -DINPUT=LighthillRobot.cpp -DOUTPUT=shared_file_name_that_will_be_generated
make
```

\ref{PYBIND11:HOW_TO_IMPORT}{このように}`import`して利用できる．

*/

#define NOMINMAX
#define _CRT_SECURE_NO_WARNINGS

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../include/minMaxOfFunctions.hpp"
#include "../../include/rootFinding.hpp"

namespace py = pybind11;

/*DOC_EXTRACT pybind11

WARNING `cmake`の`-DOUTPUT`オプションで指定した名前と同じ`shared_file_name_that_will_be_generated`を指定する．

```
PYBIND11_MODULE(shared_file_name_that_will_be_generated, m) {
   py::class_<class_name_declared_in_cpp>(m, "class_name_read_from_python")
...
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
