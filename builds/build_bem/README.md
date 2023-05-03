## main.cpp

[main.cpp#L1](main.cpp#L1):

# BEM Simulation Code

This is a C++ implementation of a BEM simulation code. Follow the instructions below to build and run the simulation.

## Prerequisites

- CMake 3.26 or higher
- AppleClang 14.0.3 or compatible C++ compiler
- LAPACK library
- Eigen 3.4.0 or higher
- Python 3 for input generation

## Building the Code

1. Clean the build directory:

```
sh clean
```

2. Configure the build using CMake:

```
cmake -DCMAKE_BUILD_TYPE=Release ../
```

3. Compile the code:

```
make
```

## Running the Simulation

1. Generate input files using the `input_generator.py` script:

```
python3 ./input_generator.py
```

2. Run the simulation with the generated input files:

```
./main ./input_files/Kramer2021_H00d03
```

## Output

The simulation results will be stored in the specified output directory.

# settingBEM.py

プログラム内でつかわfれるパラメターや，入力値や出力先は`settingBEM.py`を実行することで作られる`json`ファイルで設定される．

**💡 NOTE:** `settingBEM.py`は`settingBEM.py`と同じフォルダ内にある必要がある．

# RKのtime step毎に，Dirichlet点にはΦを与える．Neumann点にはΦnを与える

どのように境界条件を適用するか．

# remesh（再配置）の条件

## flip,divide,mergeに共通する条件

辺のフリップ，分割，削除が実行されるには，辺で繋がる２点の境界条件が同じである必要がある．

(!((p0->Neumann && p1->Dirichlet) || (p0->Dirichlet && p1->Neumann)))

がtrueである場合のみ，辺の修正を実行することができる．

## flipの条件

## divideの条件

## mergeの条件


![](https://github.com/tomoakihirakawa/cpp/blob/main/builds/build_bem/anim.gif)

![](WATCHME_settingjson.mov)

![](WATCHME_settingBEM.mov)

