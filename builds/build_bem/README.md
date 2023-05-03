## main.cpp

[main.cpp#L1](main.cpp#L1):

# BEM Simulation Code

This is a C++ implementation of a BEM simulation code. Follow the instructions below to build and run the simulation.

## Prerequisites

- CMake
- LAPACK library
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


![](https://github.com/tomoakihirakawa/cpp/blob/main/builds/build_bem/anim.gif)

![](WATCHME_settingjson.mov)

![](WATCHME_settingBEM.mov)

## input_generator.py

[input_generator.py#L1](input_generator.py#L1):

# Input Generator for BEM Simulation

This Python script generates input files for the BEM simulation code. It supports various simulation cases and handles input file generation for each case.

## Prerequisites

- Python 3

## Usage

1. Make sure the required dependencies are installed.
2. Run the script using the following command:

```
python3 input_generator.py

```

Upon running the script, it will generate input files in JSON format for the specified simulation case. The input files are saved in the `./input_files/` directory.

## Customization

To customize the input file generation for a specific case, follow these steps:

1. Locate the `SimulationCase` variable in the script and set it to the desired case name, e.g., `"Kramer2021"`.
2. Add a new `case` block in the `match SimulationCase:` section to handle the new simulation case.
3. Define the required parameters for the simulation case within the new `case` block, following the examples provided in the script.
4. Update the `inputfiles` variable with the new input objects created for the custom case.

After customizing the script, run it again to generate the input files for the new case.

## Output

The script will generate input files in JSON format for the specified simulation case. The input files will be saved in the `./input_files/` directory. The generated input files can be used to run the BEM simulation.

[input_generator.py#L55](input_generator.py#L55):

プログラムを回す際に面倒な事は，入力ファイルの設定．
入力ファイルの作り方をドキュメントで示されても，具体的な例がないとわかりにくい．
例があっても，例と違う場合どうすればいいかなど，わからないことは多い．
このように，入力ファイルを生成するプログラムを作っておけば，その面倒をだいぶ解消できる．

