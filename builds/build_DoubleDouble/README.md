DoubleDoubleを作成する前に，数値がどのタイミングで丸められるかを知っておく必要がある．

```cpp
sh clean
cmake -DCMAKE_BUILD_TYPE=Debug ../ -DSOURCE_FILE=test_DoubleDouble.cpp
make
./test_DoubleDouble
```

この実装は間違いがある．

```cpp
DoubleDouble operator+(const DoubleDouble& other) const {
DoubleDouble result;
result.a = a + other.a;//!ここで丸められてしまい誤差が生じる
result.b = b + other.b;
return result;
}
```

なので，このdouble同士の演算においても誤差を取りこぼさない工夫が必要となる．

[:material-microsoft-visual-studio-code:test_DoubleDouble.cpp#L9](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_DoubleDouble/test_DoubleDouble.cpp:9)

DoubleDoubleを作成する前に，数値がどのタイミングで丸められるかを知っておく必要がある．

```cpp
sh clean
cmake -DCMAKE_BUILD_TYPE=Debug ../ -DSOURCE_FILE=test_DoubleDouble_correct.cpp
make
./test_DoubleDouble_correct
```

[:material-microsoft-visual-studio-code:test_DoubleDouble_correct.cpp#L11](vscode://file//Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_DoubleDouble/test_DoubleDouble_correct.cpp:11)

---
