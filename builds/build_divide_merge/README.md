# Contents

- [🐋 Fusion360を使って計算用objファイルを生成](#🐋-Fusion360を使って計算用objファイルを生成)
- [🐋 計算用にメッシュの細分化](#🐋-計算用にメッシュの細分化)
    - [⛵ 実行ファイルの作成方法（`remesh.cpp`のコンパイル方法）](#⛵-実行ファイルの作成方法（`remesh.cpp`のコンパイル方法）)
    - [⛵ 実行方法](#⛵-実行方法)


---
# 🐋 Fusion360を使って計算用objファイルを生成 

<img src="sample_fusion360_step1.png" width="600px">

<img src="sample_fusion360_step2.png" width="600px">

<img src="sample_fusion360_step3.png" width="600px">

<img src="sample_fusion360_step4.png" width="600px">

# 🐋 計算用にメッシュの細分化 

## ⛵ 実行ファイルの作成方法（`remesh.cpp`のコンパイル方法） 

1. `sh clean`で古いファイルを削除する．
2. `cmake`を使って，`CMakeLists.txt`から`Makefile`を生成する．（リリースビルドタイプを指定し，ソースファイルを`remesh.cpp`に設定）
3. `make`コマンドで，`Makefile`に基づいて実行ファイルをコンパイルする.

```shell
$ sh clean
$ cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp
$ make
```

`remesh`という実行ファイルができる．

## ⛵ 実行方法 

`n`回の細分化を行う．

```
./remesh input_file output_dir output_name n
```

![./sample2.gif](sample2.gif)

出力は，`output_dir/output_name*.vtu`と`output_dir/output_name*.obj`．


![./sample.gif](sample.gif)

次は，入力ファイルを生成し，計算をする．


[./remesh.cpp#L3](./remesh.cpp#L3)


---
