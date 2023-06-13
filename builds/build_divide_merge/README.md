# Contents

- [🐋メッシュの細分化](#🐋メッシュの細分化)
    - [⛵️実行ファイルの作成方法（`remesh.cpp`のコンパイル方法）](#⛵️実行ファイルの作成方法（`remesh.cpp`のコンパイル方法）)
    - [⛵️実行方法](#⛵️実行方法)


---
# 🐋メッシュの細分化 

## ⛵️実行ファイルの作成方法（`remesh.cpp`のコンパイル方法） 

古いファイルを削除する：

```
sh clean
```

`remesh.cpp`を`cmake`を使ってコンパイルするために，`Makefile`を作る：

```
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp
```

`Makefile`を使ってコンパイルする：

```
make
```

`remesh`という実行ファイルができます．

## ⛵️実行方法 

`n`回の細分化を行う．

```
./remesh input_file output_dir output_name n
```

![./sample2.gif](sample2.gif)

出力は，`output_dir/output_name*.vtu`と`output_dir/output_name*.obj`．


![./sample.gif](sample.gif)


<div align="right>
<p align="right">
<a href="./remesh.cpp#L3">./remesh.cpp#L3</a>
</p>
</div>


---
