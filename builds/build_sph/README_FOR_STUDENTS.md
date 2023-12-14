# 基本的な使い方

#### ターミナルを開く

<table>
    <tr>
        <td>
            <figure>
                <img src="./img/README_FOR_STUDENTS_terminal1.png" width="500px" alt="Image Description">
                <figcaption>ターミナルを開いてファイルのあるディレクトリまで移動する</figcaption>
            </figure>
        </td>
        <td>
            <figure>
                <img src="./img/README_FOR_STUDENTS_terminal2.png" width="500px" alt="Image Description">
                <figcaption>または，vscodeを開き，Ctrl + jでターミナルを開き，ファイルのあるディレクトリまで移動する</figcaption>
            </figure>
        </td>
    </tr>
</table>


#### `main.cpp`が保存されているディレクトリに移動する．

```shell
cd ~/code/cpp/builds/build_sph　#　ビルドディレクトリに移動
code ./workspace # vscodeで開く
```

#### 実行する際に読み込ませる入力ファイルを作成する．

`input_generator.py`を実行すると，`input.txt`が生成される．

<img src="./img/README_FOR_STUDENTS_input_generator.png" width="500px"> 

この場合，最後に表示されている`./input_files/Lobovsky2013_PS0d018_CSML2d5_RK1`がインプットファイルのパスになる．

#### `main.cpp`を修正する
 
例えば次の箇所を修正することで，ISPHとEISPHを切り替えることができる．

```cpp
#define USE_ISPH //->#define USE_EISPH
```

#### コンパイルして実行する．

コンパイルには`cmake`を使用している．`cmake`は，`CMakeLists.txt`に書かれた内容に従って，ヘッダファイルやライブラリを探して，コンパイルを行い，実行ファイルを生成する．

```shell
sh clean # 古いCMakeFilesなどを削除する．
cmake -DCMAKE_BUILD_TYPE=Release ../　# Releaseモードでコンパイルする．
make　# コンパイル（ビルド）する．
```

実行する．

```shell
./main ./input_files/Lobovsky2013_PS0d018_CSML2d5_RK1 # 実行する．
```

その他利用するコマンド：

* 強制終了: Ctrl + C
