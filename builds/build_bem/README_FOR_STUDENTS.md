# 基本的な使い方

#### ターミナルを開く

<table>
    <tr>
        <td>
            <figure>
                <img src="./img/README_FOR_STUDENTS_terminal1.png" width="500px" alt="Image Description"><br>
                <figcaption>ターミナルを開いてファイルのあるディレクトリまで移動する</figcaption>
            </figure>
        </td>
        <td>
            <figure>
                <img src="./img/README_FOR_STUDENTS_terminal2.png" width="500px" alt="Image Description"><br>
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


# 比較対象

\cite{Lobovsky2014}と計算結果を比較しよう．

<img src="./img/README_FOR_STUDENTS_comparison_Lobovsky2014.png" width="500px">

paravisで可視化し，壁面に設置された圧力センサーの実験結果とSPHの計算結果を比較する．

* 計算方法を理解し，結果を考察する
* タイムステップや粒子数，平滑化距離などのパラメータを変更し，結果を比較する．
    * タイムステップは，CFL条件によって決まるので，CFL条件の設定を確認する．
    * 粒子数は，初期の粒子間距離によって決まるので，初期粒子間距離の設定を確認する．
    * 平滑化距離の変更に応じて，適切な核関数$`W`$も変化するので，核関数の設定を確認する．