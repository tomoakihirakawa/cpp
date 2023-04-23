[main.cpp#L1](main.cpp#L1):

# コンパイルのやり方

以下のコマンドの先頭の"$"は無視してください．
本来はコンパイルの際には，多くのヘッダーファイルをインクルードするよう長いコンパイルのコマンドを打つ必要がある．
cmakeを使えば，それをCMakeLists.txtにあらかじめ書いておくことで省くことができる．
```shell
$ cmake -DCMAKE_BUILD_TYPE=Release ../
```
次に，
```shell
$ make
```
これでコンパイル終了．後は，次のようにすればmainファイルが実行される．
```shell
$ ./main
```
ただし，古いcmake情報が今のフォルダ内に残っている場合，その情報を削除しておかないと，
cmakeの際に，エラーがでる．古いcmake関連のファイルを消したい場合．次を実行した後にcmakeする．
```shell
$ sh clean
```

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

[BEM_setBoundaryConditions.hpp#L68](BEM_setBoundaryConditions.hpp#L68):

## 多重節点
    多重節点という名前は具体性に欠ける．
    普通φnは(節点)にのみ依存する変数だが，nの変化が急なため，不連続性が著しい節点においては，(節点に加え面)にも依存する変数を複数設定する．それらは離散化などで使い分けることになる．
    BIEの離散化における，多重節点扱いについて．
    BIEを数値的に解くために，十分な数の１次方程式を作成する．これは，節点と同じ位置にBIEの原点を取ることで実現できる．
    同じ位置であるにもかかわらず，(節点に加え面)にも依存する変数φnを設定した場合，
    同じ位置であるにもかかわらず，それらを一つ一つを原点として，１次方程式を作成する．
    これらは完全に同じ方程式である．変数の数を節点の数よりも増やしたことによって，方程式の数が増えている．

[BEM_solveBVP.hpp#L650](BEM_solveBVP.hpp#L650):

このループでは，
               ある面integ_fに隣接する節点{p0,p1,p2}の列,IGIGn[origin(fixed),p0],...に値が追加されていく．
               （p0が多重接点の場合，適切にp0と同じ位置に別の変数が設定されており，別の面の積分の際にq0が参照される．）
               p0は，{面,補間添字}で決定することもできる．
               {面,補間添字0}->p0,{面,補間添字1}->p1,{面,補間添字2}->p2というように．
               //@ 多重節点：
               {面A,補間添字},{面B,補間添字},{面C,補間添字}が全て同じ節点p0を指していたとする．
               普通の節点なら，IGIGn[origin,{p0,nullptr}]を指す．
               多重節点なら，IGIGn[origin,{p0,面A}],IGIGn[origin,{p0,面B}]を指すようにする．
               この操作を言葉で言い換えると，
               「nが不連続に変化する点では，その点の隣接面にそれぞれ対してφnを求めるべきである（φは同じでも）．」
               「nが不連続に変化する点では，どの面を積分するかに応じて，参照するφnを区別し切り替える必要がある．」
               //
               //@ さて，この段階でp0が多重節点であるかどうか判断できるだろうか？
               {節点，面}-> 列ベクトルのインデックス を決めれるか？
               //
               面を区別するかどうかが先にわからないので，face*のままかnullptrとすべきかわからないということ．．．．
               //
               PBF_index[{p, Dirichlet, ある要素}]
               は存在しないだろう．Dirichlet節点は，{p, ある要素}からの寄与を，ある面に

[BEM_solveBVP.hpp#L681](BEM_solveBVP.hpp#L681):

# Example Function

 This is an example function that demonstrates how to use the keywords.

 **💡 NOTE:** This is a **💡 NOTE:**.
 **⚠️ WARNING:** This is a **⚠️ WARNING:**.
 **📝 TODO:** This is a **📝 TODO:** item.
 **❗ IMPORTANT:** This is an **❗ IMPORTANT:** point.
 **🌟 TIP:** This is a helpful **🌟 TIP:**.

