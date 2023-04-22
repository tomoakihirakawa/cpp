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

**NOTE:** `settingBEM.py`は`settingBEM.py`と同じフォルダ内にある必要がある．

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

