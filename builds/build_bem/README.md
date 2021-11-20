# コンパイルのやり方

以下のコマンドの先頭の"$"は無視してください．
本来はコンパイルの際には，多くのヘッダーファイルをインクルードするよう長いコンパイルのコマンドを打つ必要がある．
cmakeを使えば，それをCMakeLists.txtにあらかじめ書いておくことで省くことができる．

$ cmake -DCMAKE_BUILD_TYPE=Release ../

次に，

$ make

これでコンパイル終了．後は，次のようにすればmainファイルが実行される．

$ ./main

ただし，古いcmake情報が今のフォルダ内に残っている場合，その情報を削除しておかないと，
cmakeの際に，エラーがでる．古いcmake関連のファイルを消したい場合．次を実行した後にcmakeする．

$ sh clean

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


![](anim.gif)
