
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

# Fusion360でOBJファイルの作り方

単位はメートルに直して設計する．
また，Fusion360のxyzは，ここでのxzyなので注意する．

# 接触の確認

接触の有無は，まずバケツを利用して大まかに行う．

1. 互いの接触を確認したい物体（物体A）全てに対して，バケツを作成させる．メンバ関数としてバケツ作成関数が用意されている．
   * Network::makeBucketFaces
     * 点は面との接触の判定の際に，面の持つ線・点との接触も行うことができる．
   * Network::makeBucketPoints
   * Network::makeBucketParametricPoints
2. バケツクラスの持つ機能'Bucket<T>::getObject(座標)'等を利用して，指定可能な座標を中心とした正６面体状の範囲内に存在する（物体Aの）点や面のポインタを抜き出すことができる．
3. 抜き出したポインタは，物体Bの点や面のクラスに格納しておき，利用していく．
   より詳細な接触判定は，このときに行う．(例えば，networkPoint::addContactFaces)
   * addContactPoints
   * addContactFaces

物体は，自身の点や面が接触した物体の情報を取り出すことができる．また，接触対象をしていして抜き出すこともできる．
   * Network::getContactPointsOfPoints
   * Network::getContactFacesOfPoints
Pointは，即座に，自身が保持しているContactFacesに対する反射点を計算できる．


addContactFacesの修正を行なった．
ContactFacesはunordered_setからunordered_map変更した．


ポリゴン壁を利用したgradPの計算を行なったがうまくいっていない．
改めて，流体粒子が接触した壁上の座標を目で確認したい．
つまり，流体粒子の保有するContactFacesを再確認すること．



面との接触を判定．
接触面同士の接触を判定し，反射点をさらに増やす．
