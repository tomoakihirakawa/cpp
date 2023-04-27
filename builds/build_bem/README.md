## ../../include/Network.hpp

[../../include/Network.hpp#L1152](../../include/Network.hpp#L1152):

2021/09/02unordered_setを使うよう修正した
   将来的にはunordered setを返す関数に修正すべき

[../../include/Network.hpp#L1905](../../include/Network.hpp#L1905):

! 依存関係の明示方法
 下のように使う．これによって，setPointsFromLinesした後に，これを実行する必要があることを印象付けることができる．
 f->setGeometricProperties(f->setPointsFromLines())

[../../include/Network.hpp#L1913](../../include/Network.hpp#L1913):

@ networkFacesの持つ
         @ this->Points = {p0,p1,p2}
         @ this->Lines = {l0,l1,l2}
         @ の関係:
         @ 		    p2
         @         /\
         @ 		   /a2\
         @       /    \
         @  l2  /      \ l1
         @ 	   /a0    a1\
         @    -------------
         @  p0      l0      p1

[../../include/Network.hpp#L4301](../../include/Network.hpp#L4301):

p->setFaces()もf->setFaces()も,自身のp->Lines,f->Linesを元に，p->Faces,f->Facesを決定し保存する．
   p->Linesとf->Linesが正しく設定してあるかチェックする．
   特に，flipやdivideの後には，p->Lines,f->Linesが正しく設定されていない可能性があるので，要注意．

[../../include/Network.hpp#L4325](../../include/Network.hpp#L4325):

pointのFacesとLinesの関係は整合性があるか？faceのFacesとLinesの関係は整合性があるか？をチェック．
   setGeometricProperties()を実行していれば，この整合性は保たれるはずではある．

## ./BEM_setBoundaryConditions.hpp

[./BEM_setBoundaryConditions.hpp#L68](./BEM_setBoundaryConditions.hpp#L68):

## 多重節点
    多重節点という名前は具体性に欠ける．
    普通φnは(節点)にのみ依存する変数だが，nの変化が急なため，不連続性が著しい節点においては，(節点に加え面)にも依存する変数を複数設定する．それらは離散化などで使い分けることになる．
    BIEの離散化における，多重節点扱いについて．
    BIEを数値的に解くために，十分な数の１次方程式を作成する．これは，節点と同じ位置にBIEの原点を取ることで実現できる．
    同じ位置であるにもかかわらず，(節点に加え面)にも依存する変数φnを設定した場合，
    同じ位置であるにもかかわらず，それらを一つ一つを原点として，１次方程式を作成する．
    これらは完全に同じ方程式である．変数の数を節点の数よりも増やしたことによって，方程式の数が増えている．

## ./BEM_solveBVP.hpp

[./BEM_solveBVP.hpp#L414](./BEM_solveBVP.hpp#L414):

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

[./BEM_solveBVP.hpp#L445](./BEM_solveBVP.hpp#L445):

# Example Function

 This is an example function that demonstrates how to use the keywords.

 **💡 NOTE:** This is a **💡 NOTE:**.
 **⚠️ WARNING:** This is a **⚠️ WARNING:**.
 **📝 TODO:** This is a **📝 TODO:** item.
 **❗ IMPORTANT:** This is an **❗ IMPORTANT:** point.
 **🌟 TIP:** This is a helpful **🌟 TIP:**.

