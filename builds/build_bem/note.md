$`{\bf I}`$は慣性モーメントテンソル（2階のテンソル）．
実際の実験では，浮体のある基本的な姿勢における主慣性モーメント$`{\bf I}_{\rm principal}={(I_x, I_y, I_z)}`$が与えられる．
主慣性モーメント$`{\bf I}_{\rm principal}`$から，任意の浮体姿勢の慣性モーメントテンソルを求めるには，次のように考えればいい．

流体力の計算は空間に固定されたグローバルな座標系で行い，その上で流体力モーメント$`{\bf T}_{\rm G}`$が求まったとする．
これを，浮体の現在の姿勢のローカル座標系に変換する．
      
現在の姿勢を表すクォータニオンから作られる回転行列を$`{\bf R}`$とする．
この回転行列はグローバル座標系のベクトルにかけると，ローカル座標系のベクトルに変換させる．
反対に，ローカル座標系のベクトルにかけると，$`{\bf R}^{-1}={\bf R}^T`$をかけることで，グローバル座標系のベクトルに変換させる．
      
ローカルな座標系におけるモーメントは，$`{\bf T}_{\rm L}={\bf R}^{-1}\cdot{\bf T}_{\rm G}`$で計算できる．
これを主慣性モーメントで割ればローカルな座標系における角加速度となる，

$`\frac{d{\bf \Omega}_{\rm L}}{dt} = \frac{{\bf T}_{\rm L}}{{\bf I}_{\rm principal}}`$．

また，これをグローバルな座標系に戻すには，

$`\frac{d{\bf \Omega}_{\rm G}}{dt} = {\bf R}\cdot\frac{d{\bf \Omega}_{\rm L}}{dt} = {\bf R}\cdot\frac{{\bf T}_{\rm L}}{{\bf I}_{\rm principal}}= {\bf R}\cdot\frac{{\bf R}^{-1}\cdot{\bf T}_{\rm G}}{{\bf I}_{\rm principal}}`$．

これを書きかえると，

$`{\bf R}\cdot({\bf I}_{\rm principal} ({\bf R}^{-1}\cdot\frac{d{\bf \Omega}_{\rm G}}{dt})) = {\bf T}_{\rm G}`$．

これは次のようにまとめることができる．

$`\left\{{{\bf R}\cdot ({\bf I}_{\rm principal} {\bf R}^{-1})}\right\} \cdot \frac{d{\bf \Omega}_{\rm G}}{dt} = {\bf T}_{\rm G}`$．

この結果から，慣性モーメントテンソルは，$`{\bf I}={{\bf R}\cdot ({\bf I}_{\rm principal} {\bf R}^{-1})}`$であることがわかる．