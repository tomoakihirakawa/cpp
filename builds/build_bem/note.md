
## BIEの離散化について

\begin{equation}
\iiint_\Omega \left(G({\bf x},{\bf a})\nabla^2 \phi({\bf x}) - \phi({\bf x})\nabla^2 G({\bf x},{\bf a})\right)dV
= \iint_\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\end{equation}


\begin{equation}
\alpha ({\bf a})\phi({\bf a})
= \iint_\Gamma {\left({
\frac{1}{\|{\bf x}-{\bf a}\|}
\nabla \phi ({\bf{x}}) + \phi ({\bf{x}})
\frac{{\bf x}-{\bf a}}{\|{\bf x}-{\bf a}\|^3}}
\right) \cdot {\bf{n}}({\bf{x}})dS}
\end{equation}

面は面上の節点を使って補間され，面積分はこの補間された面上に沿って行われる．
補間に使うパラメタを$`{\bf \xi}=(\xi_0, \xi_1)`$として，
よく使われる３節点を使う線形補間を使うことにする．
元の面に対応する，線形補間面は，パラメタ上では$`{\xi_0 + \xi_1 = 1}`$を満たす範囲なので，
積分範囲は例えば$`0\leq \xi_0 \leq 1, 0\leq \xi_1 \leq 1-\xi_0`$となる．
しかし，数値積分につかう変数と重みの組み合わせは，コンパイルタイムに決めておき計算を効率化したいので，
この点で，変化する積分範囲は数値積分との相性が悪い．

### 積分範囲の変更と形状関数の修正

数値積分としてよく使われる，範囲$`0\leq \xi_0 \leq 1, 0\leq \xi_1 \leq 1`$のガウス求積を使うためには，
$`0\leq \xi_0 \leq 1, 0\leq \xi_1 \leq 1`$の範囲で元の面を表現する形状関数を使うのが望ましい．

\begin{equation}
\alpha ({\bf a})\phi({\bf a})
=\sum_{\rm each \, faces}\sum_{(\xi_0,w_0)} \sum_{(\xi_1,w_1)} 
w_0 w_1
{\left({
\frac{1}{\|{{\bf x}(\boldsymbol \xi)}-{\bf a}\|}
\nabla \phi ({\bf{x}}(\boldsymbol \xi)) + \phi ({\bf{x}}(\boldsymbol \xi))
\frac{{{\bf x}(\boldsymbol \xi)}-{\bf a}}{\|{{\bf x}(\boldsymbol \xi)}-{\bf a}\|^3}}
\right) \cdot {\bf{n}}({\bf{x}}(\boldsymbol \xi))
(1-\xi_0)
}
\end{equation}

ここで$`{\bf{x}}(\boldsymbol \xi)`$は，補間パラメタ$`\boldsymbol \xi`$によって表される面上の点である．

<!-- ----------------------------------------------------------------------- -->

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