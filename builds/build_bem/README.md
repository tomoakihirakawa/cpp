# Contents
- [🐋 BEM-MEL](#🐋-BEM-MEL)
    - [⛵ BEM-MEL について](#⛵-BEM-MEL-について)
        - [🪼 三角関数を使った古典的な解析手法](#🪼-三角関数を使った古典的な解析手法)
        - [🪼 周波数領域のBEM解析](#🪼-周波数領域のBEM解析)
            - [🐚 BEM周波数領域の問題点](#🐚-BEM周波数領域の問題点)
        - [🪼 時間領域のBEM解析（ほとんどの場合BEM-MEL）](#🪼-時間領域のBEM解析（ほとんどの場合BEM-MEL）)
            - [🐚 BEM-MEL の問題点](#🐚-BEM-MEL-の問題点)
        - [🪼 BEM-MEL の改良](#🪼-BEM-MEL-の改良)
        - [🪼 浮体動揺解析](#🪼-浮体動揺解析)
            - [🐚 なぜ今BEM-MELを開発するのか](#🐚-なぜ今BEM-MELを開発するのか)
    - [⛵ 入力ファイルの読み込み](#⛵-入力ファイルの読み込み)
    - [⛵ 計算プログラムの概要](#⛵-計算プログラムの概要)
        - [🪼 計算の流れ](#🪼-計算の流れ)
    - [⛵ 境界のタイプを決定する](#⛵-境界のタイプを決定する)
        - [🪼 多重節点](#🪼-多重節点)
        - [🪼 `getContactFaces()`や`getNearestContactFace()`の利用](#🪼-`getContactFaces()`や`getNearestContactFace()`の利用)
            - [🐚 `contact_angle`と`isInContact()`](#🐚-`contact_angle`と`isInContact()`)
            - [🐚 `addContactFaces()`](#🐚-`addContactFaces()`)
            - [🐚 呼び出し方法](#🐚-呼び出し方法)
        - [🪼 `uNeumann()`と`accelNeumann()`](#🪼-`uNeumann()`と`accelNeumann()`)
    - [⛵ 境界値問題](#⛵-境界値問題)
        - [🪼 基礎方程式](#🪼-基礎方程式)
        - [🪼 境界積分方程式（BIE）](#🪼-境界積分方程式（BIE）)
        - [🪼 BIEの離散化](#🪼-BIEの離散化)
        - [🪼 リジッドモードテクニック](#🪼-リジッドモードテクニック)
    - [⛵ 初期値問題](#⛵-初期値問題)
        - [🪼 流速$`\frac{d\bf x}{dt}`$の計算](#🪼-流速$`\frac{d\bf-x}{dt}`$の計算)
        - [🪼 $`\frac{d\phi}{dt}`$の計算](#🪼-$`\frac{d\phi}{dt}`$の計算)
        - [🪼 修正流速（激しい波の計算では格子が歪になりやすく，これがないと計算が難しい）](#🪼-修正流速（激しい波の計算では格子が歪になりやすく，これがないと計算が難しい）)
    - [⛵ 浮体動揺解析](#⛵-浮体動揺解析)
        - [🪼 浮体の運動方程式](#🪼-浮体の運動方程式)
        - [🪼 $`\phi _t`$と$`\phi _{nt}`$に関するBIEの解き方（と$`\phi _{nt}`$の与え方）](#🪼-$`\phi-_t`$と$`\phi-_{nt}`$に関するBIEの解き方（と$`\phi-_{nt}`$の与え方）)
            - [🐚 ディリクレ節点の$`\phi _{nt}`$の与え方(水面：圧力が既知，$`\phi`$が既知)](#🐚-ディリクレ節点の$`\phi-_{nt}`$の与え方(水面：圧力が既知，$`\phi`$が既知))
            - [🐚 ディリクレ節点の$`\phi _{t}`$の与え方($`\phi`$を与える造波装置：圧力が未知，$`\phi`$が既知)](#🐚-ディリクレ節点の$`\phi-_{t}`$の与え方($`\phi`$を与える造波装置：圧力が未知，$`\phi`$が既知))
            - [🐚 ノイマン節点での$`\phi _{nt}`$の与え方](#🐚-ノイマン節点での$`\phi-_{nt}`$の与え方)
        - [🪼 $`\phi _{nt}`$の計算で必要となる$`{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right)`$について．](#🪼-$`\phi-_{nt}`$の計算で必要となる$`{\bf-n}\cdot-\left({\nabla-\phi-\cdot-\nabla\nabla-\phi}\right)`$について．)
        - [🪼 浮体の重心位置・姿勢・速度の更新](#🪼-浮体の重心位置・姿勢・速度の更新)
        - [🪼 補助関数を使った方法](#🪼-補助関数を使った方法)
    - [⛵ 陽に与えられる境界条件に対して（造波装置など）](#⛵-陽に与えられる境界条件に対して（造波装置など）)
        - [🪼 フラップ型造波装置](#🪼-フラップ型造波装置)
        - [🪼 ピストン型造波装置](#🪼-ピストン型造波装置)
        - [🪼 正弦・余弦（`sin` もしくは `cos`）の運動](#🪼-正弦・余弦（`sin`-もしくは-`cos`）の運動)
    - [⛵ その他](#⛵-その他)
        - [🪼 境界値問題の未知変数](#🪼-境界値問題の未知変数)
        - [🪼 エネルギー保存則（計算精度のチェックに利用できる）](#🪼-エネルギー保存則（計算精度のチェックに利用できる）)
        - [🪼 内部流速の計算方法（使わなくてもいい）](#🪼-内部流速の計算方法（使わなくてもいい）)
        - [🪼 JSONファイルの出力](#🪼-JSONファイルの出力)
- [🐋 実行方法](#🐋-実行方法)
    - [⛵ ファイルのダウンロード](#⛵-ファイルのダウンロード)
    - [⛵ 入力ファイルの生成．](#⛵-入力ファイルの生成．)
    - [⛵ プログラムのコンパイルと実行](#⛵-プログラムのコンパイルと実行)
- [🐋 Input Generator](#🐋-Input-Generator)
- [🐋 Examples](#🐋-Examples)


---
# 🐋 BEM-MEL 

## ⛵ BEM-MEL について 

### 🪼 三角関数を使った古典的な解析手法 

水面がどのような微分方程式に従って運動するか調べると，
非粘性非圧縮渦なしを仮定しても，水面における境界条件は非線形である．

立てた連立偏微分方程式（境界条件と連続の式）を満たすような，関数，つまり解を，
三角関数の重ね合わせで求めようとすることは自然で賢い発想であり，この解析方法はある程度の成功を収めてきた．
この方法による線形理論はよく知られており水面波の基礎となっている．
また，摂動法を使って弱い非線形性をうまく取り込み三角関数で解を求めることもこの解析手法の延長線上にあり，よく行われている．

しかし，
複雑な形状を境界に持つ場合や，波が激しい場合においては，
この解析方法で課すことになる周期境界条件や，弱非線形性までしか考慮しないことが，
果たして結果に悪影響を及ぼさないか疑問である．
また，過渡的な現象，実際と同じように時間変化する現象に対する結果を得たい場合には，この解析手法では難しい．

### 🪼 周波数領域のBEM解析 

BEMを使った周波数領域の解析は，海洋工学の分野で標準的に行われているようだ．例えば，WAMITが有名．

周波数領域解析は以下のような流れである（[Kashiwagi et al. (2005)](https://linkinghub.elsevier.com/retrieve/pii/S0029801804002252)，[柏木 (2004)](https://www.jstage.jst.go.jp/article/technom/880/0/880_KJ00001033371/_article/-char/ja/)を参考）

速度ポテンシャルは

* 入射波ポテンシャル$`\phi _0`$
* 浮体による波の散乱ポテンシャル$`\phi _7`$
* 浮体動揺に放射ポテンシャル($`\phi _1,..\phi _6`$)．（浮体動揺の振幅が含まれており，運動方程式の解として後で求める）

の和として表す．各ポテンシャルに関して境界積分方程式を解くことで各ポテンシャルを求める．
次に，浮体にかかる流体力を浮体濡れ面上の圧力を積分し計算する．
結果として，**波浪強制力**，**付加質量**，**減衰力係数**が現れる．

$`\phi _0,\phi _7`$を含んだ力が得られる．これを**波浪強制力**とよぶ．
$`\phi _0+\phi _7`$をDiffractionポテンシャルと呼ぶこともあり，
これに関する境界値問題をDiffraction問題という（浮体の形状には関係するが浮体動揺に関係ない境界値問題でもある）

また，$`\phi _1,\phi _6`$を含んだ力は，
浮体加速度に比例する**付加質量**，浮体速度に比例する**減衰力係数**として整理できる．
これらに関する境界値問題をRadiation問題という（浮体動揺に関係する境界値問題でもある）

得られた，波浪強制力，付加質量，減衰力係数を使って，
浮体の運動方程式を解くことで，浮体動揺の振幅や浮体が受ける力を求めることができる．

#### 🐚 BEM周波数領域の問題点 

[Feng and Bai (2017)](https://linkinghub.elsevier.com/retrieve/pii/S0889974616300482)も指摘するように，
以上のBEM周波数領域解析は，線型理論であるため波が激しい場合の精度には疑問が残る．

例えば，
初めのポテンシャルの表現において，
Radiationポテンシャルは浮体動揺から切り離されている．
実際は，速度ポテンシャル自体が浮体動揺に応じて変化し，
さらにポテンシャルから計算される力自体も浮体姿勢に応じて変化する．

### 🪼 時間領域のBEM解析（ほとんどの場合BEM-MEL） 

1970 年代のコンピュータのメモリ容量は小さく，計算速度も遅かった．
当時開発された正方格子上でのシミュレーション手法を使って，
巻波砕破のシミュレーションを行おうと格子を細かくすると，
直ぐにメモリ容量を超えてしまい，また計算速度の問題もあって，正方格子を使った計算は現実的ではなかった．
これに対して，
[Longuet-Higgins and Cokelet (1976)](http://rspa.royalsocietypublishing.org/cgi/doi/10.1098/rspa.1976.0092)は，境界線上だけに計算点を設け，
その計算点の位置と速度ポテンシャルをラグランジュ的に時間発展させる方法を提案した．
水面で$`\frac{D\phi}{Dt}`$が簡単に計算できること，
流速(速度ポテンシャルの勾配)を計算するために，
$`\phi`$の接線方向微分$`\frac{\partial \phi}{\partial s}`$は節点の微分を使って，
法線方向微分$`\frac{\partial \phi}{\partial n}`$は境界積分方程式を解くことで計算できることを利用した．

<details style="background-color: rgba(144, 238, 144, 0.2);">
<summary>
💡 メモリ容量の変化
</summary>

<img src="./REVIEW/computer_memory.png" width="50%" />

現在，メモリ容量は，以前と比べて格段に大きくなり，メモリの節約を考える必要がなくなった．
また，領域型の計算手法で作られる代数連立１次方程式の係数行列は，疎行列であることが多く，
節点数は多けれども，疎行列なら反復解法を使って高速に解を求めることができる．
一方で，境界型の計算手法は，密行列を作る必要があり，密行列の生成には計算がかかる．

$`O(n _p^2)`$，$`O(n _d^3)`$
仮に$`n _p=L^3`$，$`n _d=6L^2`$としよう
$`O(L^6)`$，$`O(216 L^6)`$

つまり，当初の BEM-MEL の優位は，現在では他の手法にうばわれてしまっている．

</details>

#### 🐚 BEM-MEL の問題点 

BEM-MEL の結果に数値的な不安定が生じることは，[Longuet-Higgins and Cokelet (1976)](http://rspa.royalsocietypublishing.org/cgi/doi/10.1098/rspa.1976.0092)が既に紹介している．
計算精度を悪化させる原因は様々なものが考えられる．
例えば，係数行列を作成する際，つまり微分方程式を離散化する際に用いる，補間の精度や積分の精度．
または，時間発展の際に用いる，時間積分の精度などである．

補間と積分はセットで使うので，どちらが原因かを切り分けるのは難しい．
積分精度だけを考えても，数値積分手法の改良や，解析的な改良などが考えられる．
補間精度だけを考えても，補間手法の改良や，補間点の位置の調整などが考えられる．

### 🪼 BEM-MEL の改良 

続く研究目的は，BEM-MEL の改良に向けられた．

### 🪼 浮体動揺解析 

浮体の動揺解析を行うためには，次のようなステップを踏む．

1. 浮体に掛かる力とトルクを計算し，
2. 力と重心に関する運動方程式（トルクと角運動量に関する運動方程式）から加速度（角加速度）を求め，
3. 加速度（角加速度）を積分し速度（角速度）を更新し，
4. 速度（角速度）を積分し位置（姿勢）を更新する

浮体に掛かる圧力を面積分することで力を計算できるが，BEM-MEL では，圧力の計算で必要となる$`\phi _t`$が簡単には計算できない．
これは，FEM-MEL でも同じで，MEL を使った場合に共通雨したことである(これに関しては[Ma and Yan (2009)](http://doi.wiley.com/10.1002/nme.2505)に詳しく書かれている)．

Wu and {Eatock Taylor} (1996)や[Kashiwagi (2000)](http://journals.sagepub.com/doi/10.1243/0954406001523821)，[Wu and Taylor (2003)](www.elsevier.com/locate/oceaneng)の方法は，初めに$`\phi _t`$を計算し，次に圧力，力と計算して行くのではなく，
BIE と補助関数を使って，始めから圧力の面積分つまり力を別の変数の面積分として表した．
これと運動方程式を連立することで，直接，加速度を求めることができる．
[Feng and Bai (2017)](https://linkinghub.elsevier.com/retrieve/pii/S0889974616300482)は，この方法を発展させ２浮体の動揺解析を行っている．

本当に，複数の浮体に適用しにくい方法なのか？

#### 🐚 なぜ今BEM-MELを開発するのか 

高速多重極展開（FMM）を使った流体-物体相互作用解析は未だに達成されていないようで，
流体-物体相互作用解析においてBEM-MELが限界まで研究開発されたかというと，そうではないようだ．

##### なぜFMMを使った流体-物体相互作用解析が難しいのか？

流体-物体相互作用解析を行うには，２つの境界値問題を解く必要がある．
現在の計算手法は，一つ目の境界値問題で得られた係数行列の逆行列を再利用することで，計算コストを抑えている．

<details style="background-color: rgba(144, 238, 144, 0.2);">
<summary>
💡 逆行列を使う方法
</summary>

* **補助関数を使う方法**

補助関数（１浮体につき６つ増える）に関する境界値問題を解く必要がある（$`\phi`$-$`\phi _n`$に関するBIEの係数行列の行列を再利用することで高速化）．

[Feng and Bai (2017)](https://linkinghub.elsevier.com/retrieve/pii/S0889974616300482)から，補助関数を使う方法も係数行列が共通であるため計算コストを抑えることができるということがわかる．
言い換えれば，逆行列を使い回すことで，計算コストを抑えるということである．

> To compute the auxiliary functions, extra boundary value problems (BVPs) have to be solved. As th auxiliary functions share the same coefficient matrix with the velocity potential when proper boundary conditions are imposed, they are solved simultaneously with the potential, and not much additional computational effort is needed for solving these extra BVPs for the auxiliary functions.

しかし，FMMで利用されるGMRESのような反復解放を使うと，逆行列は求まらないし，計算が係数行列だけでなく右辺ベクトルにも依存するため，境界値問題をひとつひとつ個別に解く必要がある．

* **補助関数を使わない方法**

反復毎に境界値問題を解き$`\phi _t`$を計算する必要がある（$`\phi`$-$`\phi _n`$に関するBIEの係数行列の行列を再利用することで高速化）．

補助関数を使わない，$`\phi _t`$を反復して計算し徐々に収束させる方法も，逆行列がなければ反復毎に係数行列を求める必要がある．

</details>

そのため，逆行列に依存しない高速化手法が必要である．
その方法の指針として，境界値問題を一つにまとめることが考えられる．

[./main.cpp#L1](./main.cpp#L1)

---
## ⛵ 入力ファイルの読み込み 

1. 境界条件の設定
2. 境界値問題（BIE）を解き，$`\phi`$と$`\phi _n`$を求める
3. 三角形の線形補間を使って節点の流速を計算する

[./main.cpp#L196](./main.cpp#L196)

## ⛵ 計算プログラムの概要 

| 項目 | 詳細|
|---:|:---|
| 要素 | 線形三角要素 |
| 時間発展方法 | 4次のルンゲクッタ |
| 解析領域 | 時間領域 |
| 境界条件 | 水面の境界条件は非線形であるが，非線形のまま解く |

### 🪼 計算の流れ 

1. 境界条件の設定
2. 境界値問題（BIE）を解き，$`\phi`$と$`\phi _n`$を求める
3. 三角形の線形補間を使って節点の流速を計算する
4. 次時刻の$`\Omega(t+\Delta t)`$がわかるので，修正流速を計算する
5. 浮体の加速度を計算する．境界値問題（BIE）を解き，$`\phi _t`$と$`\phi _{nt}`$を求め，浮体面上の圧力$`p`$を計算する必要がある
6. 全境界面の節点の位置を更新．ディリクレ境界では$`\phi`$を次時刻の値へ更新

[./main.cpp#L369](./main.cpp#L369)

---
## ⛵ 境界のタイプを決定する 

0. 流体と物体の衝突を判定し，流体節点が接触する物体面を保存しておく．

* [`networkPoint::contact_angle`](../../include/networkPoint.hpp#L176)
* [`networkPoint::isInContact`](../../include/networkPoint.hpp#L192)
* [`networkPoint::addContactFaces`](../../include/networkPoint.hpp#L232)

を使って接触判定を行っている．

[流体が構造物との接触を感知する半径](../../builds/build_bem/BEM_setBoundaryTypes.hpp#L183)の設置も重要．

つぎに，その情報を使って，境界のタイプを次の順で決める．（物理量を与えるわけではない）

1. 面の境界条件：３節点全てが接触している流体面はNeumann面，それ以外はDirichlet面とする．CORNER面は設定しない．
- Neumann面$`\Gamma^{({\rm N})}`$ : 3点接触流体面
- Dirichlet面$`\Gamma^{({\rm D})}`$ : それ以外の面

2. 辺の境界条件 : 辺を含む２面がNeumann面ならNeumann辺，２面がDirichlet面ならDirichlet辺，それ以外はCORNERとする．
- Neumann辺 : 隣接面2面がNeumann面の辺
- Dirichlet辺 : 隣接面2面がDirichlet面の辺
- CORNER辺 : それ以外の辺（Neumann面とDirichlet面の間にある辺）

3. 点の境界条件：点を含む面全てがNeumann面ならNeumann点，面全てがDirichlet面ならDirichlet点，それ以外はCORNERとする．
- Neumann点 : 隣接面全てがNeumann面である点
- Dirichlet点 : 隣接面全てがDirichlet面である点
- CORNER点 : それ以外の点（Neumann面とDirichlet面の間にある点）

### 🪼 多重節点 

💡 面の向き$`\bf n`$がカクッと不連続に変わる節点には，$`\phi`$は同じでも，隣接面にそれぞれ対して異なる$`\phi _n`$を計算できるようにする

💡 $`\bf n`$が不連続に変化する節点まわりの要素は，自分のために用意された$`\phi _n`$を選択し補間に用いなければならない

これを多重節点という．

[./BEM_setBoundaryTypes.hpp#L7](./BEM_setBoundaryTypes.hpp#L7)

### 🪼 `getContactFaces()`や`getNearestContactFace()`の利用 

#### 🐚 `contact_angle`と`isInContact()` 

| `networkPoint`のメンバー関数/変数      | 説明                                                                |
|-------------------------|--------------------------------------------------------------------------------|
| [`contact_angle`](../../include/networkPoint.hpp#L176)         | ２面の法線ベクトルがこの`contact_angle`大きい場合，接触判定から除外される |
| [`isFacing()`](../../include/networkPoint.hpp#L179)       | ２面の法線ベクトルが`contact_angle`よりも小さいか判定する．ただし，角度は，向かい合う面がなす最小の角度と考える |
| [`isInContact()`](../../include/networkPoint.hpp#L192)         | 点の隣接面のいずれかが，与えられた面と接触しているか判定する．範囲内で接触しており，かつ`isFacing`が真である場合`true`を返す． |
| [`addContactFaces()`](../../include/networkPoint.hpp#L232)     | バケツに保存された面を基に，節点が接触した面を`networkPoint::ContactFaces`に登録する．   |
[../../include/networkPoint.hpp#L165](../../include/networkPoint.hpp#L165)


#### 🐚 `addContactFaces()` 

| `networkPoint`のメンバー関数/変数      | 説明                                                                |
|-------------------------|--------------------------------------------------------------------------------|
| `addContactFaces()`     | バケツに保存された面を基に，節点が接触した面を`networkPoint::ContactFaces`に登録する．   |
| `std::unordered_set<networkFace *> ContactFaces`          | 節点が接触した面が登録されている．   |
| `std::tuple<networkFace *, Tddd> nearestContactFace`    | 節点にとって最も近い面とその座標を登録されている．       |
| `std::unordered_map<networkFace *, std::tuple<networkFace *, Tddd>> f_nearestContactFaces` | この節点に隣接する各面にとって，最も近い面とその座標をこの変数に登録する．           |
[../../include/networkPoint.hpp#L285](../../include/networkPoint.hpp#L285)


#### 🐚 呼び出し方法 

* `getContactFaces()`で`ContactFaces`呼び出せる．
* `getNearestContactFace()`で`nearestContactFace`呼び出せる．
* `getNearestContactFace(face)`で`f_nearestContactFaces`呼び出せる．
[../../include/Network.hpp#L884](../../include/Network.hpp#L884)


これらは，`uNeumann()`や`accelNeumann()`で利用される．

### 🪼 `uNeumann()`と`accelNeumann()` 

接触している物体が，剛体でない場合，
`velocity_of_Body`は，物体の節点（ `networkPoint` ）の速度（加速度）を元にして速度（加速度）を計算する．
そのため，`networkPoint::velocity`や`networkPoint::accel`を設定しておく必要がある．

`uNeumann(p, const adjacent_f)`や`accelNeumann(p, const adjacent_f)`
を使う時は，必ず`adjacent_f`が`p`に**隣接面するノイマン面**であることを確認する．

[./BEM_utilities.hpp#L317](./BEM_utilities.hpp#L317)

---
## ⛵ 境界値問題 

### 🪼 基礎方程式 

```math
\begin{align}
\nabla\cdot\nabla \phi& = 0&&\text{in}&&{\bf x} \in \Omega(t),\\
\frac{\partial\phi}{\partial t} +\frac{1}{2}\nabla\phi\cdot\nabla\phi - g z &=0 &&\text{on}&&{\bf x} \in \Gamma^{(\rm D)}(t),\\
\phi _n + {{\bf u} _b}\cdot{{\bf n} _b} &=0&&\text{on}&&{\bf x}\in \Gamma^{(\rm N)}(t),
\end{align}
```

ここで，
$`{\bf x} ={(x,y,z)}`$は空間座標，$`{\bf u} _b`$は物体の流速，
$`{\bf n} _b`$は物体の外向き単位法線ベクトル，
$`\nabla=(\frac{\partial}{\partial x},\frac{\partial}{\partial y},\frac{\partial}{\partial z})`$
である．
また，$`\phi _n`$は境界面上での外向き法線方向の流速を表し，
境界面上の外向き単位法線ベクトル$`\bf n`$を使えば$`\phi _n ={\nabla\phi}\cdot {\bf n}`$で表される．

### 🪼 境界積分方程式（BIE） 

**グリーンの定理**

任意の$`\phi`$，$`G`$に対して次が成り立つ（**グリーンの定理**）．

```math
\iiint _\Omega \left(G({\bf x},{\bf a})\nabla^2 \phi({\bf x}) - \phi({\bf x})\nabla^2 G({\bf x},{\bf a})\right)dV
= \iint _\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
```


$`\phi`$がラプラス方程式$`\nabla^2\phi=0`$を満たし，$`G=1/\|{\bf x}-{\bf a}\|`$とすると，
グリーンの定理から$`\phi`$と$`\phi _n`$の関係式，BIEが得られる．

```math
\alpha ({\bf{a}})\phi ({\bf{a}}) = \iint _\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi ({\bf{x}}) - \phi ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
```

ここで，$`{\bf a}`$は境界面上の位置ベクトルであり，この原点$`{\bf a}`$を固定し$`{\bf x}`$について面積分される．
$`G`$は任意のスカラー関数で$`G=1/\|{\bf x}-{\bf a}\|`$とすることで，グリーンの定理の体積積分が消え，BIEの左辺のように，
原点での立体角$`\alpha\left( {\bf{a}} \right)`$とポテンシャル$`\phi( {\bf{a}})`$の積だけが残る．

<img src="schematic_BIE.png" width="400px">

この式は，流体内部では，$`\alpha ({\bf{a}})`$は$`1`$とできる．
この式は，$`\bf{a}`$におけるポテンシャル$`\phi ({\bf{a}})`$が，右辺の１重層ポテンシャルと２重層ポテンシャルの和で表されることを示している．
$`G=1/\|{\bf x}-{\bf a}\|`$がラプラス法廷式の基本解であり，$`\phi`$は境界におけるポテンシャルの分布である．

[./BEM_solveBVP.hpp#L7](./BEM_solveBVP.hpp#L7)

### 🪼 BIEの離散化 

BIEを線形三角要素とGauss-Legendre積分で離散化すると，

```math
\sum\limits _{k _\vartriangle}\sum\limits _{{\xi _1},{w _1}} {\sum\limits _{{\xi _0},{w _0}} {\left( {{w _0}{w _1}\left( {\sum\limits _{j=0}^2 {{{\left( {{\phi _n}} \right)} _{k _\vartriangle,j }}{N _{j }}\left( \pmb{\xi } \right)} } \right)\frac{1}{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x} _{i _\circ}}} \|}}\left\|\frac{{\partial{\bf{x}}}}{{\partial{\xi _0}}} \times \frac{{\partial{\bf{x}}}}{{\partial{\xi _1}}}\right\|} \right)} }=
```
```math
\alpha _{i _\circ}(\phi) _{i _\circ}-\sum\limits _{k _\vartriangle}\sum\limits _{{\xi _1},{w _1}} \sum\limits _{{\xi _0},{w _0}} {\left( {{w _0}{w _1}\left({\sum\limits _{j =0}^2{{{\left( \phi  \right)} _{k _\vartriangle,j }}{N _{j}}\left( \pmb{\xi } \right)} } \right)\frac{\bf{x}(\pmb{\xi})-{{\bf x} _{i _\circ} }}{{{{\| {{\bf{x}}\left( \pmb{\xi } \right) - {{\bf x} _{i _\circ}}}\|}^3}}} \cdot\left(\frac{{\partial {\bf{x}}}}{{\partial {\xi _0}}}\times\frac{{\partial {\bf{x}}}}{{\partial {\xi _1}}}\right)}\right)}
```

ここで，$`\phi _{k _\vartriangle,j}`$における$`k _\vartriangle`$は三角形要素の番号，$`j`$は三角形要素の頂点番号．
$`N _j`$は三角形要素の形状関数，$`\pmb{\xi}`$は三角形要素の内部座標，$`w _0,w _1`$はGauss-Legendre積分の重み，$`\alpha _{i _\circ}`$は原点$`i _\circ`$における立体角，$`\phi`$はポテンシャル，$`\phi _n`$は法線方向のポテンシャル，$`\bf{x}`$は空間座標，$`{\bf x} _{i _\circ}`$は原点の空間座標である．

形状関数$`{\pmb N} _j({\pmb \xi}),{\pmb \xi}=(\xi _0,\xi _1)`$は，$`\xi _0,\xi _1`$が$`0`$から$`1`$動くことで，範囲で三角要素全体を動くように定義している．

```math
{\pmb N}({\pmb \xi}) = (N _0({\pmb \xi}),N _1({\pmb \xi}),N _2({\pmb \xi})) = (\xi _0, - \xi _1 (\xi _0 - 1), (\xi _0-1)(\xi _1-1))
```

[./BEM_solveBVP.hpp#L195](./BEM_solveBVP.hpp#L195)

このループでは，BIEの連立一次方程式の係数行列`IGIGn`を作成する作業を行なっている．
`IGIGn`は，ある節点$`i _\circ`$（係数行列の行インデックス）に対する
他の節点$`j _\circ`$（係数行列の列インデックス）の影響度合いのようなものである．
その影響度合いは，他の節点$`j _\circ`$の所属する要素までの距離や向きによって決まることが離散化された式からわかる．

| Variable | Description |
|:--------:|:-----------:|
| `origin` | 原点となる節点$`i _\circ`$ |
| `integ_f` | Element $`k _{\triangle}`$ |
| `t0, t1, ww` | Gaussian points and thier wieghts $`\xi _0, \xi _1, w _0 w _1`$ |
| `p0, p1, p2` | Node of the element $`k _{\triangle}`$ |
| `N012` | Shape function $`\pmb{N} _j`$ |
| `IGIGn` | Coefficient matrices of the left and right sides |
| `nr` | $`\| \pmb{x} - \pmb{x} _{i\circ } \|`$ |
| `tmp` | $`w _0 w _1 \frac{1 - \xi _0}{\| \pmb{x} - \pmb{x} _{i\circ } \|}`$ |
| `cross` | $`\frac{\partial \pmb{x}}{\partial \xi _0} \times \frac{\partial \pmb{x}}{\partial \xi _1}`$ |

[./BEM_solveBVP.hpp#L259](./BEM_solveBVP.hpp#L259)

### 🪼 リジッドモードテクニック 

全て$`\phi=1`$とすると，$`\alpha({\bf a}) = -\int\int{\nabla G({\bf x},{\bf a})\cdot{\bf n}({\bf x})dS}`$となり，これを離散化すると，数値積分による評価が難しかった係数行列の対角成分がより精確に計算できる．
これはリジッドモードテクニックと呼ばれている．
$`{\bf x} _{i\circ}`$が$`{\bf x}({\pmb \xi})`$に近い場合，$`G`$は急激に特異的に変化するため，数値積分精度が悪化するが，リジッドモードテクニックによって積分を回避できる．

[./BEM_solveBVP.hpp#L334](./BEM_solveBVP.hpp#L334)

係数行列`IGIGn`は，左辺の$`I _G \phi _n`$，右辺の$`I _{G _n}\phi`$の係数．

```math
(I _G) _{i _\circ,j _\circ} (\phi _n) _{j _\circ} = (I _{Gn}) _{i _\circ,j _\circ}  \phi _{j _\circ}
```

境界条件に応じて，未知変数は$`\phi,\phi _n`$のどちらかに決まる．
未知変数が$`\phi`$の場合（Dirichlet境界条件の場合），
係数行列`IGIGn`中で対応する列を符号変えて入れ替えることで移項したことになる．


移項前:
```math
\begin{bmatrix}I _{G0} & I _{G1} & I _{G2} & I _{G3}\end{bmatrix} \begin{bmatrix}\phi _{n0} \\ \phi _{n1} \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}I _{Gn0} & I _{Gn1} & I _{Gn2} & I _{Gn3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _1 \\ \phi _2 \\ \phi _3\end{bmatrix}
```

移項後:
```math
\begin{bmatrix}I _{G0} & -I _{Gn1} & I _{G2} & I _{G3}\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}I _{Gn0} & -I _{G1} & I _{Gn2} & I _{Gn3}\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}
```

多重節点(1と3が多重節点の場合):
```math
\begin{bmatrix}0 & 1 & 0 & 0\end{bmatrix}\begin{bmatrix}\phi _{n0} \\ \phi _1 \\ \phi _{n2} \\ \phi _{n3}\end{bmatrix} =\begin{bmatrix}0 & 0 & 0 & 1\end{bmatrix}\begin{bmatrix}\phi _0 \\ \phi _{n1} \\ \phi _2 \\ \phi _3\end{bmatrix}
```

[./BEM_solveBVP.hpp#L372](./BEM_solveBVP.hpp#L372)

---
## ⛵ 初期値問題 

節点の位置と速度ポテンシャル$`\phi`$に関する初期値問題を解いて行くことが，シミュレーションである．
言い換えると，節点位置$`\frac{d\bf x}{dt}`$と速度ポテンシャル$`\frac{d\phi}{dt}`$を少しずつ$`\Delta t`$ずつ時間積分することが，シミュレーションである．
ちなみに，$`\frac{d\bf x}{dt}`$や$`\frac{d\phi}{dt}`$を計算するには，境界値問題を解く必要がある．

ある時刻において，境界値問題が解けたら，$`\frac{d\bf x}{dt}`$と$`\frac{d\phi}{dt}`$はどのように計算できるだろうか．

### 🪼 流速$`\frac{d\bf x}{dt}`$の計算 

ある三角要素上の接線流速$`\nabla \phi _{\parallel}`$は，線形三角要素補間を使って次のように計算する．

```math
\nabla \phi _{\parallel} = \frac{\bf n}{2A} \times (({\bf x} _2 - {\bf x} _1) \phi _0 +({\bf x} _0 - {\bf x} _2) \phi _1 + ({\bf x} _1 - {\bf x} _0) \phi _2)
```

三角要素上の流速$`\nabla \phi`$は，次のように計算する．

```math
\nabla \phi = \frac{(\phi _n) _0+(\phi _n) _1+(\phi _n) _2}{3} {\bf n} + \nabla \phi _{\parallel}
```

### 🪼 $`\frac{d\phi}{dt}`$の計算 

ある流体粒子に乗ってみたときの，速度ポテンシャルの時間変化$`\frac{D \phi}{D t}`$は，次のように計算できる．

```math
\frac{D \phi}{D t} = \frac{\partial \phi}{\partial t} + \nabla \phi \cdot \nabla \phi
```

<details style="background-color: rgba(144, 238, 144, 0.2);">
<summary>
💡 オイラー的記述
</summary>

$`\phi=\phi(t,{\bf x})`$のように書き表し，位置と空間を独立させ分けて考える方法を，オイラー的記述という．こう書くと，$`\frac{d \phi}{d t}`$は，$`\frac{\partial \phi}{\partial t}`$であり，これは，速度ポテンシャルの純粋な時間変化ではない．純粋な，ある流体粒子の速度ポテンシャルの時間変化を表すためには，位置が時間によって変わると考え，つまり$`\phi=\phi(t,{\bf x}(t))`$と一時的に考えなおし，そして，時間微分する．そうすると$`\frac{d\phi}{dt} = \frac{\partial \phi}{\partial t} + \frac{d\bf x}{dt}\cdot \nabla \phi`$となる．

</details>

ここの$`\frac{\partial \phi}{\partial t}`$の計算は簡単ではない．そこで，ベルヌーイの式（大気圧と接する水面におけるベルヌーイの式は圧力を含まず簡単）を使って，$`\frac{\partial \phi}{\partial t}`$を消去する．

[./BEM_utilities.hpp#L495](./BEM_utilities.hpp#L495)

---
### 🪼 修正流速（激しい波の計算では格子が歪になりやすく，これがないと計算が難しい） 

ディリクレ節点（水面）：

求めた流速から，次の時刻の境界面$`\Omega(t+\Delta t)`$を見積もり，その面上で節点を移動させ歪さを解消する．
修正ベクトルは，$`\Delta t`$で割り，求めた流速$`\nabla \phi`$に足し合わせて，節点を時間発展させる．

ノイマン節点：

ノイマン節点も修正流速を加え時間発展させる．
ただし，ノイマン節点の修正流速に対しては，節点が水槽の角から離れないように，工夫を施している．

[`calculateVecToSurface`](../../builds/build_bem/BEM_calculateVelocities.hpp#L228)で$`\Omega(t+\Delta t)`$上へのベクトルを計算する．

1. まず，[`vectorTangentialShift`](../../builds/build_bem/BEM_calculateVelocities.hpp#L134)で接線方向にシフトし，
2. [`vectorToNextSurface`](../../builds/build_bem/BEM_calculateVelocities.hpp#L143)で近くの$`\Omega(t+\Delta t)`$上へのベクトルを計算する．

[./BEM_calculateVelocities.hpp#L207](./BEM_calculateVelocities.hpp#L207)

---
## ⛵ 浮体動揺解析 

BEM-MELで浮体動揺解析ができるようにするのは簡単ではない．
浮体に掛かる圧力の計算に必要な$`\phi _t`$が簡単には求まらないためである．
これに関しては，[Wu and Taylor (2003)](www.elsevier.com/locate/oceaneng)が参考になる．

### 🪼 浮体の運動方程式 

<img src="schematic_float.png" width="400px" />

浮体の重心の運動方程式：

```math
m \frac{d {\boldsymbol U} _{\rm c}}{d t} = \boldsymbol{F} _{\text {ext }}+\boldsymbol{F} _{\text {hydro }}, \quad
\boldsymbol{I} \frac{d {\boldsymbol \Omega} _{\rm c}}{d t} = \boldsymbol{T} _{\text {ext }}+\boldsymbol{T} _{\text {hydro }}
```

$`{\boldsymbol U} _{\rm c}`$は浮体の移動速度．
$`\boldsymbol{F} _{\text {ext }}`$は重力などの外力，$`\boldsymbol{F} _{\text {hydro }}`$は水の力，$`\boldsymbol{T} _{\text {ext }}`$は外力によるトルク，$`\boldsymbol{T} _{\text {hydro }}`$は水の力によるトルク．
浮体が流体から受ける力$`\boldsymbol{F} _{\text {hydro }}`$は，浮体表面の圧力$`p`$を積分することで得られ，
また圧力$`p`$は速度ポテンシャル$`\phi`$を用いて，以下のように書ける．

[圧力積分](../../builds/build_bem/BEM_solveBVP.hpp#L117)と
[トルクの積分](../../builds/build_bem/BEM_solveBVP.hpp#L104)：

```math
\boldsymbol{F} _{\text {hydro }}=\iint _{\Gamma _{\rm float}} p\boldsymbol{n}  d S, \quad
\boldsymbol{T} _{\text {hydro }}=\iint _{\Gamma _{\rm float}} ({\bf x}-{\bf x} _{\rm c})\times (p\boldsymbol{n})  d S, \quad
p= p({\bf x}) =-\rho\left(\frac{\partial \phi}{\partial t}+\frac{1}{2} (\nabla \phi)^{2}+g z\right)
```

$`\frac{\partial \phi}{\partial t}`$を$`\phi _t`$と書くことにする．この$`\phi _t`$は陽には求められない．
そこで，$`\phi`$と似た方法，BIEを使った方法で$`\phi _t`$を求める．$`\phi`$と$`\phi _n`$の間に成り立つ境界積分方程式と全く同じ式が，$`\phi _t`$と$`\phi _{nt}`$の間にも成り立つ：

```math
\alpha ({\bf{a}})\phi _t ({\bf{a}}) = \iint _\Gamma {\left( {G({\bf{x}},{\bf{a}})\nabla \phi _t ({\bf{x}}) - \phi _t ({\bf{x}})\nabla G({\bf{x}},{\bf{a}})} \right) \cdot {\bf{n}}({\bf{x}})dS}
\quad\text{on}\quad{\bf x} \in \Gamma(t).
```

[./BEM_solveBVP.hpp#L567](./BEM_solveBVP.hpp#L567)

### 🪼 $`\phi _t`$と$`\phi _{nt}`$に関するBIEの解き方（と$`\phi _{nt}`$の与え方） 

$`\phi _t`$と$`\phi _{nt}`$に関するBIEを解くためには，ディリクレ境界には$`\phi _t`$を，ノイマン境界には$`\phi _{nt}`$を与える．

#### 🐚 ディリクレ節点の$`\phi _{nt}`$の与え方(水面：圧力が既知，$`\phi`$が既知) 

このディリクレ境界では，圧力が与えられていないので，このBiEにおいては，ノイマン境界条件を与える．
ただし，壁が完全に固定されている場合，$`\phi _{nt}`$は0とする．

#### 🐚 ディリクレ節点の$`\phi _{t}`$の与え方($`\phi`$を与える造波装置：圧力が未知，$`\phi`$が既知) 

ディリクレ境界では$`\phi _t`$は，圧力が大気圧と決まっているので，ベルヌーイの圧力方程式から$`\phi _t`$を求めることができる．

#### 🐚 ノイマン節点での$`\phi _{nt}`$の与え方 

境界面が静止しているかどうかに関わらず，流体と物体との境界では，境界法線方向速度が一致する．
境界面上の点の位置ベクトルを$`\boldsymbol r`$とする．
表面上のある点の移動速度$`\frac{d\boldsymbol r}{dt}`$と流体粒子の流速$`\nabla \phi`$の間には，次の境界条件が成り立つ．

```math
{\bf n}\cdot\frac{d\boldsymbol r}{dt} =  {\bf n} \cdot \nabla \phi,\quad \frac{d\boldsymbol r}{dt} = \boldsymbol U _{\rm c} + {\boldsymbol \Omega} _{\rm c} \times \boldsymbol r
```

物体上のある点ではこれが常に成り立つ．

これを微分することで，$`\phi _{nt}`$を$`\phi`$と加速度$`\frac{d{\boldsymbol U} _{\rm c}}{dt}`$と角加速度$`\frac{d{\boldsymbol \Omega} _{\rm c}}{dt}`$を使って表すことができる．
[Wu (1998)](https://www.sciencedirect.com/science/article/pii/S088997469890158X)

```math
\begin{aligned}
&\rightarrow& 0& =\frac{d}{dt}\left({\bf n}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)\right) \\
&\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \frac{d}{dt}\left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)\\
&\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2}-\left(\frac{\partial}{\partial t}+\frac{d{\boldsymbol r}}{dt}\cdot\nabla\right)\nabla \phi\right)\\
&\rightarrow& 0& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2}- {\nabla \phi _t - \left(\frac{d\boldsymbol r}{dt} \cdot \nabla\right)\nabla \phi}\right)\\
&\rightarrow& \phi _{nt}& =\frac{d{\bf n}}{dt}\cdot \left(\frac{d\boldsymbol r}{dt}-\nabla \phi\right)+ {\bf n}\cdot \left(\frac{d^2\boldsymbol r}{dt^2} - \frac{d\boldsymbol r}{dt} \cdot (\nabla\otimes\nabla \phi) \right)
\end{aligned}
```

ここの$`\frac{d{\bf n}}{dt}`$と$`\frac{d^2\boldsymbol r}{dt^2}`$は，$`{\boldsymbol U} _{\rm c}`$と$`\boldsymbol \Omega _{\rm c}`$を用いて，

```math
\frac{d^2\boldsymbol r}{dt^2} = \frac{d}{dt}\left({\boldsymbol U} _{\rm c} + \boldsymbol \Omega _{\rm c} \times \boldsymbol r\right),\quad \frac{d{\bf n}}{dt} = {\boldsymbol \Omega} _{\rm c}\times{\bf n}
```

[`phin_Neuamnn`](../../builds/build_bem/BEM_utilities.hpp#L675)で$`\phi _{nt}`$を計算する．これは[`setPhiPhin_t`](../../builds/build_bem/BEM_solveBVP.hpp#L780)で使っている．

$`\frac{d^2\boldsymbol r}{dt^2}`$を上の式に代入し，$`\phi _{nt}`$を求め，
次にBIEから$`\phi _t`$を求め，次に圧力$p$を求める．
そして，浮体の重さと慣性モーメントを考慮して圧力から求めた$`\frac{d^2\boldsymbol r}{dt^2}`$は，
入力した$`\frac{d^2\boldsymbol r}{dt^2}`$と一致しなければならない．

現状を整理すると，この浮体動揺解析において，知りたい未知変数は，浮体の加速度と角加速度だけ．
しかし，浮体の没水面上にある節点での圧力$`p`$が得られないと，$`\boldsymbol{F} _{\text {hydro }}`$が得られず，運動方程式から浮体加速度が計算できない．
圧力を計算するためには，$`\phi _t`$が必要で，$`\phi _t`$は簡単には得られない，という状況．

物体の加速度は， 節点における$`\{\phi _{nt0},\phi _{nt1},\phi _{nt2},..\} = \Phi _{nt}`$が分かれば求まるが，
逆に$`\phi _{nt}`$は$`\frac{d\boldsymbol U _{\rm c}}{dt}`$と$\frac{d {\boldsymbol \Omega} _{\rm c}}{d t}$が分かれば求まる．また，物体の角加速度に関しても同様である．

```math
m \frac{d\boldsymbol U _{\rm c}}{dt} = \boldsymbol{F} _{\text {ext }}+ F _{\text {hydro}}\left(\Phi _{nt}\left(\frac{d\boldsymbol U _{\rm c}}{dt},\frac{d {\boldsymbol \Omega} _{\rm c}}{d t}\right)\right),\quad
\boldsymbol{I} \frac{d {\boldsymbol \Omega} _{\rm c}}{d t} = \boldsymbol{T} _{\text {ext }}+\boldsymbol{T} _{\text {hydro }}\left(\Phi _{nt}\left(\frac{d\boldsymbol U _{\rm c}}{dt},\frac{d {\boldsymbol \Omega} _{\rm c}}{d t}\right)\right)
```

これを満たすように，$`\Phi _{nt}`$を求める．これは次のように書き換えて，根探し問題として解く．
このプログラムでは，[Broyden法](../../builds/build_root_finding/example1_Broyden.cpp#L22)を使って，根探している．

```math
\boldsymbol{0} = m \frac{d\boldsymbol U _{\rm c}}{dt} - \boldsymbol{F} _{\text {ext }} - F _{\text {hydro}}\left(\Phi _{nt}\left(\frac{d\boldsymbol U _{\rm c}}{dt},\frac{d {\boldsymbol \Omega} _{\rm c}}{d t}\right)\right),\quad
\boldsymbol{0} = \boldsymbol{I} \frac{d {\boldsymbol \Omega} _{\rm c}}{d t} - \boldsymbol{T} _{\text {ext }} - \boldsymbol{T} _{\text {hydro }}\left(\Phi _{nt}\left(\frac{d\boldsymbol U _{\rm c}}{dt},\frac{d {\boldsymbol \Omega} _{\rm c}}{d t} \right)\right)
```

この式を，$`{\boldsymbol Q}\left(\dfrac{d {\boldsymbol U} _{\rm c}}{d t}, \dfrac{d {\boldsymbol \Omega} _{\rm c}}{d t}\right)=(0,0,0,0,0,0)`$
として，これを満たすような$`\dfrac{d {\boldsymbol U} _{\rm c}}{d t}`$と$`\dfrac{d {\boldsymbol \Omega} _{\rm c}}{d t}`$を求める．
$`\phi _{nt}`$はこれを満たした$`\dfrac{d {\boldsymbol U} _{\rm c}}{d t}`$と$`\dfrac{d {\boldsymbol \Omega} _{\rm c}}{d t}`$を用いて求める．

$`\phi _{nt}`$は，[ここ](../../builds/build_bem/BEM_solveBVP.hpp#L794)で与えている．

[./BEM_solveBVP.hpp#L610](./BEM_solveBVP.hpp#L610)

```math
\nabla\otimes{\bf u} = \nabla \otimes \nabla \phi =
\begin{bmatrix} \phi _{xx} & \phi _{xy} & \phi _{xz} \\
\phi _{yx} & \phi _{yy} & \phi _{yz} \\
\phi _{zx} & \phi _{zy} & \phi _{zz}
\end{bmatrix}
```

ヘッセ行列の計算には，要素における変数の勾配の接線成分を計算する[`HessianOfPhi`](../../builds/build_bem/BEM_utilities.hpp#L647)を用いる．
節点における変数を$`v`$とすると，$`\nabla v-{\bf n}({\bf n}\cdot\nabla v)`$が計算できる．
要素の法線方向$`{\bf n}`$が$`x`$軸方向$`{(1,0,0)}`$である場合，$`\nabla v - (\frac{\partial}{\partial x},0,0)v`$なので，
$`(0,\frac{\partial v}{\partial y},\frac{\partial v}{\partial z})`$が得られる．

[./BEM_solveBVP.hpp#L691](./BEM_solveBVP.hpp#L691)

### 🪼 $`\phi _{nt}`$の計算で必要となる$`{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right)`$について． 

$`\nabla`$を，$`(x,y,z)`$の座標系ではなく，
面の法線方向$`{\bf n}`$を$`x`$の代わりにとり，
面に水平な方向を$`t _0,t _1`$とする座標系で考えることにして，$`\nabla^\ast`$と書くことにする．
$`{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right)`$では，$`{\bf n}`$方向成分だけをとる操作をしているので，
新しい座標系でも同じようにすれば，結果は変わらない．

```math
{\bf n}\cdot \left({\nabla \phi \cdot \nabla\nabla \phi}\right) =  {(1,0,0)}\cdot\left({\nabla^* \phi \cdot \nabla^* \nabla^* \phi}\right).
\quad
\nabla^* \phi = \left(\phi _n, \phi _{t _0}, \phi _{t _1}\right),
\quad \nabla^* \nabla^* \phi =
\begin{bmatrix}
\phi _{nn} & \phi _{nt _0} & \phi _{nt _1} \\
\phi _{t _0n} & \phi _{t _0t _0} & \phi _{t _0t _1} \\
\phi _{t _1n} & \phi _{t _1t _0} & \phi _{t _1t _1}
\end{bmatrix}
```

最後に第１成分だけが残るので，

```math
{(1,0,0)}\cdot\left({\nabla^* \phi \cdot \nabla^* \nabla^* \phi}\right) = \nabla^* \phi \cdot (\phi _{nn}, \phi _{t _0n}, \phi _{t _1n})
```

$`\phi _{nn}`$は，直接計算できないが，ラプラス方程式から$`\phi _{nn}=- \phi _{t _0t _0}- \phi _{t _1t _1}`$となるので，水平方向の勾配の計算から求められる．

[./BEM_utilities.hpp#L632](./BEM_utilities.hpp#L632)

### 🪼 浮体の重心位置・姿勢・速度の更新 

浮体の重心位置は，重心に関する運動方程式を解くことで求める．
姿勢は，角運動量に関する運動方程式などを使って，各加速度を求める．姿勢はクオータニオンを使って表現する．

[./main.cpp#L487](./main.cpp#L487)

---
### 🪼 補助関数を使った方法 

浮体動揺解析で問題となったのは，圧力の計算に使う$`\phi _t\,{\rm on}\,🚢`$が簡単には求まらないことであったが，
$`\iint _{\Gamma _{🚢}} \phi _t{\bf n}dS`$と$`\iint _{\Gamma _{🚢}}\phi _{t}({\bf x}-{\bf x} _c)\times{\bf n}dS`$がわかればある場所の圧力はわからないが，
🚢にかかる力は計算できるのでそれでも問題ない．

体積積分がゼロとなるように，領域内でラプラス方程式を満たすような$`\varphi`$，
そして$`\Gamma _{🚢}`$上ではこちらが望む$`\varphi _n`$となり，また$`\Gamma \rm other`$上では$`\varphi=0`$となる
そんな$`\varphi`$をBIEを使って計算する．この$`\varphi`$を使うと次の式が成り立つ．
（NOTE：境界上の全ての節点上で$`\varphi`$と$`\varphi _n`$が求まったとする）

```math
\begin{align*}
0 &= \iint _\Gamma {\left( {\varphi\nabla {\phi _t} ({\bf{x}}) - {\phi _t} ({\bf{x}})\nabla \varphi} \right) \cdot {\bf{n}}({\bf{x}})dS}\\
\rightarrow 0 &= \iint _{\Gamma _{🚢}+\Gamma _{🌊}+\Gamma _{\rm wall}} \varphi {\phi _{nt}} dS - \iint _{\Gamma _{🚢}+\Gamma _{🌊}+\Gamma _{\rm wall}} {\phi _t} \varphi _n dS\\
\rightarrow 0 &= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} \varphi {\phi _{nt}} dS - \iint _{\Gamma _{🚢}+\Gamma _{🌊}} {\phi _t} \varphi _n dS\\
\rightarrow \iint _{\Gamma _{🚢}} {\phi _t} \varphi _n dS &= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} \varphi {\phi _{nt}} dS - \iint _{\Gamma _{🌊}} {\phi _t} \varphi _n dS\\
\rightarrow \iint _{\Gamma _{🚢}} \phi _t
\begin{bmatrix}
\boldsymbol{n} \\
(\boldsymbol{x} - \boldsymbol{x} _c) \times \boldsymbol{n}
\end{bmatrix} dS
&= \iint _{\Gamma _{🚢}+\Gamma _{\rm wall}} {\boldsymbol{\varphi} _{1-6}} {\phi _{nt}} dS - \iint _{\Gamma _{🌊}} {\phi _t} {\boldsymbol{\varphi} _n} _{1-6} dS\\
\end{align*}
```

つまり，$`\varphi _n`$を適当に選べば，左辺は知りたかった積分となり，右辺の積分で計算できることになる．

もし浮体がもう一つあると

```math
\begin{align*}
\iint _{\Gamma _{🚢}} \phi _t
\begin{bmatrix}
\boldsymbol{n} \\
(\boldsymbol{x} - \boldsymbol{x} _c) \times \boldsymbol{n}
\end{bmatrix} dS
& = \iint _{\Gamma _{🚢}+\Gamma _{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi} _{1-6}} {\phi _{nt}} dS - \iint _{\Gamma _{🚤}+\Gamma _{🌊}} {\phi _t} {\boldsymbol{\varphi} _n} _{1-6} dS\\
\rightarrow \iint _{\Gamma _{🚢}} \phi _t
\begin{bmatrix}
\boldsymbol{n} \\
(\boldsymbol{x} - \boldsymbol{x} _c) \times \boldsymbol{n}
\end{bmatrix} dS
& = \iint _{\Gamma _{🚢}+\Gamma _{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi} _{1-6}} {\phi _{nt}} dS - \iint _{\Gamma _{🌊}} {\phi _t} {\boldsymbol{\varphi} _n} _{1-6} dS
\end{align*}
```

同じように

```math
\begin{align*}
\iint _{\Gamma _{🚤}} \phi _t
\begin{bmatrix}
\boldsymbol{n} \\
(\boldsymbol{x} - \boldsymbol{x} _c) \times \boldsymbol{n}
\end{bmatrix} dS
& = \iint _{\Gamma _{🚢}+\Gamma _{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi} _{7-12}} {\phi _{nt}} dS - \iint _{\Gamma _{🌊}} {\phi _t} {\boldsymbol{\varphi} _n} _{7-12} dS
\end{align*}
```

$`\iint _{\Gamma _{🚢}+\Gamma _{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi} _{1-6}} {\phi _{nt}} dS`$や
$`\iint _{\Gamma _{🚢}+\Gamma _{🚤}+\Gamma _{\rm wall}} {\boldsymbol{\varphi} _{7-12}} {\phi _{nt}} dS`$
は加速度行列とある既知変数から成る行列の積で表される．こうして，運動方程式の$`\boldsymbol{F} _{\text {hydro }}`$と$`\boldsymbol{T} _{\text {hydro }}`$を加速度によって表すことができ，
運動方程式は加速度だけに関する連立方程式となる．

この方法は，Wu and {Eatock Taylor} (1996)，[Kashiwagi (2000)](http://journals.sagepub.com/doi/10.1243/0954406001523821)，[Wu and Taylor (2003)](www.elsevier.com/locate/oceaneng)で使用されている．
この方法は，複数の浮体を考えていないが，[Feng and Bai (2017)](https://linkinghub.elsevier.com/retrieve/pii/S0889974616300482)はこれを基にして２浮体の場合でも動揺解析を行っている．

[./BEM_solveBVP.hpp#L708](./BEM_solveBVP.hpp#L708)

---
## ⛵ 陽に与えられる境界条件に対して（造波装置など） 

造波理論については，[Dean et al. (1991)](http://books.google.co.uk/books/about/Water_Wave_Mechanics_for_Engineers_and_S.html?id=9-M4U_sfin8C&pgis=1)のp.170に書いてある．

造波板となるobjectに速度を与えることで，造波装置などを模擬することができる．
[強制運動を課す](../../builds/build_bem/main.cpp#L337)

[ここ](../../builds/build_bem/BEM_utilities.hpp#L297)では，Hadzic et al. 2005の造波板の動きを模擬している．
角速度の原点は，板の`COM`としている．

[`setNeumannVelocity`](../../builds/build_bem/BEM_setBoundaryTypes.hpp#L118)で利用され，$\phi _{n}$を計算する．

[./BEM_utilities.hpp#L15](./BEM_utilities.hpp#L15)

### 🪼 フラップ型造波装置 

|   | name   |  description  |
|:-:|:-------:|:-------------:|
| 0 | `flap`|    name       |
| 1 | `start` | start time    |
| 2 | `A`     | wave amplitude|
| 3 | `T`     | wave period   |
| 4 | `h`     | water depth   |
| 5 | `l`     | length from hinge to flap end |
| 6 | `axis`  | x       |
| 7 | `axis`  | y       |
| 8 | `axis`  | z       |

[./BEM_utilities.hpp#L158](./BEM_utilities.hpp#L158)

### 🪼 ピストン型造波装置 

|   | name   |  description  |
|:-:|:-------:|:-------------:|
| 0 | `piston`|    name       |
| 1 | `start` | start time    |
| 2 | `A`     | wave amplitude|
| 3 | `T`     | wave period   |
| 4 | `h`     | water depth   |
| 5 | `axis`  | x       |
| 6 | `axis`  | y       |
| 7 | `axis`  | z       |

ピストン型の造波特性関数：

```math
F(f,h) = \frac{H}{S}=\frac{4\sinh^2(kh)}{2kh+\sinh(2kh)}=\frac{2 (\cosh(2kh) - 1)}{2kh+\sinh(2kh)}
```

$`S`$は造波版のストロークで振幅の２倍である．例えば，振幅が$`A=1`$mの波を発生させたい場合，
$`S = \frac{H}{F}= \frac{2A}{F} = \frac{1}{F(f,h)}`$となり，
これを造波板の変位：$`s(t) = \frac{S}{2} \cos(wt)`$と速度：$`\frac{ds}{dt}(t) = \frac{S}{2} w \sin(wt)`$に与えればよい．(see [Dean et al. (1991)](http://books.google.co.uk/books/about/Water_Wave_Mechanics_for_Engineers_and_S.html?id=9-M4U_sfin8C&pgis=1))

[./BEM_utilities.hpp#L195](./BEM_utilities.hpp#L195)

### 🪼 正弦・余弦（`sin` もしくは `cos`）の運動 

|   | name        |  description  |
|:-:|:-----------:|:-------------:|
| 0 | `sin`/`cos` |    name       |
| 1 | `start`     | start time    |
| 2 | `a`         | amplitude     |
| 3 | `T`         | period        |
| 4 | `axis`      | x             |
| 5 | `axis`      | y             |
| 6 | `axis`      | z             |
| 7 | `axis`      | rotation in x axis  |
| 8 | `axis`      | rotation in y axis  |
| 9 | `axis`      | rotation in z axis  |

名前が$`\cos`$の場合、$`{\bf v}={\rm axis}\, A w \sin(w (t - \text{start}))`$ と計算されます．
名前が$`\sin`$の場合、$`{\bf v}={\rm axis}\, A w \cos(w (t - \text{start}))`$ と計算されます．

[./BEM_utilities.hpp#L243](./BEM_utilities.hpp#L243)

---
## ⛵ その他 

### 🪼 境界値問題の未知変数 

`isNeumannID_BEM`と`isDirichletID_BEM`は，節点と面の組みが，境界値問題の未知変数かどうかを判定する．
多重節点でない場合は，`{p,nullptr}`が変数のキーとなり，多重節点の場合は，`{p,f}`が変数のキーとなる．

[./BEM_utilities.hpp#L570](./BEM_utilities.hpp#L570)

---
### 🪼 エネルギー保存則（計算精度のチェックに利用できる） 

流体全体の運動エネルギーは，ラプラス方程式と発散定理を使うと，次のように境界面に沿った積分で表される．

```math
E _K =\frac{\rho}{2} \iint _\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
```

また，流体の位置エネルギーは，次のように表される．

```math
E _P = \frac{\rho}{2} \iint _\Gamma (0,0,g(z - z _0)^2) \cdot {\bf n} d\Gamma
```

<details>

---

<summary>
💡 なぜか？
</summary>

テンソルを使って考えてみると

```math
\begin{align*}
\nabla \cdot (\phi\nabla\phi) &= \frac{\partial\phi}{\partial x _i} \frac{\partial\phi}{\partial x _i} + \phi \frac{\partial^2\phi}{\partial x _i \partial x _i}\\
&= \nabla \phi \cdot \nabla \phi + \phi \nabla^2 \phi\\
&= \nabla \phi \cdot \nabla \phi
\end{align*}
```

よって，

```math
\iiint _\Omega \nabla\phi\cdot\nabla\phi d\Omega = \iiint _\Omega \nabla \cdot (\phi\nabla\phi) d\Omega = \iint _\Gamma \phi\nabla\phi\cdot {\bf n} d\Gamma
```

---

```math
E _P = \rho g \iiint _\Omega (z - z _0) d\Omega
= \rho g \iiint _\Omega \frac{1}{2} \nabla \cdot (0,0,(z - z _0)^2) d\Omega
= \rho g \iint _\Gamma \frac{1}{2} (0,0,(z - z _0)^2) \cdot {\bf n} d\Gamma
= \frac{1}{2}\rho g \iint _\Gamma (z - z _0)^2 n _z d\Gamma
```

---

</details>

[./BEM_calculateVelocities.hpp#L334](./BEM_calculateVelocities.hpp#L334)

### 🪼 内部流速の計算方法（使わなくてもいい） 

[Fochesato2005](https://onlinelibrary.wiley.com/doi/10.1002/fld.838)にあるように，
流体内部の流速$`\nabla \phi`$は，BIEを微分して求めることができる．

```math
u({\bf a}) = \nabla\phi({\bf a}) = \int _{\partial \Omega} \frac{\partial Q}{\partial n} ({\bf x})Q({\bf x}, {\bf a}) - \phi({\bf x}) \frac{\partial Q}{\partial n} ({\bf x}, {\bf a}) d\Gamma
```

```math
Q({\bf x},{\bf a}) = \frac{{\bf r}}{4\pi r^3}, \quad \frac{\partial Q}{\partial n} ({\bf x},{\bf a}) = \frac{1}{4\pi r^3} (3 \mathbf{n} - (\mathbf{r} \cdot \mathbf{n}) \frac{\mathbf{r}}{r^2})
```

[./BEM_calculateVelocities.hpp#L421](./BEM_calculateVelocities.hpp#L421)

---
### 🪼 JSONファイルの出力 

JSONファイルには，計算結果を出力する．

流体の場合

| 項目 | 詳細|
|---:|:---|
| `simulation_time` | シミュレーション上の時間 |
| `cpu_time` | CPU時間(CPUがプログラムを実行していた時間の合計) |
| `wall_clock_time` | 実時間 |
| `***_volume` | 流体の体積 |
| `***_EK` | 流体の運動エネルギー |
| `***_EP` | 流体の位置エネルギー |
| `***_E` | 流体の全エネルギー |

剛体などで，浮体か`output`に`json`が指定されている場合

| 項目 | 詳細|
|---:|:---|
| `simulation_time` | シミュレーション上の時間 |
| `cpu_time` | CPU時間(CPUがプログラムを実行していた時間の合計) |
| `wall_clock_time` | 実時間 |
| `***_pitch` | 浮体のピッチ角 |
| `***_yaw` | 浮体のヨー角 |
| `***_roll` | 浮体のロール角 |
| `***_force` | 浮体に働く力 |
| `***_torque` | 浮体に働くトルク |
| `***_accel` | 浮体の加速度 |
| `***_velocity` | 浮体の速度 |
| `***_COM` | 浮体の重心位置 |
| `***_area` | 浮体の面積 |
| `***_EK` | 浮体の運動エネルギー |
| `***_EP` | 浮体の位置エネルギー |

[./main.cpp#L615](./main.cpp#L615)

---
# 🐋 実行方法 

## ⛵ ファイルのダウンロード 

上書きされるので注意．ダウンロードしたら，`build_bem`ディレクトリに移動．

```sh
git clone https://github.com/tomoakihirakawa/cpp.git
cd ./cpp/builds/build_bem
```

## ⛵ 入力ファイルの生成． 

```sh
python3 input_generator.py
```

例えば，`./input_files/Hadzic2005`が生成される．

## ⛵ プログラムのコンパイルと実行 

`clean`でCMake関連のファイルを削除して（ゴミがあるかもしれないので），
`cmake`で`Makefile`を生成して，`make`でコンパイルする．

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

実行

```sh
./main ./input_files/Hadzic2005
```

[./main.cpp#L758](./main.cpp#L758)

---
# 🐋 Input Generator 

This file is used to generate the input files for the BEM-MEL.

[./input_generator.py#L1](./input_generator.py#L1)

---
<img src="schematic_Hadzic2005.png" width="400px"/>

This case based on [Had{\v{z}}i{\'{c}} et al. (2005)](https://linkinghub.elsevier.com/retrieve/pii/S0307904X05000417) is for the validation of the floating body motion analysis using the BEM-MEL.        
The floating body is a rectangular box with the dimension of L10 cm x H5 cm x W29 cm.        
The density of the floating body is 0.68x1000 kg/m^3, therefore the mass of the floating body is 0.68x0.05x0.1x0.29x1000 kg.
The moment of inertia of the floating body is 14 kg cm^2.

[CAD data](https://a360.co/46CisV7)

[spheric Test 12](https://www.spheric-sph.org/tests/test-12)

[Youtube Nextflow](https://www.youtube.com/watch?v=H92xupH9508)

[./input_generator.py#L244](./input_generator.py#L244)

---
<img src="schematic_Ren2015.png" width="400px" />

This case based on [Ren et al. (2015)](https://linkinghub.elsevier.com/retrieve/pii/S0141118714001175) is for the validation of the floating body motion analysis using the BEM-MEL.
The floating body is a rectangular box with the dimension of $`(l _x,l _y,l _z)=(0.3,0.42,0.2) {\rm m}`$ 
The density of the floating body is $`0.5\times1000 {\rm kg/m^3}`$.
The moment of inertia of the floating body is $`(I _{xx},I _{yy},I _{zz}) = (\frac{m}{12}(l _y^2+l _z^2),\frac{m}{12}(l _x^2+l _z^2),\frac{m}{12}(l _x^2+l _y^2))`$.

You can find numerical results compared with this case from Cheng and Lin (2018) and \cite{Bihs2017}.

[Youtube DualSPHysics](https://www.youtube.com/watch?v=VDa4zcMDjJA)

[./input_generator.py#L123](./input_generator.py#L123)

---
This case is for the validation of the floating body motion analysis using the BEM-MEL.

<img src="schematic_Kramer2021.png" width="400px" />

The floating body is a sphere with the diameter of 0.3 m.
The mass of the floating body is 7.056 kg.
The moment of inertia of the floating body is set to be almost infinite to ignore the effect of the rotation.

The sphere is dropped from the height of 0.03 m above the water surface.

[./input_generator.py#L322](./input_generator.py#L322)

---
# 🐋 Examples 

**[See the Examples here!](EXAMPLES.md)**

[./main.cpp#L798](./main.cpp#L798)

---
