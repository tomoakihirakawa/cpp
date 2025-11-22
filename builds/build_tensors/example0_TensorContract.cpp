
#include "basic_arithmetic_array_operations.hpp"

/*DOC_EXTRACT

(
   cd builds/build_tensors
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)


```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example0_TensorContract.cpp
make
./example0_TensorContract
```

## テンソル積 `TensorProduct`

テンソル積$(a_i b_j)$は，記号$\otimes$を使って$\boldsymbol{a} \otimes \boldsymbol{b}$と表される(p.16)．
$a_i b_j c_k$は，$\boldsymbol{a} \otimes \boldsymbol{b} \otimes \boldsymbol{c}$．
かけるテンソルの階数が異なってもよい．積の結果，階数は添字の種類の数になる．

テンソル積$(a_i b_j)$と$(b_j a_i)$は同じものであるが，$(a_j b_i)$は異なるものである．
$i$と$j$を入れ替えても，
総和を取る場合は結果は同じになるかもしれないが，
他のテンソルとの計算結果が変わってしまう可能性がある．
なので，計算に関わるテンソル間で，添字は一貫していなければならない．

$k$が$ｊ$との場合，$a_i b_j c_j$の場合，
これはベクトルを使って，$(\boldsymbol{a} \otimes \boldsymbol{b}) \cdot \boldsymbol{c}$と表せる．
$a_i b_j c_i$なら，$\boldsymbol{c} \cdot (\boldsymbol{a} \otimes \boldsymbol{b})$と表せる．

```Mathematica
a = {1, 2, 3};
b = {4, 5, 6};
c = {7, 8, 9};
Dot[TensorProduct[a, b], c]
Dot[c, TensorProduct[a, b]]
{122, 244, 366}
{200, 250, 300}
```

なので，$Dot[TensorProduct[a, b], c]$は，$b$と$c$の内積が計算されなければならない．
また，添字の表記を眺めてもわかるように，
$(\boldsymbol{a} \otimes \boldsymbol{b}) \cdot \boldsymbol{c}
= (\boldsymbol{a} \otimes (\boldsymbol{b}) \cdot \boldsymbol{c})
= (\boldsymbol{a} (\boldsymbol{b}) \cdot \boldsymbol{c})$
$

```Mathematica
TensorProduct[a, Dot[b, c]]
a*Dot[b, c]
{122, 244, 366}
{122, 244, 366}
```

## テンソルの縮約 `TensorContract`

$a_i b_j c_j$は，
テンソル積をとって，$j$について総和を取る操作，とみれるし，
${\boldsymbol{b}} \cdot {\boldsymbol{c}}$を計算して，${\boldsymbol{a}}$にスカラーをかける操作ともみれる．
さらに，
遠回りな方法として，まず一度${\boldsymbol{a}} \otimes {\boldsymbol{b}} \otimes {\boldsymbol{c}}$を計算して，
その後，$j$について総和を取る操作ともみれる．最終結果が，２階のテンソルであるとわかるにも関わらず，
最後の計算は，一度３階まで次元が増えてしまうため，無駄が多い．
ただ，最終結果がわからない場合もあるだろう．このような計算は，`TensorContract`という関数で行うことができる．

```Mathematica
TensorContract[TensorProduct[a, b, c], {{2, 3}}]
TensorContract[TensorProduct[a, b, c], {{1, 3}}]
{122, 244, 366}
{218, 175, 132}
```

`{{2, 3}}`は，$a_ib_jc_k$の$ijk$という添字の2番目と3番目を総和を取ることを意味する．
ちなみに，２番目の計算からわかるように，
始めに計算した`Dot[c, TensorProduct[a, b]]`
は，`TensorContract[TensorProduct[a, b, c], {{1, 3}}]`と同じ結果であり，つまり，$a_ib_jc_i$の計算をしていたことになる．

## プログラムにおけるloopとの関係

テンソルを使った計算は，loopを使った計算と同じようなものだが，
違う点は，

*/

int main() {

   std::array<double, 3> A = {1, 2, 3};
   std::array<double, 3> B = {4, 5, 6};

   std::cout << Dot(A, B) << std::endl;
   std::cout << TensorProduct(A, B) << std::endl;
   std::cout << TensorProduct(A, TensorProduct(A, B)) << std::endl;
}