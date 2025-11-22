/*DOC_EXTRACT sparse_crs_demo

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test2_CRS.cpp
make; ./test2_CRS
```

目的
- 疎行列 A とベクトル b に対して、CRS 行オブジェクト（各行が隣接（列）ポインタと係数を持つグラフ表現）で構築した A の行列-ベクトル積が、密行列の計算 Dot(A, b) と一致することをデモし、CRS ベース実装の利点を確認する。

何をしているか
- 4x4 の行列 A とベクトル b を用意。
- 行ごとに CRS ノード（CRS*）を生成し、CRS::value に b[i] を格納（b をノード値として持たせる）。
- 全ペア (i, j) を走査して A[i][j] ≠ 0 のとき crs_i->increment(crs_j, A[i][j]) で隣接と係数を登録（行 i の非零列 j と値 A[i][j] をエッジとして張る）。
- Dot(V_CRS_vec) を呼ぶと、各行 i で selfDot() = Σ_j A[i][j] · value_j = (A·b)[i] を並列で計算する。
- Dot(A, b) は密行列版の計算。両者の結果が一致することを確認する。

出力の意味
- "Dot : {35,15,15,35}" は A·b の結果（各行の線形結合）。二行とも同じ値が出るのは、前者が CRS ベース、後者が密行列ベースで、同一の A·b を計算しているため。

メリット（CRS ベースの利点）
- 計算量: O(nnz) の行列-ベクトル積。密行列の O(n^2) に比べ、疎な場合に有利。
- 並列化: Dot は std::execution::par_unseq を用いて行ごとに並列化。スレッド/ベクトル化の両方を活かせる。
- メモリ局所性/アクセス削減: PreparedDot は各 CRS に入力ベクトル要素へのポインタを一時保持し、ランダムアクセスを削減してスループットを向上。
- 動的組立て: increment により行の非零パターンをインクリメンタルに構築でき、有限要素/グラフからの組立に向く。
- インデックス安全性: 計算は ret[crs->__index__] へ書き出すため、コンテナ順序（unordered_set など）に依存しない。

使用上の注意
- Dot(const std::vector<T*>&) を使うため、unordered_set から ToVector で変換している。
- CRS のライフサイクル管理（new/delete）は実コードではスマートポインタ等で安全に扱うことを推奨。
- __index__ は 0..n-1 の連番で設定すること。

実行方法（例）
```sh
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test2_CRS.cpp
make; ./test2_CRS
```

期待される出力（要旨）
```
 Dot : {35,15,15,35}
 Dot : {35,15,15,35}
```

この先へ
- GMRES/Arnoldi の A·v を、CRS ベースの Dot/PreparedDot に差し替えることで大規模疎問題に対応。プリコンディショナ（例: ILU）とも親和性が高い。
*/

#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

int main() {

  /* -------------------------------------------------------------------------- */
  const VV_d A = {{4., -1, 0, -1}, {-1., 4, -1, 0}, {0., -1, 4, -1}, {-1., 0, -1, 4}};
  const V_d b = {15., 10., 10, 15};
  // 初期値が大事
  V_d x0(b.size(), 0.);
  V_d ans = {6.875, 5.625, 5.625, 6.875};
  /* -------------------------------------------------------------------------- */

  std::unordered_set<CRS *> V_CRS;

  {
    auto crs = new CRS();
    crs->setIndexCRS(0);
    crs->value = 15.;
    V_CRS.emplace(crs);
  }
  {
    auto crs = new CRS();
    crs->setIndexCRS(1);
    crs->value = 10.;
    V_CRS.emplace(crs);
  }
  {
    auto crs = new CRS();
    crs->setIndexCRS(2);
    crs->value = 10.;
    V_CRS.emplace(crs);
  }
  {
    auto crs = new CRS();
    crs->setIndexCRS(3);
    crs->value = 15.;
    V_CRS.emplace(crs);
  }
  for (const auto &crs0 : V_CRS)
    for (const auto &crs1 : V_CRS) {
      int i = crs0->getIndexCRS();
      int j = crs1->getIndexCRS();
      if (A[i][j] != static_cast<double>(0.))
        crs0->increment(crs1, A[i][j]);
    }

  // std::cout << " Dot : " << Dot(V_CRS, V_CRS) << std::endl;
  // std::cout << " Dot : " << Dot(V_CRS) << std::endl;

  // unordered_set -> vector に変換してから Dot
  auto V_CRS_vec = ToVector(V_CRS);
  std::cout << " Dot : " << Dot(V_CRS_vec) << std::endl;
  std::cout << " Dot : " << Dot(A, b) << std::endl;
};
