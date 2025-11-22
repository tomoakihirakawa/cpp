#include <array>
#include <memory>

bool _PSEUDO_QUADRATIC_ELEMENT_ = false;

#include "/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_bem/BEM_setBoundaryTypes.hpp"
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT 0_0_0_translation_of_a_multipole_expansion

# 多重極展開の実装

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Debug ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion_with_tree_20240818.cpp
make
./test_translation_of_a_multipole_expansion_with_tree_20240818 ./pumpkin.obj
paraview check_M2L.pvsm
```

*/

/*

# Fast Multipole Method

## pole class

pole class has the following attributes:

- position
- weights
- normal vector
- updater function (to update the intensity, that is the potential, of the pole)

## Buckets class

Buckets class stores specified objects as `Buckets<T>`, and generates tree structure until the number of objects in a bucket is less than or equal to the specified number of objects per bucket.

The step to generate the tree structure should be as follows:

1. add objects to the bucket
2. set the maximum level of the tree using `setLevel`
3. generate the tree structure using `generateTree` while specifying the condition to stop the generation of the tree structure

# Fast Multipole Method

The Fast Multipole Method (FMM) is an algorithm for the efficient calculation of the integration of the pole/potential using the tree structure, the multipole expansion, shifting expansion, and the local expansion. Since FMM calculates integration/summation, such as BIE and does not make the coefficient matrix, solver for the simultaneous linear equations should be iterative methods. GMRES is commonly used for the solver with FMM.

| First steps | GRMES iterative step | description | | |
| --- | --- | --- | --- | --- |
| 1 | | add poles to the root bucket | | |
| 2 | | generate the tree structure from the root bucket | | |
| 3 (before M2M) | | expansion of the poles | | |
| 4 | 1 | **update the intensity of the poles** | | |
| 5 | 2 | Multipole to Multipole (M2M): shift the multipole expansion at each center, from the deeper level to the upper level | about 8 🪣 -> 1 parent 🪣 | use pre-computed SPH |
| 6 | 3 |  Multipole to Local (M2L)| every 🪣 -> (only same level) -> many local 🪣 | use pre-computed SPH |
| 7 | 4 | Local to Local (L2L) | 1 🪣 -> about 8 children 🪣 | use pre-computed SPH |
| 8 | 5 | Add direct integration for the near field and the integration using the local expansion for the far field | | |

Many part of process are dependent on relative position of the poles and the buckets. Therefore, many part of the first steps are saved and reused in the following iterative steps. Remaining part for iterative steps are the update of the intensity of the poles, and simple incrementatation in four-fold for-loops. However, the number of incrementation is not negligible, and the direct integration for the near field also takes time. FMM is surely faster than the direct summation when the number of poles is more than about 10000, but the calculation time is already long when the number of poles is about 10000.

## 要素法特有の話

そもそものガウス点がすくなければツリーを伸ばしていけて，直接積分の量を減らせるため，早くなる
直接積分の部分のガウス点は減らしたくない．これを両立することが大事だ．

このような議論はされていない．

\insert{Multipole_Expansion}

*/

#define _BEM_

#include "Network.hpp"
#include "vtkWriter.hpp"

// auto obj = std::make_unique<Network>("./bunny.obj");
// auto obj = std::make_unique<Network>("./pumpkin.obj");
// auto obj = std::make_unique<Network>("./HumanFace.obj");
// auto obj = std::make_unique<Network>("./HumanBrain.obj");
// auto obj = std::make_unique<Network>("./Armadillo.obj");

// Function to write polygons to a file
void writePolygonToFile(const std::string& name, auto polygon) {
   std::ofstream ofs(name);
   vtkPolygonWrite(ofs, polygon);
   ofs.close();
}

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[]) {

   if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
      return 1;
   }

   auto obj = std::make_unique<Network>(std::string(argv[1]));

   /* -------------------------------------------------------------------------- */
   /*                  境界面の読み込み．境界面を原点に移動する．                       */
   /* -------------------------------------------------------------------------- */

   std::array<double, 3> center = {0., 0., 0.};  //{5., 5., 5.};
   double num = 0.;
   for (const auto& p : obj->getPoints()) {
      center += p->X;
      num += 1.;
   }
   center /= num;
   obj->translate(-center);
   // obj->scale(10.);

   std::ofstream ofs("./bunny_obj.vtp");
   vtkPolygonWrite(ofs, obj->getFaces());
   ofs.close();

   TimeWatch tw;

   /*

   1. add poles to the bucket
   2. generate tree
   3. multipole expansion

   */

   /* -------------- BEM ----------------- */
#ifdef _BEM_
   obj->setGeometricProperties();
   for (const auto& f : obj->getFaces())
      f->setIntegrationInfo();
   setBoundaryTypes(obj.get(), {});
   setPhiPhinOnFace(obj.get());
   auto size = setNodeFaceIndices(obj.get());
   std::cout << "size = " << size << std::endl;
#endif
   /* ------------------------------------- */

   //@ -------------------------------------------------------------------------- */
   //@                       バケットの作成．極の追加 add                             */
   //@ -------------------------------------------------------------------------- */

   //! 実際に問題を解く場合は，計算の流れにそって，phiOnFaceやphinOnFaceを更新する．
   //! ここでは，適当な値を入れている．
   auto targets_points = obj->getPoints();
   for (const auto& p : targets_points) {
      for (const auto& [f, i] : p->f2Index) {
         p->phiphin[0] = p->phiOnFace[f] = RandomReal({-1, 1});
         p->phiphin[1] = p->phinOnFace[f] = RandomReal({-1, 1});
      }
   }

   /* -------------------------------------------------------------------------- */

   auto bounds = obj->getUniformBounds(0.6);
   Buckets<std::shared_ptr<source4FMM<target4FMM>>> B_poles(bounds, (bounds[0][1] - bounds[0][0]) / 3.);
   std::cout << "バケットの作成．極の追加" << std::endl;
   /* phiOnFaceやphinOnFaceを更新することで，get_valuesの結果は自動的に更新される */
   for (auto& F : obj->getSurfaces()) {
      auto [p0, p1, p2] = F->getPoints();
      auto closest_p_to_origin = p0;
      auto X012 = ToX(F->getPoints());
      auto cross = Cross(X012[1] - X012[0], X012[2] - X012[0]);
      auto key0 = std::get<1>(pf2ID(p0, F));
      auto key1 = std::get<1>(pf2ID(p1, F));
      auto key2 = std::get<1>(pf2ID(p2, F));
      std::array<double*, 2> pair_pointer_phiphin0 = {&p0->phiOnFace.at(key0), &p0->phinOnFace.at(key0)};
      std::array<double*, 2> pair_pointer_phiphin1 = {&p1->phiOnFace.at(key1), &p1->phinOnFace.at(key1)};
      std::array<double*, 2> pair_pointer_phiphin2 = {&p2->phiOnFace.at(key2), &p2->phinOnFace.at(key2)};
      for (const auto& [t0t1, quadrature_weight, shape3, X, cross, J_det] : F->map_Point_LinearIntegrationInfo_vector[0].at(closest_p_to_origin)) {

         /*DOC_EXTRACT 0_0_0_translation_of_a_multipole_expansion

         ## 最適な実装を考える

         Fast Multipole Method（FMM）は，**反復解法（例：GMRES）**と組み合わせて使われることが多い．
         FMMも工夫して実装しなければ，大幅な計算時間の短縮は難しい．

         実装のポイント・設計指針

         1. **強度更新の効率化**
            各極の「強度」（$\phi$や$\partial\phi/\partial n$×重み）は毎回計算し直す必要がある．
            ポインタで値を参照する仕組みにし，毎回`source4FMM`オブジェクトを作り直さずに済むようにする．
            具体的には**「ラムダ関数（クロージャ）」で強度計算ロジックを保持し，必要なときだけ呼び出して最新の値を得る**．

         2. **遠方場（FMM）・近接場（直接積分）の効率化**
            M2L（Multipole-to-Local）などコストの高い計算を効率化（例えばFFTによる高速畳み込みなど）．
            近接場の計算も**「影響係数」（積分カーネル）を事前計算**し，ループ中はその結果と最新の強度を掛けるだけにする．

         3. **ソースの定義・拡張性**
            物理問題ごとに異なるソース（$\phi$, $\partial\phi/\partial n$ など）を柔軟に扱えるように，ラムダ関数で初期化できる形式にしておく．
            これにより，**コード全体の柔軟性と拡張性が高まる**．

         4. **BEM独自の要素（リジッドモード対応など）**
            リジッドモード法を使いやすい形で，ミスが起こりにくいように分かりやすく，組み込みやすい形でロジックを注入できるように設計．

         ``` mermaid
         classDiagram
            direction RL
            %% --- クラス定義 ---
            class source4FMM {
               +X: 位置
               +normal: 法線
               +weighted_source_densities: get_weighted_source_densitiesで更新
               +get_weighted_source_densities: ラムダ関数
               +use_this_source_when_set_direct_integration: ラムダ関数
               +update()
            }
            class target4FMM {
               +Xtarget: 位置
               +setL2P()
               +setDirectIntegration()
               +integrate()
               +integrateFMM()
               +integrate(Buckets)
            }
            class Buckets {
               +struct Moments
               +MomentsMultipoleExpansion : Moments
               +MomentsLocalExpansion : Moments
               +all_stored_objects_vector: source4FMM
               +buckets_for_M2L
               +buckets_near
               +add()
               +generateTree()
               +setLevel()
               +forEachAtLevel()
               +forEachAtDeepest()
            }
            class GMRES {
               +dot_product()
            }

            %% --- 関係 ---
            Buckets o-- source4FMM : (極, n)
            target4FMM ..> Buckets : 参照
            GMRES ..> target4FMM : dot_product
            GMRES ..> potential : 更新
            source4FMM ..> potential : get_weighted_source_densities

         ```

         `get_weighted_source_densities`は，更新された強度を取得し，`source4FMM`の
         `weighted_source_densities`に保存するためのラムダ関数である．

         三角要素１つに対して，ソース１つを保存することにした．
         三角要素とソースは，1対1の関係にする．これはtarget点から十分離れた要素では大きな問題にはならない．
         targetがソース点を参照した際に，十分離れていない場合は，`use_this_source_when_set_direct_integration`を使って要素を細かく分割し積分する．

         - `get_weighted_source_densities()`:　`source4FMM::weighted_source_densities`を更新するためのラムダ関数
         - `use_this_source_when_set_direct_integration(const target4FMM* origin)`:　ラムダ関数

         ### M2L（２D畳み込み）

         \ref{Moments::set_m2l}{この関数}で，`Moments::m2l_cache`を作成する．

         `set_m2l()`関数は，何度も更新される**多重極モーメント**（Multipole Moments）と，定数の**変換係数`m2lFunction()`**をキャッシュ (m2l_cache) へ保存する関数．
         両方ともconstantであり（モーメントはポインタとしてはconst.），繰り返される積和の計算において更新する必要がない．

         $$
         {L}_{j}^{k}=\mathop{\sum}\limits_{{n}={0}}\limits^{N}{\mathop{\sum}\limits_{{m}=-{n}}\limits^{n}{\left\{{\left\{{\mathop{\sum}\limits_{{k}_{\blacksquare}}\limits^{}{{\left\{{{O}_{n}^{m}}\right\}}_{{k}_{\blacksquare}}}}\right\}\frac{{\mathcal{A}}_{j}^{k}{\mathcal{A}}_{n}^{m}{Y}_{{n}+{j}}^{{m}-{k}}\left({\mathit{\alpha}{,}\mathit{\beta}}\right)}{{\left({-{1}}\right)}^{n}{\mathcal{A}}_{{n}+{j}}^{{m}-{k}}{\mathit{\rho}}^{{n}+{j}+{1}}}}\right\}}}
         \label{eq:m2l}
         $$

         `m2l_cache`は，次のような形で保存される．`{*Ljk, {{Ajkmn, *MMjkmn}, ...}}`のような形で保存される：

         ```
         m2l_cache =
         {*L00,{{A0000, *MM0000},{A0001, *MM0001},{A0010, *MM0010},{A0011, *MM0011}}},
         {*L01,{{A0100, *MM0100},{A0101, *MM0101},{A0110, *MM0110},{A0111, *MM0111}}},
         {*L10,{{A1000, *MM1000},{A1001, *MM1001},{A1010, *MM1010},{A1011, *MM1011}}},
         {*L11,{{A1100, *MM1100},{A1101, *MM1101},{A1110, *MM1110},{A1111, *MM1111}}},
         ...}
         ```

         ### 離散フーリエ変換を使ったM2L（２D畳み込み）の高速化

         \eqref{eq:m2l}は，次のように変形できる．

         $$
         \begin{align}
         {L}_{j}^{k}&=\mathop{\sum}\limits_{{n}={0}}\limits^{N}{\mathop{\sum}\limits_{{m}=-{n}}\limits^{n}{\left\{{\left\{{\mathop{\sum}\limits_{{k}_{\blacksquare}}\limits^{}{{\left\{{{O}_{n}^{m}}\right\}}_{{k}_{\blacksquare}}}}\right\}\frac{{\mathcal{A}}_{j}^{k}{\mathcal{A}}_{n}^{m}}{{\left({-{1}}\right)}^{n}{\mathcal{A}}_{{n}+{j}}^{{m}-{k}}}\frac{{Y}_{{n}+{j}}^{{m}-{k}}\left({\mathit{\alpha}{,}\mathit{\beta}}\right)}{{\mathit{\rho}}^{{n}+{j}+{1}}}}\right\}}}\\
                    &={\mathcal{A}}_{j}^{k}\mathop{\sum}\limits_{{n}={0}}\limits^{N}{\mathop{\sum}\limits_{{m}=-{n}}\limits^{n}{\left\{{\left\{{\mathop{\sum}\limits_{{k}_{\blacksquare}}\limits^{}{{\left\{{\frac{{\mathcal{A}}_{n}^{m}}{{\left({-{1}}\right)}^{n}}{O}_{n}^{m}}\right\}}_{{k}_{\blacksquare}}}}\right\}\frac{{Y}_{{n}+{j}}^{{m}-{k}}\left({\mathit{\alpha}{,}\mathit{\beta}}\right)}{{\mathcal{A}}_{{n}+{j}}^{{m}-{k}}{\mathit{\rho}}^{{n}+{j}+{1}}}}\right\}}}\\
                    &={\mathcal{A}}_{j}^{k}\mathop{\sum}\limits_{{n}={0}}\limits^{N}{\mathop{\sum}\limits_{{m}=-{n}}\limits^{n}{\left\{{\left\{{\mathop{\sum}\limits_{{k}_{\blacksquare}}\limits^{}{{\left\{{{\mathcal{O}}_{n}^{m}}\right\}}_{{k}_{\blacksquare}}}}\right\}{\mathcal{Y}}_{{n}+{j}}^{-\left({{m}-{k}}\right)}}\right\}}}\\
                    &={\mathcal{A}}_{j}^{k}{\left\{{{\mathcal{F}}_{2D}^{-{1}}\left[{\left\{{\mathop{\sum}\limits_{{k}_{\blacksquare}}\limits^{}{{\left\{{{\mathcal{F}}_{2D}\left[{{\mathcal{O}}_{n}^{m}}\right]}\right\}}_{{k}_{\blacksquare}}}}\right\}\odot{\mathcal{F}}_{2D}\left[{{\mathcal{Y}}_{q}^{p}}\right]}\right]}\right\}}_{j}^{k}
         \end{align}
         \label{eq:m2l_fft}
         $$

         where

         $$
         \begin{gathered}
         {A}_{j}^{k}=\frac{{\left({-{1}}\right)}^{k}}{\sqrt{\left({{j}-{k}}\right){!}\left({{j}+{k}}\right){!}}},\quad{{\mathcal{A}}_{j}^{k}={i}^{-\left|{k}\right|}{A}_{j}^{k}},\quad
         {{\mathcal{O}}_{n}^{m}={\left({-{1}}\right)}^{n}{\mathcal{A}}_{n}^{m}{O}_{n}^{m}},\quad{{\mathcal{Y}}_{q}^{p}=\frac{1}{{\mathcal{A}}_{q}^{p}}\frac{{Y}_{q}^{-{p}}\left({\mathit{\alpha}{,}\mathit{\beta}}\right)}{{\mathit{\rho}}^{{q}+{1}}}}
         \end{gathered}
         $$

         注意：DFTによる畳み込み積分または相互相関の計算では，パッディングが必要．相互相関の場合は，相互相関のシフト方向にデータを反転させた上でパッディングを行う必要がある．
         各バケツにおけるモーメントの計算後，パッディングや反転や累積の処理は，\ref{Convolver2D}{`Convolver2D`}クラスに任せる．

         \eqref{eq:m2l_fft}から，

         * ${\mathcal{F}}_{2D}\left[{{\mathcal{O}}_{n}^{m}}\right]$はバケツ毎に独立して計算できること（計算は，ツリー全体で１度だけで高速）．
         * そして，あるlocalとなるバケツにおいて，対象となるバケツの${\mathcal{F}}_{2D}\left[{{\mathcal{O}}_{n}^{m}}\right]$を足し合わせること．$O(k_{buckets}N^2)$で軽量．
         * 各バケツにおいて，${\mathcal{F}}_{2D}\left[{{\mathcal{O}}_{n}^{m}}\right]$と${\mathcal{F}}_{2D}\left[{{\mathcal{Y}}_{q}^{p}}\right]$を準備しておく．

         非常に高速に計算できることがわかる．

         ### 準備`set_m2l(const std::vector<TYPE> &buckets)`

         * $\mathcal{A}_{j}^{k}$は，`include/i_absk_A_FMM.hpp`に定義されている．
         * $\mathcal{O}_n^m$は，ME,M2M,M2M...と段階的に計算してきたモーメントに係数を掛けたものであり，
         * $\phi$,$\phi_n$の係数それぞれに対して存在する複素行列．
         * $\mathcal{Y}_q^p$は，$\phi$と$\phi_n$のそれぞれに対して共通の係数行列であり，バケツ毎に計算される．

         \ref{Onm0_Onm1_Yqp}{}で`Yqp`, `Onm0`, `Onm1`を計算し保存する．

         次に，`Fourier2D<std::complex<double>>`へ`Onm0`, `Onm1`, `Yqp`を変換し保存する．\ref{DFT2D_Onm0_Onm1_Yqp}{}

         DFT/FFTを使った畳み込みは，\ref{Fourier2D}{`Convolver2D`}クラスに任せる．

         ```cpp
         Convolver2D<std::complex<double>> convolver2D;
         convolver2D.reset(img_shift.size(), img_shift[0].size(), img_kernel.size(), img_kernel[0].size());
         convolver2D.addKernel(img_kernel);
         convolver2D.addShift(img_shift, {true, false});
         convolver2D.convolve();
         ```

         FMMへの応用を考えた，\ref{Fourier2D}{`Fourier2D`}クラスの設計のポイント：

         重要な処理は，パッディングと反転と累積の処理．
         Convolver2Dの外部でFFT済みの行列が作成されていることを前提とすべき．
         FFTの前に決定しておくしかない処理は，パッディングと反転の処理

         * Convolver2DKernelクラス
         * Convolver2DShiftクラス:係数行列を読み込み，パッディングや反転処理を行う，Convolver2DShift同士は足し合わせることができる

         を作成し，それらを使ってconvolveを計算するようにする．Convolver2DShiftは係数行列を読み込み，パッディングや反転処理を行う．また，
         Convolver2DShift同士は足し合わせることができるようにする．

         ### 実行`m2l()`



         */

         auto get_weighted_source_densities = [shape3,
                                               W = quadrature_weight * J_det,
                                               pair_pointer_phiphin0,
                                               pair_pointer_phiphin1,
                                               pair_pointer_phiphin2]() -> std::array<double, 2> { return {W * Dot(shape3, Tddd{*pair_pointer_phiphin0[0], *pair_pointer_phiphin1[0], *pair_pointer_phiphin2[0]}),
                                                                                                           W * Dot(shape3, Tddd{*pair_pointer_phiphin0[1], *pair_pointer_phiphin1[1], *pair_pointer_phiphin2[1]})}; };

         auto use_this_soruce_when_set_direct_integration = [p0, p1, p2, shape3, cross, X012,
                                                             W = quadrature_weight * J_det,
                                                             pair_pointer_phiphin0,
                                                             pair_pointer_phiphin1,
                                                             pair_pointer_phiphin2](const target4FMM* origin) -> std::vector<std::tuple<std::array<double*, 2>, std::array<double, 2>>> {
            std::array<double, 3> cross = Cross(X012[1] - X012[0], X012[2] - X012[0]), WGN = {0., 0., 0.}, WGnN = {0., 0., 0.}, N012_geometry, R;
            double norm_cross = Norm(cross);
            double tmp, nr;
            for (const auto& [t0, t1, ww] : __array_GW1xGW1__) {
               N012_geometry = ModTriShape<3>(t0, t1);
               nr = Norm(R = (Dot(N012_geometry, X012) - origin->Xtarget));
               tmp = (ww * (1. - t0) / nr);
               WGN += (norm_cross * tmp) * N012_geometry;
               WGnN += (-Dot(R / nr, cross) * tmp / nr) * N012_geometry;
            }
            if (p0 == origin)
               std::get<0>(WGnN) = 0.;  //! リジッドモードテクニック
            if (p1 == origin)
               std::get<1>(WGnN) = 0.;  //! リジッドモードテクニック
            if (p2 == origin)
               std::get<2>(WGnN) = 0.;  //! リジッドモードテクニック
            return {{pair_pointer_phiphin0, std::array<double, 2>{WGN[0], WGnN[0]}},
                    {pair_pointer_phiphin1, std::array<double, 2>{WGN[1], WGnN[1]}},
                    {pair_pointer_phiphin2, std::array<double, 2>{WGN[2], WGnN[2]}}};
         };

         //! ソースの追加

         auto pole = std::make_shared<source4FMM<target4FMM>>(X, F->normal,
                                                              get_weighted_source_densities,
                                                              use_this_soruce_when_set_direct_integration);
         B_poles.add(X, pole);
      }
   }
   std::cout << Magenta << "Add poles" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   //^ -------------------------------------------------------------------------- */
   //^                                ツリー構造を生成                               */
   //^ -------------------------------------------------------------------------- */

   std::cout << "ツリー構造を生成" << std::endl;
   int max_level = 5;
   B_poles.setLevel(0, max_level);
   B_poles.generateTree([](auto bucket) {
      if (bucket->all_stored_objects_vector.empty())
         return false;
      else
         return bucket->all_stored_objects_vector.size() > 500 && bucket->level < bucket->max_level;
   });
   std::cout << Magenta << "Tree" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   // show info of tree
   for (auto i = 0; i < B_poles.level_buckets.size(); ++i) {
      int mean_M2L_size = 0;
      for (auto m2l : B_poles.level_buckets[i])
         mean_M2L_size += m2l->buckets_for_M2L.size();
      mean_M2L_size /= B_poles.level_buckets[i].size();
      std::cout << "level = " << i << ", size = " << B_poles.level_buckets[i].size() << ", mean M2L size = " << mean_M2L_size << std::endl;
   }

   //@ -------------------------------------------------------------------------- */
   //@                                   FMM                                      */
   //@ -------------------------------------------------------------------------- */

   //! 繰り返し計算によって問題が生じないかもチェック
   initializeFMM(B_poles, targets_points);

   for (auto i = 0; i < 1; ++i) {
      updateFMM(B_poles);
      TimeWatch twFMM2;
#pragma omp parallel
      for (auto& p : targets_points)
#pragma omp single nowait
      {
         for (const auto& [f, i] : p->f2Index) {
            // auto [IgPhin_IgnPhi_near, IgPhin_IgnPhi_far] = integrate(B_poles, p->X, eps);
            // std::get<1>(IgPhin_IgnPhi_near) -= p->getSolidAngle() * p->phiOnFace.at(f);  //! ここにphiOnFaceは使わない独立させる
            p->wGPhin_wGnPhi_near = p->integrate();
            p->wGPhin_wGnPhi_far = p->integrateFMM();
            p->wGPhin_wGnPhi_FMM = p->wGPhin_wGnPhi_near + p->wGPhin_wGnPhi_far;
            p->diagonal_coefficient = p->wGPhin_wGnPhi_FMM[0] - p->wGPhin_wGnPhi_FMM[1];
         }
      }

      std::cout << Magenta << "integration including L2P" << Green << ", Elapsed time : " << twFMM2() << colorReset << std::endl;
   }

   //! -------------------------------------------------------------------------- */
   //!                                   直接積分                                   */
   //! -------------------------------------------------------------------------- */

   if (true) {
      TimeWatch twDirect;
      std::cout << "Direct Integration..." << std::endl;
#pragma omp parallel
      for (auto& p : obj->getPoints())
#pragma omp single nowait
      {
         p->igign = p->integrate(B_poles);
      }

      std::cout << Magenta << "Direct Integration" << Green << ", Elapsed time : " << twDirect() << colorReset << std::endl;
   }

   tw.reset();

   // % -------------------------------------------------------------------------- */
   // %                                     出力                                    */
   // % -------------------------------------------------------------------------- */
   {

      //! バケツの可視化のための出力

      std::vector<T8Tddd> cube_level, cube_level_deepest, cube_near;
      std::vector<std::vector<T8Tddd>> cube_M2L;
      std::vector<std::vector<Tddd>> poles_in_bucket;
      std::vector<std::vector<T2Tddd>> line_M2L;
      PVDWriter cubePVD("./output/cubes.pvd");
      PVDWriter lineM2LPVD("./output/line_M2L.pvd");

      std::cout << "output" << std::endl;

      {
         B_poles.forEachAtDeepest([&](Buckets<std::shared_ptr<source4FMM<target4FMM>>>* B) {
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            cube_level_deepest.push_back(t8tddd);
         });
         std::ofstream ofs("./output/cube_level_deepest.vtp");
         vtkPolygonWrite(ofs, cube_level_deepest);
         ofs.close();
      }

      std::cout << "paraview ./output/cube_level_deepest.vtp" << std::endl;

      //! Oはランダムに選ぶのではなく，例としてふさわしい，最も深いレベルの中にある極を選ぶ
      Tddd sample_Origin;
      int sample_max_level = 0;

      for (int i = 0; i < B_poles.max_level; ++i) {
         cube_level.clear();
         B_poles.forEachAtLevel(i, [&](Buckets<std::shared_ptr<source4FMM<target4FMM>>>* B) {
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            cube_level.push_back(t8tddd);
            if (sample_max_level < i && !B->all_stored_objects_vector.empty())
               sample_Origin = B->all_stored_objects_vector[0]->X;
         });
         auto name = "./output/cube_level" + std::to_string(i) + ".vtp";
         std::ofstream ofs(name);
         std::cout << "name = " << name << std::endl;
         vtkPolygonWrite(ofs, cube_level);
         ofs.close();
      }

      std::cout << "paraview ./output/cube_level*.vtp" << std::endl;

      const auto& target_bucket = B_poles.getBucketAtDeepest(sample_Origin);
      //
      cube_M2L.resize(B_poles.max_level + 1);
      line_M2L.resize(B_poles.max_level + 1);
      poles_in_bucket.resize(B_poles.max_level + 1);
      int l = B_poles.getBucketAtDeepest(sample_Origin)->level;
      std::cout << "l = " << l << std::endl;

      for (int i = 0; i <= l; i++) {
         auto A = B_poles.getBucketAtLevel(i, sample_Origin);
         // std::cout << "A->level = " << A->level << std::endl;
         for (auto& B : A->buckets_for_M2L) {
            line_M2L[B->level].push_back(T2Tddd{B->X, A->X});
            for (auto& pole : B->all_stored_objects_vector)
               poles_in_bucket[B->level].push_back(pole->X);
            T8Tddd t8tddd = (CoordinateBounds)(B->bounds);
            std::cout << "B->level = " << B->level << std::endl;
            cube_M2L[B->level].push_back(t8tddd);
         }
      }

      //! 出力
      for (int i = 0; i <= max_level; i++) {
         {
            std::ofstream ofs("./output/cube_M2L" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, cube_M2L[i]);
            cubePVD.push("cube_level" + std::to_string(i) + ".vtp", i);
            ofs.close();
         }
         {
            std::ofstream ofs("./output/line_M2L" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, line_M2L[i]);
            lineM2LPVD.push("line_M2L" + std::to_string(i) + ".vtp", i);
            ofs.close();
         }
         {
            std::ofstream ofs("./output/poles_in_bucket" + std::to_string(i) + ".vtp");
            vtkPolygonWrite(ofs, poles_in_bucket[i]);
            ofs.close();
         }
      }

      cubePVD.output();
      lineM2LPVD.output();

      std::cout << "paraview ./output/cube_M2L*.vtp" << std::endl;
      std::cout << "paraview ./output/line_M2L*.vtp" << std::endl;

      {
         std::unordered_map<networkPoint*, double> data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13, data14, data15;
#ifdef _BEM_
         std::unordered_map<networkPoint*, double> data1_BEM, data2_BEM, data3_BEM;
#endif
         double max_abs_ign = 0;
         for (const auto& p : obj->getPoints()) {
            auto d_igign = p->igign - p->wGPhin_wGnPhi_FMM;
            auto d_igign_near = p->igign - p->wGPhin_wGnPhi_far;
            auto d_igign_far = p->igign - p->wGPhin_wGnPhi_near;
            data1[p] = d_igign[0];
            data2[p] = d_igign[1];
            data3[p] = p->igign[0];
            data4[p] = p->igign[1];
            data5[p] = std::abs(d_igign[0] / p->igign[0]);
            data6[p] = std::abs(d_igign[1] / p->igign[1]);
            data7[p] = p->wGPhin_wGnPhi_FMM[0];
            data8[p] = p->wGPhin_wGnPhi_FMM[1];
            data9[p] = p->wGPhin_wGnPhi_near[0];
            data10[p] = p->wGPhin_wGnPhi_near[1];
            data11[p] = p->wGPhin_wGnPhi_far[0];
            data12[p] = p->wGPhin_wGnPhi_far[1];
            data14[p] = d_igign_far[0];
            data15[p] = d_igign_far[1];
            data1_BEM[p] = p->CORNER ? 0 : (p->Neumann ? 1 : (p->Dirichlet ? 2 : 3));
            max_abs_ign = std::max(max_abs_ign, std::abs(p->igign[1]));
         }

         for (const auto& p : obj->getPoints()) {
            auto d_igign = p->igign - p->wGPhin_wGnPhi_FMM;
            data13[p] = std::abs(d_igign[1] / max_abs_ign);
         }

         std::vector<std::tuple<std::string, std::unordered_map<networkPoint*, double>>> data = {{"diff_ig", data1},
                                                                                                 {"diff_ign", data2},
                                                                                                 {"ig", data3},
                                                                                                 {"ign", data4},
                                                                                                 {"rel_err_ig", data5},
                                                                                                 {"rel_err_ign", data6},
                                                                                                 {"ig_FMM", data7},
                                                                                                 {"ign_FMM", data8},
                                                                                                 {"ig_near", data9},
                                                                                                 {"ign_near", data10},
                                                                                                 {"ig_far_FMM", data11},
                                                                                                 {"ign_far_FMM", data12},
                                                                                                 {"rel_err_ign_max", data13},
                                                                                                 {"ig_far_true", data14},
                                                                                                 {"ign_far_true", data15},
                                                                                                 {"boundary_type", data1_BEM}};

         std::ofstream ofs("./output/nodes_igign.vtp");
         vtkPolygonWrite(ofs, obj->getPoints(), data);
         ofs.close();
      }

      {
         poles_in_bucket[0].clear();
         for (auto& B : target_bucket->buckets_near)
            for (auto& pole : B->all_stored_objects_vector)
               poles_in_bucket[0].push_back(pole->X);

         std::ofstream ofs("./output/poles_in_bucket_near.vtp");
         vtkPolygonWrite(ofs, poles_in_bucket[0]);
         ofs.close();
      }

      {
         cube_near.clear();
         for (auto& B : target_bucket->buckets_near)
            cube_near.push_back((CoordinateBounds)(B->bounds));
         std::ofstream ofs("./output/cube_near.vtp");
         vtkPolygonWrite(ofs, cube_near);
         ofs.close();
      }
   }
   std::cout << Magenta << "出力" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
}
