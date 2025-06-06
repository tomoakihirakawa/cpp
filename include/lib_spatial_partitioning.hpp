#pragma once
// #define _DEBUG_FMM_

#include <execution>
#include <ranges>
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT lib_spatial_partitioning

## `Bucket`クラス

`Bucket`クラスは，オブジェクトを３次元空間内に配置し，効率的に検索できるようにするための「バケツ（Bucket）」構造を提供します．

WARNING: テンプレート型`T`のオブジェクトは，予め`getX()`を使ってxyz座標を取得できるようにしておく必要がある．

### メンバ変数

| 変数名                           | 説明                                          |
|:-------------------------------:|:--------------------------------------------:|
| `xbounds`, `ybounds`, `zbounds` | X,Y,Z 座標の境界値．`Tdd` 型                    |
| `xsize`, `ysize`, `zsize`       | 各座標のサイズ．`ST` 型                         |
| `bounds`                        | 全体のバウンド．`T3Tdd` 型                      |
| `center`                        | 空間の中心座標．`Tddd` 型                       |
| `dn`                            | 各座標のサイズ（ST3 型）                        |
| `data`                          | バケツ（3D配列）                               |
| `data_vector`                   | バケツのベクター版                              |
| `data_bool`                     | バケツが空であるかどうかを示すbool値の3D配列       |
| `vector_is_set`                 | ベクター版が設定されているか                     |
| `all_stored_objects`            | 保存されている全てのオブジェクト                  |
| `map_to_ijk`                    | オブジェクトからインデックスへのマッピング          |
| `dL`                            | バケツの１辺の長さ                              |

### メソッド

#### 初期化関連

| メソッド | 説明 |
|:--------:|:-------------------------------------------------:|
| `initialize(const T3Tdd &boundingboxIN, const double dL_IN)` | バケツを初期化する． |

#### インデックス変換

| メソッド | 説明 |
|:--------:|:-------------------------------------------------:|
| `itox(const ST i, const ST j, const ST k) const` | インデックス `(i, j, k)` から座標 `(x, y, z)` へ変換． |
| `indices(const Tddd &x) const` | 座標 `x` からインデックス `(i, j, k)` へ変換． |

#### データ追加・削除

| メソッド | 説明 |
|:--------:|:-------------------------------------------------:|
| `add(const Tddd &x, const T p)` | オブジェクト `p` を座標 `x` に追加． |
| `erase(T const p)` | オブジェクト `p` を削除． |

#### データ取得

| メソッド | 説明 |
|:--------:|:-------------------------------------------------:|
| `getData(const std::array<double,3> &x)` | 座標 `x` を含むバケツのデータを取得． |
| `getData(const std::array<double,3> &x, const double d)` | 座標 `x` を中心とする半径 `d` の球状範囲のバケツデータを取得． |

#### 範囲指定検索

| メソッド | 説明 |
|:--------:|:-------------------------------------------------:|
| `apply(const Tddd &x, const double d, const std::function<bool(const T &)> &func)` | バケツの範囲を指定し，その範囲内のオブジェクトに対して `func` を適用． |
| `any_of` | 範囲内のオブジェクトに対して，関数適用結果が1つでも `true` なら `true` を返す． |
| `all_of` | 範囲内のオブジェクトに対して，関数適用結果がすべて `true` なら `true` を返す． |
| `none_of` | 範囲内のオブジェクトに対して，関数適用結果がすべて `false` なら `true` を返す． |

---

### TODO
- `apply` メソッドは，現在**球状の範囲指定**のみ対応している．
  これを拡張し，以下のような**直線上や平面上の範囲指定**もできるようにする．
  - **直線上の範囲指定**: 細かい分割によってインデックス変換する．
  - **平面上の範囲指定**: 平面の方程式を用い，バケツセルとの交差判定を高速に行う．

*/

// Buckets is derived from CoordinateBounds
// N <= 8, 10だと破綻することがあった．
//! ./fast ./input_files/Goring1979_DT0d05_MESHwater0d04refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
template <typename T, int N = 15>
struct Buckets : public CoordinateBounds {

   Moments<N> MomentsMultipoleExpansion;
   Moments<N> MomentsLocalExpansion;

   //$ buckets_for_M2Lは，このバケツが局所展開する対象を保存している．
   std::vector<Buckets<T> *> buckets_for_M2L;
   //$ buckets_for_L2Mは，このバケツが局所展開を受け取る対象を保存している．
   std::vector<Buckets<T> *> buckets_for_L2M;
   std::vector<Buckets<T> *> buckets_near;  //@ 近傍バケツかつ，（重要）子バケツを持たないバケツ．途中までしか育たないツリーを使ったFMMにおいてこの性質が重要．;

   std::vector<std::vector<std::vector<Buckets<T> *>>> buckets;  //! child buckets
   Buckets<T> *parent = nullptr;                                 // 親バケツへのポインタを追加

   std::vector<Buckets<T> *> getAllBucket() {
      std::vector<Buckets<T> *> all_buckets;
      for (auto i = 0; i < this->buckets.size(); ++i)
         for (auto j = 0; j < this->buckets[i].size(); ++j)
            for (auto k = 0; k < this->buckets[i][j].size(); ++k)
               all_buckets.emplace_back(this->buckets[i][j][k]);
      return all_buckets;
   }

   void forEachBuckets(const std::function<void(Buckets<T> *)> &func_for_bucket) {
      for (auto &i : this->buckets)
         for (auto &j : i)
            for (auto &k : j)
               func_for_bucket(k);
   }

   bool has_child = false;  //! このバケツはツリー構造を持っているかどうか
   int level = 0;           //! ツリー全体における，このバケツのレベル
   int max_level = 1;       //! ツリーの最深レベル

   void setLevel(const int levelIN, const int max_levelIN) {
      this->level = levelIN;
      this->max_level = max_levelIN;
      // std::cout << "level : " << Green << this->level << "," << Blue << "max_level : " << this->max_level << colorReset << std::endl;
   }

   //! NOTE
   // !% level_bucketsは，levelがゼロの根っこバケツのみがもつもの．nullptrを含まない．
   // !%すべてのツリーが同じ深さまで成長していない場合，level_bucketsの最後が最も深いバケツとは限らない．
   // !%なのでDeepestは別に実装する必要がある．
   std::vector<std::vector<Buckets<T> *>> level_buckets;
   std::vector<Buckets<T> *> deepest_level_buckets;

   //! この方法で，`data[i][j][k]`を８分割しバケツを作成する．\label{buckets_generateTree}
   bool generateTree(const std::function<bool(const Buckets<T> *)> &condition = [](const Buckets<T> *) { return true; }) {
      this->MomentsMultipoleExpansion.initialize(this->X);
      this->MomentsLocalExpansion.initialize(this->X);
      buckets.clear();
      this->has_child = false;
      if (condition(this)) {
         buckets.resize(this->xsize, std::vector<std::vector<Buckets<T> *>>(this->ysize, std::vector<Buckets<T> *>(this->zsize, nullptr)));
         for (auto [i, j, k] : this->vector_ijk) {
            auto &bucket = buckets[i][j][k];                                          //! Create a local reference to the current bucket
            bucket = new Buckets<T>(getBounds({i, j, k}), this->dL * (0.5 + 1e-13));  //! Initialize the bucket
            bucket->parent = this;                                                    //! 親バケツを設定
            for (auto &p : this->data[i][j][k])
               bucket->add(p->X, p);                             //! Add elements to the bucket if they satisfy the condition
            bucket->setLevel(this->level + 1, this->max_level);  //! Set the level and generate the tree for the bucket
            bucket->generateTreeSingle(condition);
         }
         this->has_child = true;
      }

      // fill level_buckets if the level is 0
      if (this->level == 0) {
         this->level_buckets.clear();
         this->deepest_level_buckets.clear();
         this->level_buckets.resize(this->max_level + 1);  // zero を含むので+1
         this->forEachAll([&](Buckets<T> *b) -> void {
            if (b != nullptr) {
               if (!b->all_stored_objects.empty()) {
                  this->level_buckets[b->level].emplace_back(b);
                  if (!b->has_child)
                     this->deepest_level_buckets.emplace_back(b);
               }
            }
         });
      }

      return this->has_child;
   }

   bool generateTreeSingle(const std::function<bool(const Buckets<T> *)> &condition = [](const Buckets<T> *) { return true; }) {
      this->MomentsMultipoleExpansion.initialize(this->X);
      this->MomentsLocalExpansion.initialize(this->X);
      buckets.clear();
      this->has_child = false;
      if (condition(this)) {
         buckets.resize(this->xsize, std::vector<std::vector<Buckets<T> *>>(this->ysize, std::vector<Buckets<T> *>(this->zsize, nullptr)));
         for (auto [i, j, k] : this->vector_ijk) {
            auto &bucket = buckets[i][j][k];  //! Create a local reference to the current bucket
            // bucket = std::make_shared<Buckets<T>>(getBounds({i, j, k}), this->dL * (0.5 + 1e-13));
            bucket = new Buckets<T>(getBounds({i, j, k}), this->dL * (0.5 + 1e-13));
            bucket->parent = this;  // 親バケツを設定
            for (auto &p : this->data[i][j][k])
               bucket->add(p->X, p);
            bucket->setLevel(this->level + 1, this->max_level);
            bucket->generateTreeSingle(condition);
         }
         this->has_child = true;
      }
      return this->has_child;
   }

   /* =================== Access to the bucket at the level of the tree =================== */

   //@ applyToBucket(1,func)
   //@ bucketsは，level==1のバケツに適用される

   //@ this is traversing using the tree structure.これにlevel_bucketsを使ってはならない．なぜなら，level_bucketsの生成にこの関数が使われるから．
   void forEachAll(const std::function<void(Buckets<T> *)> &func_for_bucket) {
      func_for_bucket(this);
      if (this->has_child) {
         std::for_each(std::execution::unseq, this->buckets.begin(), this->buckets.end(), [&](auto &i) {
            std::for_each(std::execution::unseq, i.begin(), i.end(), [&](auto &j) {
               std::for_each(std::execution::unseq, j.begin(), j.end(), [&](auto &k) {
                  if (k != nullptr) {
                     k->forEachAll(func_for_bucket);
                  }
               });
            });
         });
      }
   }

   void forEachAllParallel(const std::function<void(Buckets<T> *)> &func_for_bucket) {
      func_for_bucket(this);
      if (this->has_child) {
         std::for_each(std::execution::par_unseq, this->buckets.begin(), this->buckets.end(), [&](auto &i) {
            std::for_each(std::execution::par_unseq, i.begin(), i.end(), [&](auto &j) {
               std::for_each(std::execution::par_unseq, j.begin(), j.end(), [&](auto &k) {
                  if (k != nullptr) {
                     k->forEachAllParallel(func_for_bucket);
                  }
               });
            });
         });
      }
   }

   void forEachAtLevel(const int level, const std::function<void(Buckets<T> *)> &func_for_bucket) {
      if (level >= this->level_buckets.size())
         return;

      for (auto &b : this->level_buckets[level])
         func_for_bucket(b);
   }

   void forEachAtLevel(const std::vector<int> &levels, const std::function<void(Buckets<T> *)> &func_for_bucket) {
      for (const auto &lv : levels)
         forEachAtLevel(lv, func_for_bucket);
   }

   void forEachAtLevelParallel(const int lv, const std::function<void(Buckets<T> *)> &func_for_bucket) {
      if (lv >= this->level_buckets.size())
         return;

#pragma omp parallel
      for (auto &b : this->level_buckets[lv])
#pragma omp single nowait
         func_for_bucket(b);
   }

   void forEachAtLevelParallel(const std::vector<int> &levels, const std::function<void(Buckets<T> *)> &func_for_bucket) {
      //! this->level_bucketsは，{0,1,2,3}のように，levelに対応するバケツのリストを持っている．(size()=4)
      //! この場合，level=4を指定するとエラー，これをチェックするには，lv < this->level_buckets.size()とする．
      for (const auto &lv : levels)
         forEachAtLevelParallel(lv, func_for_bucket);
   }

   //    template <typename Func>
   //    void forEachAtDeepestParallel(Func func_for_bucket) {
   // #pragma omp parallel
   //       for (auto &b : this->deepest_level_buckets)
   // #pragma omp single nowait
   //          func_for_bucket(b);
   //    }

   template <typename Func>
   void forEachAtDeepestParallel(const Func &func_for_bucket) {
      // for (auto &b : this->deepest_level_buckets)
      //    func_for_bucket(b);
      std::for_each(std::execution::par_unseq, this->deepest_level_buckets.begin(), this->deepest_level_buckets.end(), func_for_bucket);
   }

   template <typename Func>
   void forEachAtDeepest(const Func &func_for_bucket) {
      for (auto &b : this->deepest_level_buckets)
         func_for_bucket(b);
   }

   Buckets<T> *getBucketAtLevel(const int level, const Tddd &x) const {
      if (!this->has_child || level > this->max_level)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "The tree structure does not exist or the level is too large.");

      auto ijk = this->indices(x);
      auto ret = this->buckets[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
      int count = 0;
      while (ret->level < level) {
         ijk = ret->indices(x);
         ret = ret->buckets[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
         if (count++ > 10)
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "The level is too large.");
      }
      return ret;
   };

   Buckets<T> *getBucketAtDeepest(const Tddd &x) const {
      if (!this->has_child)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "The tree structure does not exist.");
      auto ret = this->getBucket(x);
      while (ret->has_child)
         ret = ret->getBucket(x);
      return ret;
   };

   /* ===================================================================================== */

   using sizeType = int;
   using ST = sizeType;
   using ST2 = std::array<sizeType, 2>;
   using ST3 = std::array<sizeType, 3>;
   using ST6 = std::array<sizeType, 6>;
   //! getX()でxyz座標を取得できるオブジェクトTのための，バケツ
   // Tdd xbounds, ybounds, zbounds;
   // T3Tdd bounds;
   ST xsize, ysize, zsize;
   Tddd center;
   ST3 dn;
   std::vector<std::vector<std::vector<std::unordered_set<T>>>> data;
   std::vector<std::vector<std::vector<std::vector<T>>>> data_vector;
   std::vector<std::vector<std::vector<bool>>> data_bool;
   bool vector_is_set;
   std::unordered_set<T> all_stored_objects;
   std::vector<T> all_stored_objects_vector;
   std::unordered_map<T, ST3> map_to_ijk;
   double dL;
   double bucketVolume() const { return std::pow(this->dL, 3.); };
   //
   Buckets() = default;
   Buckets(const CoordinateBounds &c_bounds) : CoordinateBounds(c_bounds) {
      auto [xmin, xmax] = c_bounds.bounds[0];
      auto [ymin, ymax] = c_bounds.bounds[1];
      auto [zmin, zmax] = c_bounds.bounds[2];
      double max_dL = std::max(std::abs(xmax - xmin), std::max(std::abs(ymax - ymin), std::abs(zmax - zmin)));
      initialize(this->bounds, max_dL);
      this->center = this->X;
   };
   Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : CoordinateBounds(c_bounds) {
      initialize(this->bounds, dL_IN);
      this->center = this->X;
   };
   Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : CoordinateBounds(boundingboxIN) {
      initialize(this->bounds, dL_IN);
      this->center = this->X;
   };

   ~Buckets() {
      // std::cout << "Deleting Buckets at level: " << this->level << std::endl;
      for (auto &i : this->buckets)
         for (auto &j : i)
            for (auto &k : j)
               if (k != nullptr)
                  delete k;
      // std::cout << "Buckets at level " << this->level << " deleted." << std::endl;
   }

   std::vector<std::array<int, 3>> vector_ijk;
   void initialize(const auto &boundingboxIN, const double dL_IN) {
      this->all_stored_objects.clear();
      this->all_stored_objects_vector.clear();
      CoordinateBounds::setBounds(boundingboxIN);
      this->dL = dL_IN;

      const std::size_t max_size = 1000000;

      auto calculate_size = [&](const auto &bounds) {
         double size = std::ceil((std::get<1>(bounds) - std::get<0>(bounds)) / this->dL);
         if (size <= 0.) size = 1.;
         if (size > max_size) {
            std::cerr << "The size of the bucket is too large (" << size << "). Please reduce the bucket size." << std::endl;
            exit(1);
         }
         return static_cast<ST>(size);
      };

      this->xsize = calculate_size(this->xbounds());
      this->ysize = calculate_size(this->ybounds());
      this->zsize = calculate_size(this->zbounds());

      this->dn = {xsize, ysize, zsize};

      this->data.assign(xsize, std::vector<std::vector<std::unordered_set<T>>>(ysize, std::vector<std::unordered_set<T>>(zsize)));
      this->data_bool.assign(xsize, std::vector<std::vector<bool>>(ysize, std::vector<bool>(zsize, false)));
      this->vector_is_set = false;

      vector_ijk.clear();
      vector_ijk.reserve(this->xsize * this->ysize * this->zsize);
      for (auto i = 0; i < this->xsize; ++i)
         for (auto j = 0; j < this->ysize; ++j)
            for (auto k = 0; k < this->zsize; ++k)
               vector_ijk.emplace_back(std::array<int, 3>{i, j, k});

      this->buckets.reserve(this->xsize);
   }

   const std::vector<T> &getAll() const { return this->all_stored_objects_vector; };
   //@ -------------------------------- インデックス変換 -------------------------------- */
   // x座標を内包するバケツのインデックスを返す
   constexpr Tddd itox(const ST i, const ST j, const ST k) const {
      return {this->dL * 0.5 + this->dL * i + std::get<0>(this->xbounds()),
              this->dL * 0.5 + this->dL * j + std::get<0>(this->ybounds()),
              this->dL * 0.5 + this->dL * k + std::get<0>(this->zbounds())};
   };
   constexpr Tddd itox(const ST3 &ijk) const { return itox(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk)); };
   constexpr ST3 indices_no_clamp(const Tddd &x) const {
      //! floorは必要！もし直接intキャストを使うと，-0.**が0になってしまい，isInsideがfalseなはずがtrueが返ってしまう．
      return {static_cast<ST>(std::floor((std::get<0>(x) - std::get<0>(this->xbounds())) / this->dL)),
              static_cast<ST>(std::floor((std::get<1>(x) - std::get<0>(this->ybounds())) / this->dL)),
              static_cast<ST>(std::floor((std::get<2>(x) - std::get<0>(this->zbounds())) / this->dL))};
   };
   constexpr ST3 indices(const Tddd &x) const {
      /*
      //! intキャストはゼロ方向へ実数を切り捨てた結果を返すので，static_cast<int>によって正しくセルのインデックスに変換できる．
       0.**   1.**   2.**   3.**
      <-dL-> <-dL-> <-dL-> <-dL->
      *-----*------*------*------*
      |  0  |   1  |   2  |   3  |
      *-----*------*------*------*
      */
      const auto ijk = indices_no_clamp(x);
      return {std::clamp(std::get<0>(ijk), static_cast<ST>(0), this->xsize - 1),
              std::clamp(std::get<1>(ijk), static_cast<ST>(0), this->ysize - 1),
              std::clamp(std::get<2>(ijk), static_cast<ST>(0), this->zsize - 1)};
   };
   constexpr ST6 indices_ranges(const Tddd &x, const double d) const {
      auto [i_min, j_min, k_min] = indices(x - d);
      auto [i_max, j_max, k_max] = indices(x + d);
      return {i_min, i_max, j_min, j_max, k_min, k_max};
   };

   constexpr ST6 indices_ranges(const Tddd &x, const Tddd &dx_dy_dz) const {
      auto [dx, dy, dz] = dx_dy_dz;
      auto minmax_x = indices_ranges(x, dx);
      auto minmax_y = indices_ranges(x, dy);
      auto minmax_z = indices_ranges(x, dz);
      return {std::get<0>(minmax_x), std::get<1>(minmax_x),
              std::get<2>(minmax_y), std::get<3>(minmax_y),
              std::get<4>(minmax_z), std::get<5>(minmax_z)};
   };
   //@ ------------------------ インデックスがboundsに収まっているかどうか ------------------------ */
   bool InsideQ(const ST i, const ST j, const ST k) const { return (i >= 0 && j >= 0 && k >= 0 && i < this->xsize && j < this->ysize && k < this->zsize); };
   bool InsideQ(const ST3 &ijk) const { return InsideQ(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk)); };
   bool InsideQ(const Tddd &x) const { return InsideQ(indices_no_clamp(x)); };
   void erase(T const p) {
      auto it = this->map_to_ijk.find(p);
      if (it != this->map_to_ijk.end()) {
         auto [i, j, k] = it->second;
         this->all_stored_objects.erase(p);
         auto p_it = std::find(this->all_stored_objects_vector.begin(), this->all_stored_objects_vector.end(), p);
         if (p_it != this->all_stored_objects_vector.end())
            this->all_stored_objects_vector.erase(p_it);
         this->data[i][j][k].erase(p);
         this->map_to_ijk.erase(it);
      }
   };
   //@ -------------------------------------------------------------------------- */
   auto getBounds(const ST3 &ijk) const {
      const auto [i, j, k] = ijk;
      return CoordinateBounds(T3Tdd{{{this->dL * i + std::get<0>(this->xbounds()), this->dL * (i + 1) + std::get<0>(this->xbounds())},
                                     {this->dL * j + std::get<0>(this->ybounds()), this->dL * (j + 1) + std::get<0>(this->ybounds())},
                                     {this->dL * k + std::get<0>(this->zbounds()), this->dL * (k + 1) + std::get<0>(this->zbounds())}}});
   };
   //@ -------------------------------------------------------------------------- */
   // unordered_setをvectorにコピー
   void setVector() {
      this->data_vector.clear();
      this->data_vector.resize(this->xsize, std::vector<std::vector<std::vector<T>>>(this->ysize, std::vector<std::vector<T>>(this->zsize)));

      // Iterate through the 3D array and populate the 3D vector
      ST i, j, k;
      for (i = 0; i < this->xsize; ++i) {
         for (j = 0; j < this->ysize; ++j) {
            for (k = 0; k < this->zsize; ++k) {
               auto &bucket = this->data[i][j][k];
               auto &target_vector = this->data_vector[i][j][k];
               target_vector.reserve(bucket.size());
               target_vector.assign(bucket.begin(), bucket.end());
            }
         }
      }

      std::cout << "      bounds : " << this->bounds << std::endl;
      std::cout << "          dL : " << this->dL << std::endl;
      std::cout << "bucket sizes : {" << this->xsize << ", " << this->ysize << ", " << this->zsize << "}" << std::endl;
      this->vector_is_set = true;
   };
   //@ -------------------------------------------------------------------------- */
   //! explicitly specify coordinates
   bool add(const ST3 &ijk, const T p) {
      this->vector_is_set = false;
      const auto [i, j, k] = ijk;
      this->map_to_ijk[p] = ijk;
      const bool bucket_inserted = this->data[i][j][k].emplace(p).second;
      const bool all_objects_inserted = this->all_stored_objects.emplace(p).second;
      if (all_objects_inserted)
         this->all_stored_objects_vector.emplace_back(p);
      return bucket_inserted && all_objects_inserted;
   };
   bool add(const Tddd &x, const T p) {
      this->vector_is_set = false;
      return add(indices(x), p);
   };

   std::tuple<bool, std::array<int, 3>> addAndGetIndices(const Tddd &x, const T p) {
      this->vector_is_set = false;
      auto ijk = indices(x);
      const auto [i, j, k] = ijk;
      this->map_to_ijk[p] = ijk;
      const bool bucket_inserted = this->data[i][j][k].emplace(p).second;
      const bool all_objects_inserted = this->all_stored_objects.emplace(p).second;
      if (all_objects_inserted)
         this->all_stored_objects_vector.emplace_back(p);
      return {bucket_inserted && all_objects_inserted, ijk};
   };

   //! automatically specify coordinates
   bool add(const std::unordered_set<T> &P) {
      this->vector_is_set = false;
      bool ret = true;
      for (const auto p : P) {
         if (!add(indices(ToX(p)), p))
            ret = false;
      }
      return ret;
   };

   bool add(const std::unordered_set<T> &P, const std::function<Tddd(const T &)> &ToX) {
      this->vector_is_set = false;
      bool ret = true;
      for (const auto p : P) {
         if (!add(indices(ToX(p)), p))
            ret = false;
      }
      return ret;
   };

   /* -------------------------------------------------------------------------- */
   auto getBucket(const Tddd &x) const {
      const auto ijk = this->indices(x);
      return this->buckets[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
   };

   auto getData(const Tddd &x) const {
      const auto ijk = this->indices(x);
      return this->data[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
   };

   auto getData(const Tddd &x, const double &range) const {
      const auto ijk = this->indices(x);
      std::unordered_set<T> ret;
      this->apply(x, range, [&](auto &a) { ret.emplace(a); });
      return ret;
   };

   void clear() {
      this->vector_is_set = false;
      this->data.clear();
      this->all_stored_objects.clear();
      this->all_stored_objects_vector.clear();
   };
   // b@ -------------------------------------------------------------------------- */
   // b@                             STL like functions                             */
   // b@ -------------------------------------------------------------------------- */

   //! none_of
   bool none_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (this->data.empty()) {
         return true;
      }
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
      if (!this->vector_is_set) {
         return !std::any_of(std::execution::unseq, this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::any_of(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::any_of(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  // return std::any_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), func);
                  return std::any_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); });
               });
            });
         });
      } else {
         return !std::any_of(std::execution::unseq, this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::any_of(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::any_of(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  // return std::any_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), func);
                  return std::any_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); });
               });
            });
         });
      }
   }

   //! all_of
   bool all_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (this->data.empty()) {
         return true;
      }
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);

      if (!this->vector_is_set) {
         return std::all_of(std::execution::unseq, this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::all_of(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::all_of(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  // return std::all_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), func);
                  return std::all_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); });
               });
            });
         });
      } else {
         return std::all_of(std::execution::unseq, this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::all_of(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::all_of(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  // return std::all_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), func);
                  return std::all_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); });
               });
            });
         });
      }
   }

   //! any_of
   bool any_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (this->data.empty()) {
         return false;
      }
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);

      if (!this->vector_is_set) {
         return std::any_of(std::execution::unseq, this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::any_of(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::any_of(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  // return std::any_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), func);
                  return std::any_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); });
               });
            });
         });
      } else {
         return std::any_of(std::execution::unseq, this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::any_of(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::any_of(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  // return std::any_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), func);
                  return std::any_of(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { return func(p); });
               });
            });
         });
      }
   }

   //! apply
   void apply(const Tddd &x, const double d, const std::function<void(const T &)> &func) const {
      // if (this->data.empty())
      //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");

      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);

      if (!this->vector_is_set) {
         std::for_each(std::execution::unseq, this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
            std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
               std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
                  // for (const auto &p : Bijk) func(p);
                  std::for_each(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { func(p); });
               });
            });
         });
      } else {
         std::for_each(std::execution::unseq, this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
            std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
               std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
                  // for (const auto &p : Bijk) func(p);
                  std::for_each(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { func(p); });
               });
            });
         });
      }
   }
   //
   //! apply
   void applyVec1D(const Tddd &x, const double d, const std::function<void(const std::vector<T> &)> &func) const {
      if (this->data.empty() || !this->vector_is_set)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
      std::for_each(std::execution::unseq, this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
         std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
            std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
               func(Bijk);
            });
         });
      });
   }

   //! apply
   void apply(const Tddd &x, const double d, const std::function<void(const int, const int, const int)> &func) const {
      // if (this->data.empty())
      //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
      int i, j, k;
      for (i = i_min; i <= i_max; ++i)
         for (j = j_min; j <= j_max; ++j)
            for (k = k_min; k <= k_max; ++k)
               func(i, j, k);
   }

   //! apply
   void apply_to_the_nearest_bound(const Tddd &x, const std::function<void(const int, const int, const int)> &func) const {
      if (this->data.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");
      auto [I, J, K] = indices(x);

      int values[] = {I, J, K, this->xsize - I, this->ysize - J, this->zsize - K};
      int smallest = *std::min_element(values, values + 6);
      int smallestIndex = std::distance(values, std::find(values, values + 6, smallest));

      int i, j, k;

      if (smallestIndex == 0) {
         for (i = I; i >= 0; --i)
            func(i, J, K);
      } else if (smallestIndex == 1) {
         for (j = J; j >= 0; --j)
            func(I, j, K);
      } else if (smallestIndex == 2) {
         for (k = K; k >= 0; --k)
            func(I, J, k);
      } else if (smallestIndex == 3) {
         for (i = I; i < this->xsize; ++i)
            func(i, J, K);
      } else if (smallestIndex == 4) {
         for (j = J; j < this->ysize; ++j)
            func(I, j, K);
      } else if (smallestIndex == 5) {
         for (k = K; k < this->zsize; ++k)
            func(I, J, k);
      }
   }

   //! apply
   void apply(const Tddd &x, const Tddd dx_dy_dz, const std::function<void(const T &)> &func) const {
      // if (this->data.empty())
      //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");

      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, dx_dy_dz);

      if (!this->vector_is_set) {
         std::for_each(std::execution::unseq, this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
            std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
               std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
                  // for (const auto &p : Bijk) func(p);
                  std::for_each(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { func(p); });
               });
            });
         });
      } else {
         std::for_each(std::execution::unseq, this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
            std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
               std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
                  // for (const auto &p : Bijk) func(p);
                  std::for_each(std::execution::unseq, Bijk.cbegin(), Bijk.cend(), [&](const auto p) { func(p); });
               });
            });
         });
      }
   };

   /* ----------------------------- FOR LINE SEARCH ---------------------------- */

   std::vector<ST3> line2indices(const Tddd &A, const Tddd &B) const {
      std::unordered_set<ST3> uniqueIndexSet;
      for (const auto &X : Subdivide(A, B, std::ceil(Norm(A - B) / (this->dL / 2.))))
         this->apply(X, this->dL, [&uniqueIndexSet](const int i, const int j, const int k) { uniqueIndexSet.insert({i, j, k}); });
      std::vector<ST3> uniqueIndices(uniqueIndexSet.begin(), uniqueIndexSet.end());
      return uniqueIndices;
   }

   //! apply using indices
   void apply(const std::vector<ST3> &V_ijk, const std::function<void(const T &)> &func) const {
      if (!this->vector_is_set) {
         for (const auto &ijk : V_ijk)
            for (const auto &p : this->data[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)])
               func(p);
      } else {
         for (const auto &ijk : V_ijk)
            for (const auto &p : this->data_vector[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)])
               func(p);
      }
   }
};

template <typename T>
Buckets<T> copyPartition(const auto &buckets) {
   Buckets<T> ret(buckets.bounds, buckets.dL);
   ret.buckets.resize(buckets.xsize, std::vector<std::vector<Buckets<T> *>>(buckets.ysize, std::vector<Buckets<T> *>(buckets.zsize, nullptr)));
   ret.level = buckets.level;
   ret.max_level = buckets.max_level;
   ret.has_child = buckets.has_child;
   ret.xsize = buckets.xsize;
   ret.ysize = buckets.ysize;
   ret.zsize = buckets.zsize;
   ret.center = buckets.center;
   ret.dn = buckets.dn;
   ret.data.resize(buckets.xsize, std::vector<std::vector<std::unordered_set<T>>>(buckets.ysize, std::vector<std::unordered_set<T>>(buckets.zsize, std::unordered_set<T>{})));
   ret.data_vector.resize(buckets.xsize, std::vector<std::vector<std::vector<T>>>(buckets.ysize, std::vector<std::vector<T>>(buckets.zsize)));
   ret.data_bool.resize(buckets.xsize, std::vector<std::vector<bool>>(buckets.ysize, std::vector<bool>(buckets.zsize, false)));
   ret.vector_is_set = false;
   return ret;
}

// template <typename T>
// struct Buckets : public BaseBuckets<T> {
//    Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<T>(c_bounds, dL_IN){};
//    Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<T>(boundingboxIN, dL_IN){};
// };

/* -------------------------------------------------------------------------- */

const double scale = 3;
bool isFar(auto *A, auto *B) { return !InsideQ(B->X, A->scaledBounds(scale)); };
bool isNear(auto *A, auto *B) { return InsideQ(B->X, A->scaledBounds(scale)); };
// bool isFar(const std::shared_ptr<auto> &A, const std::shared_ptr<auto> &B) { return !InsideQ(B->X, A->scaledBounds(scale)); };
// bool isNear(const std::shared_ptr<auto> &A, const std::shared_ptr<auto> &B) { return InsideQ(B->X, A->scaledBounds(scale)); };
bool isFar(const auto *A, const auto *B);

template <typename T>
void checkAndAddBuckets(T A, T B) {
   if (isNear(A, B)) {

      if (!A->all_stored_objects_vector.empty() && !B->all_stored_objects_vector.empty())
         if (A != B && isNear(A, B)) {
            A->buckets_near.emplace_back(B);
            // isnear = true;
         }

      for (auto &A_c : A->getAllBucket()) {
         if (!A_c->all_stored_objects_vector.empty() || A_c->has_child) {
            for (auto &B_c : B->getAllBucket()) {
               if (!B_c->all_stored_objects_vector.empty()) {
                  // bool isnear = false;

                  if (A_c != B_c && isFar(A_c, B_c)) {
                     A_c->buckets_for_M2L.emplace_back(B_c);
                  } else {
                     checkAndAddBuckets(A_c, B_c);
                  }
               }
            }
         }
      }
   }
}

//% -------------------------------------------------------------------------- */
//%        Functions Related to Multipole Expansion Utilizing Buckets          */
//% -------------------------------------------------------------------------- */

void setM2M(auto &B_poles) {
   TimeWatch tw, tw_all;
   for (int level = B_poles.max_level - 1; level >= 0; level--) {
      B_poles.forEachAtLevel({level}, [&](auto *B) {
         B->MomentsMultipoleExpansion.set_m2m(B->getAllBucket());
      });
#if defined(_DEBUG_FMM_)
      std::cout << yellow << "setM2M" << ", level=" << level << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
#endif
   }
   std::cout << Yellow << "setM2M" << Green << ", Elapsed time : " << tw_all() << colorReset << std::endl;
}

void M2M(auto &B_poles) {
   TimeWatch tw, tw_all;
   for (int level = B_poles.max_level - 1; level >= 0; level--) {
      B_poles.forEachAtLevel({level}, [&](auto *B) {
         B->MomentsMultipoleExpansion.m2m();
      });
#if defined(_DEBUG_FMM_)
      std::cout << magenta << "M2M" << ", level=" << level << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
#endif
   }
   // std::cout << Magenta << "M2M" << Green << ", Elapsed time : " << tw_all() << colorReset << std::endl;
};

/* -------------------------------------------------------------------------- */

void setL2L(auto &B_poles) {
   TimeWatch tw, tw_all;
   int level = 0;
   for (auto &buckets_from_top_level : B_poles.level_buckets) {
      for (auto &parent : buckets_from_top_level) {
         parent->forEachBuckets([&](auto *child) {
            child->MomentsLocalExpansion.set_l2l(&parent->MomentsLocalExpansion);
         });
      }
#if defined(_DEBUG_FMM_)
      std::cout << yellow << "setL2L" << ", level=" << level << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
#endif
      level++;
   }
   std::cout << Yellow << "setL2L" << Green << ", Elapsed time : " << tw_all() << colorReset << std::endl;
}

void L2L(auto &B_poles) {
   TimeWatch tw, tw_all;
   int level = 0;
   for (auto &buckets_from_top_level : B_poles.level_buckets) {
      for (auto &parent : buckets_from_top_level) {
         parent->forEachBuckets([&](auto *child) {
            child->MomentsLocalExpansion.l2l();
         });
      }
#if defined(_DEBUG_FMM_)
      std::cout << magenta << "L2L" << ", level=" << level << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
#endif
      level++;
   }
   // std::cout << Magenta << "L2L" << Green << ", overall Elapsed time : " << tw_all() << colorReset << std::endl;
}

/* -------------------------------------------------------------------------- */

// ツリーはとりあえず生成しなければならない．
// 極が少ないバケツは，M2Lの時に省略する．

void setBucketsForM2L(auto &B_poles) {
   B_poles.forEachAll([&](auto *A) {
      A->buckets_for_M2L.clear();
      A->buckets_near.clear();
      A->buckets_for_L2M.clear();
   });
   //! store M2L buckets for the bucket at level 1
   //! ここは，Aの展開係数を，M2Lすべきバケツに保存する．Aが空なら，M2LすべきバケツはAにとってないことになる．
   B_poles.forEachAtLevel({1}, [&](auto *A) {
      //! (1) この設定では，Mする場所が，ソース点があるバケツに限られる．これは，常に妥当な設定である．
      if (!A->all_stored_objects_vector.empty())
         B_poles.forEachAtLevel({1}, [&](auto *B) {  //$ check if another bucket is inside the bucket. If it is not, add it to the list of buckets for M2L
            if (!B->all_stored_objects_vector.empty()) {
               //! (2) この設定は，Lできる場所が，節点などに限られる．しかも，nearにすら含めない．．．．
               /*
                  +----+----+----+----+----+
                  |    |    |    |    |    |
                  +----+----+----+----+----+
                  |    | X  | X  | X  |    |
                  +----+----+----+----+----+
                  |    | X  | A  | X  |    |
                  +----+----+----+----+----+
                  |    | X  | X  | X  |    |
                  +----+----+----+----+----+
                  |    |    |    |    |    |
                  +----+----+----+----+----+
               */
               if (isFar(A, B))
                  A->buckets_for_M2L.emplace_back(B);
               else
                  checkAndAddBuckets(A, B);  // % 次のレベルに移動する．
            }
         });
   });

   for (auto &buckets_from_top_level : B_poles.level_buckets) {
      for (auto &A : buckets_from_top_level)
         for (auto &B : A->buckets_for_M2L)
            B->buckets_for_L2M.emplace_back(A);
   }
}

void setM2L(auto &B_poles) {
   TimeWatch tw, tw_all;
   // std::cout << "各レベルの各セルのM2Lの相手を保存する" << std::endl;
   setBucketsForM2L(B_poles);
   // A -> M2L -> B
   int level = 0;
   for (auto &buckets_at_a_level : B_poles.level_buckets) {
#pragma omp parallel
      for (auto &A : buckets_at_a_level)
#pragma omp single nowait
         A->MomentsLocalExpansion.set_m2l(A->buckets_for_L2M);
      std::cout << yellow << "setM2L" << ", level=" << level << ", buckets_at_a_level.size()=" << buckets_at_a_level.size() << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
      level++;
   }
   std::cout << Yellow << "setM2L" << Green << ", Elapsed time : " << tw_all() << colorReset << std::endl;
};

void M2L(auto &B_poles) {
   for (auto &buckets_at_a_level : B_poles.level_buckets) {
#pragma omp parallel
      for (auto &A : buckets_at_a_level)
#pragma omp single nowait
         A->MomentsLocalExpansion.m2l();
   }
};

/* -------------------------------------------------------------------------- */

/*
source4FMMは，
Tddd X;
Tdd weights;
Tddd normal;
std::function<Tdd()> getValues;
を持つ．これらを使うと，直接数値面積分ができる．
*/

/*

BIE：alpha_phi = w_Gn_phin - w_Gn_phiの，w_Gn_phin，w_Gn_phiをまず計算する．

ここで，w_Gn_phi，w_Gn_phiはそれぞれ，

```math
\begin{align*}
\int\int_{\Gamma} G \nabla \phi \cdot {\boldsymbol n} dS\\
\int\int_{\Gamma} \phi \nabla G \cdot {\boldsymbol n} dS
\end{align*}
```

w_Gn_phiのincrementにおいて，符号が負になっている．
これは，勘違いして修正してしまいそうだが，w_Gn_phi自体が負になっているためで，
BIEの符号と混同しないように注意する．

BIEは次の形である．

alpha_phi = w_Gn_phin - w_Gn_phi

多重極展開（see L2P）も，w_Gn_phin，w_Gn_phiの近似を計算しており，直接積分しているものと対応している．

*/

/* -------------------------------------------------------------------------- */

void initializeFMM(auto &B_poles, const auto &targets) {

   for (auto pole : B_poles.all_stored_objects_vector)
      pole->update();

   B_poles.forEachAll([&](auto *b) {
      b->MomentsMultipoleExpansion.initialize();  //! MとM_の初期化
      b->MomentsLocalExpansion.initialize();      //! MとM_の初期化
   });

   B_poles.forEachAtDeepestParallel([&](auto *b) {
      b->MomentsMultipoleExpansion.increment_M(b->all_stored_objects_vector);  // !vector_weightedSourceDensities_RRn_MMの初期化
   });

   setM2M(B_poles);
   setM2L(B_poles);

   //! check number of couping
   int level = 0;
   for (auto &buckets_at_a_level : B_poles.level_buckets) {
      std::cout << "level=" << level++ << std::endl;
      if (buckets_at_a_level.size() == 0) continue;
      double size = 0.;
      for (auto &A : buckets_at_a_level)
         size += A->buckets_for_M2L.size();
      std::cout << "mean A->buckets_for_M2L=" << size / buckets_at_a_level.size() << std::endl;
   }

   setL2L(B_poles);

   TimeWatch tw;
   std::cout << "各レベルの各セルのL2Pの相手を保存する" << std::endl;
#pragma omp parallel
   for (auto &t : targets)
#pragma omp single nowait
      t->setL2P(B_poles);  // !L2Pの初期化
   std::cout << Yellow << "setL2P" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   std::cout << "各レベルの各セルのDirect integrationの相手を保存する" << std::endl;
#pragma omp parallel
   for (auto &t : targets)
#pragma omp single nowait
      t->setDirectIntegration(B_poles);  // !Direct integrationの初期化
   std::cout << Yellow << "setDirectIntegration" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   double size = 0.;
   for (auto &t : targets) {
      if (t->vec_phiphin_WGNWGnN.size() == 0)
         continue;
      size = size + t->vec_phiphin_WGNWGnN.size();
   }
   std::cout << "mean vec_phiphin_WGNWGnN=" << size / targets.size() << std::endl;
};

/* -------------------------------------------------------------------------- */

void updateFMM(auto &B_poles) {

   // TimeWatch tw;

   for (auto pole : B_poles.all_stored_objects_vector)
      pole->update();

   // std::cout << Yellow << "update poles" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   // !この操作は，Momentsのvector_weightedSourceDensities_RRn_MM_を削除しない．
   // !vector_weightedSourceDensities_RRn_MM_の保存するポインタ変数を初期化するだけで，配列内のポインタは生きたまま残る
   // !increment_Mを使うと，vector_weightedSourceDensities_RRn_MM_は初期化され再計算される
   B_poles.forEachAllParallel([&](auto *b) {
      b->MomentsMultipoleExpansion.initialize();
      b->MomentsLocalExpansion.initialize();
   });

   // std::cout << Yellow << "reset Moments" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   B_poles.forEachAtDeepestParallel([&](auto *b) {
      b->MomentsMultipoleExpansion.increment_M_reuse();
   });

   // std::cout << Yellow << "increment_M_reuse" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;

   M2M(B_poles);
#if defined(_DEBUG_FMM_)
   std::cout << Magenta << "M2L ..." << colorReset << std::endl;
#endif
   M2L(B_poles);
#if defined(_DEBUG_FMM_)
   std::cout << Magenta << "L2L ..." << colorReset << std::endl;
#endif
   L2L(B_poles);

   // std::cout << Yellow << "updateFMM" << Green << ", Elapsed time : " << tw() << colorReset << std::endl;
};

/* -------------------------------------------------------------------------- */