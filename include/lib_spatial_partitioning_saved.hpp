#ifndef lib_spatial_partitioning_H
#define lib_spatial_partitioning_H

#include <execution>
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT lib_spatial_partitioning

## `Bucket`クラス

`Bucket`クラスは，オブジェクトを３次元空間内に配置し，効率的に検索できるようにするための「バケツ（Bucket）」構造を提供します．

WARNING: テンプレート型`T`のオブジェクトは，予め`getX()`を使ってxyz座標を取得できるようにしておく必要がある．

### メンバ変数

| 変数名                | 説明                                              |
|:---------------------:|:-------------------------------------------------:|
| `xbounds`, `ybounds`, `zbounds` | X,Y,Z 座標の境界値．`Tdd` 型                   |
| `xsize`, `ysize`, `zsize` | 各座標のサイズ．`ST` 型                          |
| `bounds`              | 全体のバウンド．`T3Tdd` 型                        |
| `center`              | 空間の中心座標．`Tddd` 型                          |
| `dn`                  | 各座標のサイズ（ST3 型）                           |
| `data`             | バケツ（3D配列）                                  |
| `data_vector`      | バケツのベクター版                                 |
| `data_bool`        | バケツが空であるかどうかを示すbool値の3D配列       |
| `vector_is_set`       | ベクター版が設定されているか                       |
| `all_stored_objects`  | 保存されている全てのオブジェクト                   |
| `map_to_ijk`          | オブジェクトからインデックスへのマッピング         |
| `dL`                  | バケツの１辺の長さ                                 |

### メソッド

#### 初期化関連

- `initialize(const T3Tdd &boundingboxIN, const double dL_IN)`: バケツを初期化する．
œa
#### インデックス変換¸

- `itox(const ST i, const ST j, const ST k) const`: インデックスから座標へ変換．
- `indices(const Tddd &x) const`: 座標からインデックスへ変換．

#### データ追加・削除

- `add(const Tddd &x, const T p)`: オブジェクトを追加．
- `erase(T const p)`: オブジェクトを削除．

#### その他

`apply(const Tddd &x, const double d, const std::function<bool(const T &)> &func)`は，バケツの範囲を指定して，その範囲内のオブジェクトに対して関数を適用する．
これと似た関数として，

* `any_of`
* `all_of`
* `none_of`

があり，それぞれ，バケツの範囲内のオブジェクトに対して，関数を適用し，その結果が，それぞれ，`true`，`false`，`false`であれば，`true`を返す．

TODO: これらの関数は，`apply`はある点を中心として半径`d`の球状の範囲を指定することができる．これは球状の範囲を指定していることになる．このような範囲指定以外に，直線上の範囲指定や，平面上の範囲指定などもできるようにしたい．
そのためには，バケツのセルと，線分や平面の交差判定を高速に行う関数が必要になる．
ラフに行っても問題ない．
線に関しては細かい分割によってインデックス変換できる．
平面に関しては，平面の方程式を使って，バケツのセルとの交差判定を行う．

*/

// Buckets is derived from CoordinateBounds

template <typename T>
struct Buckets : public CoordinateBounds {

   /* =============================== Tree structure ================================= */

   ExpCoeffs<10> multipole_expansion;
   ExpCoeffs<10> local_expansion;

   std::vector<std::vector<std::vector<std::shared_ptr<Buckets<T>>>>> buckets;

   std::vector<std::shared_ptr<Buckets<T>>> getAllBucket() {
      std::vector<std::shared_ptr<Buckets<T>>> all_buckets;
      for (auto i = 0; i < this->buckets.size(); ++i)
         for (auto j = 0; j < this->buckets[i].size(); ++j)
            for (auto k = 0; k < this->buckets[i][j].size(); ++k)
               all_buckets.emplace_back(this->buckets[i][j][k]);
      return all_buckets;
   }

   bool has_tree = false;  //! このバケツはツリー構造を持っているかどうか
   int level = 0;          //! ツリー全体における，このバケツのレベル
   int max_level = 1;      //! ツリーの最深レベル

   void setLevel(const int level, const int max_level) {
      this->level = level;
      this->max_level = max_level;
   }

   void generateTree(const std::function<bool(const T &)> &condition = [](const T &) { return true; }) {
      if (this->level >= this->max_level)
         return;

      buckets.resize(this->xsize, std::vector<std::vector<std::shared_ptr<Buckets<T>>>>(this->ysize, std::vector<std::shared_ptr<Buckets<T>>>(this->zsize, nullptr)));

      for (auto i = 0; i < this->data.size(); ++i)
         for (auto j = 0; j < this->data[i].size(); ++j)
            for (auto k = 0; k < this->data[i][j].size(); ++k) {

               //! この方法で，`data[i][j][k]`を８分割しバケツを作成する．\label{buckets_generateTree}
               /*
               1. data[i][j][k]のcoordinate boundsを取得し，その領域を８分割したバケツとしてbuckets[i][j][k]を初期化する．
               2. data[i][j][k]に含まれるオブジェクトのうち，conditionを満たすものだけをbuckets[i][j][k]にコピーする．自動で８つのバケツに分配される．
               3. 再帰的にgenerateTreeを呼び出し，バケツのツリーを生成する．

               ```cpp
               if (this->level >= this->max_level) {
                  this->has_tree = false;
                  return;
               } else
                  this->has_tree = true;
               ```

               によって，バケツの最大レベルに達したら，再帰を終了する．
               */

               auto bounds = getBounds({i, j, k});
               buckets[i][j][k] = std::make_shared<Buckets<T>>(bounds, this->dL * 0.5 + 1e-10);
               buckets[i][j][k]->setLevel(this->level + 1, this->max_level);
               for (auto &p : this->data[i][j][k])
                  if (condition(p))
                     buckets[i][j][k]->add(p->X, p);

               buckets[i][j][k]->generateTree(condition);
               buckets[i][j][k]->has_tree = true;
            }
      this->has_tree = true;
   }

   /* =================== Access to the bucket at the level of the tree =================== */

   void applyToBucket(const int levelIN, const std::function<void(std::shared_ptr<Buckets<T>>)> &func_for_bucket) {

      if (!this->has_tree) {
         // std::cout << Blue << "max_level : " << this->max_level << std::endl;
         // std::cout << Blue << "level : " << level << std::endl;
         // std::cout << Blue << "This bucket does not have a tree structure." << colorReset << std::endl;
         return;
      }

      if (this->max_level < levelIN) {
         // std::cout << Magenta << "max_level : " << this->max_level << std::endl;
         // std::cout << Magenta << "level : " << level << std::endl;
         // std::cout << Magenta << "The level of the tree is less than the specified level." << colorReset << std::endl;
         return;
      }

      for (auto i = 0; i < this->buckets.size(); ++i)
         for (auto j = 0; j < this->buckets[i].size(); ++j)
            for (auto k = 0; k < this->buckets[i][j].size(); ++k) {
               //! this->bucketsにはlevel+1のバケツが格納されている．
               auto &B = this->buckets[i][j][k];
               std::cout << "i, j, k : " << i << ", " << j << ", " << k << std::endl;
               if (!B->has_tree)
                  continue;
               else if (B->level == levelIN) {
                  std::cout << "level : " << levelIN << std::endl;
                  func_for_bucket(B);
               } else
                  B->applyToBucket(level, func_for_bucket);
            };
   }

   void applyToData(const int levelIN, const std::function<void(std::vector<T>)> &func_for_bucket) {

      if (!this->has_tree) {
         // std::cout << Blue << "max_level : " << this->max_level << std::endl;
         // std::cout << Blue << "level : " << level << std::endl;
         // std::cout << Blue << "This bucket does not have a tree structure." << colorReset << std::endl;
         return;
      }

      if (this->max_level < levelIN) {
         // std::cout << Magenta << "max_level : " << this->max_level << std::endl;
         // std::cout << Magenta << "level : " << level << std::endl;
         // std::cout << Magenta << "The level of the tree is less than the specified level." << colorReset << std::endl;
         return;
      }

      for (auto i = 0; i < this->buckets.size(); ++i)
         for (auto j = 0; j < this->buckets[i].size(); ++j)
            for (auto k = 0; k < this->buckets[i][j].size(); ++k) {
               //! this->bucketsにはlevel+1のバケツが格納されている．
               auto &B = this->buckets[i][j][k];
               std::cout << "i, j, k : " << i << ", " << j << ", " << k << std::endl;
               if (!B->has_tree)
                  continue;
               else if (B->level == levelIN) {
                  std::cout << "level : " << levelIN << std::endl;
                  func_for_bucket(this->data[i][j][k]);
               } else
                  B->applyToData(level, func_for_bucket);
            };
   }

   void applyToDataAndBucket(const int levelIN, const std::function<void(std::unordered_set<T>, std::shared_ptr<Buckets<T>>)> &func_for_bucket) {

      if (!this->has_tree) {
         // std::cout << Blue << "max_level : " << this->max_level << std::endl;
         // std::cout << Blue << "level : " << level << std::endl;
         // std::cout << Blue << "This bucket does not have a tree structure." << colorReset << std::endl;
         return;
      }

      if (this->max_level < levelIN) {
         // std::cout << Magenta << "max_level : " << this->max_level << std::endl;
         // std::cout << Magenta << "level : " << level << std::endl;
         // std::cout << Magenta << "The level of the tree is less than the specified level." << colorReset << std::endl;
         return;
      }

      for (auto i = 0; i < this->buckets.size(); ++i)
         for (auto j = 0; j < this->buckets[i].size(); ++j)
            for (auto k = 0; k < this->buckets[i][j].size(); ++k) {
               //! this->bucketsにはlevel+1のバケツが格納されている．
               auto &B = this->buckets[i][j][k];
               std::cout << "i, j, k : " << i << ", " << j << ", " << k << std::endl;
               if (!B->has_tree)
                  continue;
               else if (B->level == levelIN) {
                  std::cout << "level : " << levelIN << std::endl;
                  func_for_bucket(this->data[i][j][k], this->buckets[i][j][k]);
               } else
                  B->applyToDataAndBucket(level, func_for_bucket);
            };
   }

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
   std::unordered_map<T, ST3> map_to_ijk;
   double dL;
   double bucketVolume() const { return std::pow(this->dL, 3.); };
   //
   Buckets() = default;
   Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : CoordinateBounds(c_bounds) {
      initialize(this->bounds, dL_IN);
      this->center = this->X;
   };
   Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : CoordinateBounds(boundingboxIN) {
      initialize(this->bounds, dL_IN);
      this->center = this->X;
   };

   void initialize(const auto &boundingboxIN, const double dL_IN) {
      // Clear existing data
      this->all_stored_objects.clear();

      // Set bounds and dimensions
      CoordinateBounds::setBounds(boundingboxIN);
      this->dL = dL_IN;

      // check if the size is valid (lesss than 100000)
      const std::size_t max_size = 1000000;
      {
         double s = std::ceil((std::get<1>(this->xbounds()) - std::get<0>(this->xbounds())) / this->dL);
         if (s <= 0.)
            s = 1.;
         if (s > max_size) {
            std::cout << "xsize : " << s << std::endl;
            std::cout << "this->xbounds() : " << this->xbounds() << std::endl;
            std::cout << "this->ybounds() : " << this->ybounds() << std::endl;
            std::cout << "this->zbounds() : " << this->zbounds() << std::endl;
            std::cout << "this->dL : " << this->dL << std::endl;
            std::cout << "The size of the bucket is too large. Please reduce the size of the bucket." << std::endl;
            std::cout << "The size of the bucket is " << s << std::endl;
            std::cout << "The size of the bucket should be less than " << max_size << "." << std::endl;
            std::cout << "The program will be terminated." << std::endl;
            exit(1);
         } else
            this->xsize = s;
      }

      {
         double s = std::ceil((std::get<1>(this->ybounds()) - std::get<0>(this->ybounds())) / this->dL);
         if (s <= 0.)
            s = 1.;
         if (s > max_size) {
            std::cout << "ysize : " << s << std::endl;
            std::cout << "this->xbounds() : " << this->xbounds() << std::endl;
            std::cout << "this->ybounds() : " << this->ybounds() << std::endl;
            std::cout << "this->zbounds() : " << this->zbounds() << std::endl;
            std::cout << "this->dL : " << this->dL << std::endl;
            std::cout << "The size of the bucket is too large. Please reduce the size of the bucket." << std::endl;
            std::cout << "The size of the bucket is " << s << std::endl;
            std::cout << "The size of the bucket should be less than " << max_size << "." << std::endl;
            std::cout << "The program will be terminated." << std::endl;
            exit(1);
         } else
            this->ysize = s;
      }
      {
         double s = std::ceil((std::get<1>(this->zbounds()) - std::get<0>(this->zbounds())) / this->dL);
         if (s <= 0.)
            s = 1.;  // 以前，0になったためエラーになった．
         if (s > max_size) {
            std::cout << "zsize : " << s << std::endl;
            std::cout << "this->xbounds() : " << this->xbounds() << std::endl;
            std::cout << "this->ybounds() : " << this->ybounds() << std::endl;
            std::cout << "this->zbounds() : " << this->zbounds() << std::endl;
            std::cout << "The size of the bucket is too large. Please reduce the size of the bucket." << std::endl;
            std::cout << "The size of the bucket is " << s << std::endl;
            std::cout << "The size of the bucket should be less than " << max_size << "." << std::endl;
            std::cout << "The program will be terminated." << std::endl;
            exit(1);
         } else
            this->zsize = s;
      }

      this->dn = {xsize, ysize, zsize};

      // Clear and resize data
      this->data.assign(xsize, std::vector<std::vector<std::unordered_set<T>>>(ysize, std::vector<std::unordered_set<T>>(zsize, std::unordered_set<T>{})));
      this->data_bool.assign(xsize, std::vector<std::vector<bool>>(ysize, std::vector<bool>(zsize, false)));

      // Set vector_is_set to false to indicate data are not yet set
      this->vector_is_set = false;

      // Output information for debugging
      std::cout << "1 bucket width = " << this->dL << std::endl;  // Changed the comment to English
      std::cout << "Bounds = " << this->bounds << std::endl;
      std::cout << "Size = " << this->dn << std::endl;
   }

   const std::unordered_set<T> &getAll() const { return this->all_stored_objects; };
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
   bool isInside(const ST i, const ST j, const ST k) const { return (i >= 0 && j >= 0 && k >= 0 && i < this->xsize && j < this->ysize && k < this->zsize); };
   bool isInside(const ST3 &ijk) const { return isInside(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk)); };
   bool isInside(const Tddd &x) const { return isInside(indices_no_clamp(x)); };
   void erase(T const p) {
      auto it = this->map_to_ijk.find(p);
      if (it != this->map_to_ijk.end()) {
         auto [i, j, k] = it->second;
         this->all_stored_objects.erase(p);
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

      std::cout << "this->xbounds() : " << this->xbounds() << std::endl;
      std::cout << "this->ybounds() : " << this->ybounds() << std::endl;
      std::cout << "this->zbounds() : " << this->zbounds() << std::endl;
      std::cout << "this->dL : " << this->dL << std::endl;
      std::cout << "this->xsize : " << this->xsize << std::endl;
      std::cout << "this->ysize : " << this->ysize << std::endl;
      std::cout << "this->zsize : " << this->zsize << std::endl;
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
      return bucket_inserted && all_objects_inserted;
   };
   bool add(const Tddd &x, const T p) {
      this->vector_is_set = false;
      return add(indices(x), p);
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
      for (const auto &X : Subdivide(A, B, std::ceil(Norm(A - B) / this->dL)))
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
   ret.buckets.resize(buckets.xsize, std::vector<std::vector<std::shared_ptr<Buckets<T>>>>(buckets.ysize, std::vector<std::shared_ptr<Buckets<T>>>(buckets.zsize, nullptr)));
   ret.level = buckets.level;
   ret.max_level = buckets.max_level;
   ret.has_tree = buckets.has_tree;
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

#endif