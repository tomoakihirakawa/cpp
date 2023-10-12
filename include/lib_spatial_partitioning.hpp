#ifndef lib_spatial_partitioning_H
#define lib_spatial_partitioning_H

/*DOC_EXTRACT 0_3_space_partitioning

## Bucket クラスの説明

このクラスは，オブジェクトを３次元空間内に配置し，効率的に検索できるようにするための「バケツ（Bucket）」構造を提供します．

WARNING: テンプレート型`T`のオブジェクトは，`getX()`でxyz座標を取得できる必要があります．

### 型エイリアス

| エイリアス   | 説明                                         |
|:------------:|:--------------------------------------------:|
| `sizeType`   | サイズ型（int）                              |
| `ST`         | sizeTypeの別名                               |
| `ST2`        | 2次元サイズ配列（std::array<sizeType, 2>）   |
| `ST3`        | 3次元サイズ配列（std::array<sizeType, 3>）   |
| `ST6`        | 6次元サイズ配列（std::array<sizeType, 6>）   |

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

#### インデックス変換

- `itox(const ST i, const ST j, const ST k) const`: インデックスから座標へ変換．
- `indices(const Tddd &x) const`: 座標からインデックスへ変換．

#### データ追加・削除

- `add(const Tddd &x, const T p)`: オブジェクトを追加．
- `erase(T const p)`: オブジェクトを削除．

#### その他

- `none_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const`: 条件に合うオブジェクトがないか確認．

### 使用例

```cpp
// 座標の境界を定義
T3Tdd bounding_box = {{{0.0, 10.0}, {0.0, 10.0}, {0.0, 10.0}}};

// バケツの一辺の長さを定義
double bucket_edge_length = 1.0;

// Bucket インスタンスを初期化
BaseBuckets<MyObject> my_buckets(bounding_box, bucket_edge_length);

// オブジェクトの座標と初期値を定義
Tddd obj1_coordinates = {5.0, 5.0, 5.0};
MyObject obj1 = MyObject(some initializers);

Tddd obj2_coordinates = {6.0, 6.0, 6.0};
MyObject obj2 = MyObject(some initializers);

// オブジェクトを追加
my_buckets.add(obj1_coordinates, obj1);
my_buckets.add(obj2_coordinates, obj2);

// オブジェクトを削除
my_buckets.erase(obj1);
```

*/

// Buckets is derived from CoordinateBounds

template <typename T>
struct Buckets : public CoordinateBounds {
   /* -------------------------------------------------------------------------- */
   std::vector<std::vector<std::vector<std::shared_ptr<Buckets<T>>>>> buckets;

   int level = 0;
   int max_level = 2;
   bool has_tree = false;

   void generateTree(const std::function<bool(const std::unordered_set<T> &)> &condition = [](const std::unordered_set<T> &) { return true; }) {
      if (this->level >= this->max_level)
         return;

      buckets.resize(this->xsize, std::vector<std::vector<std::shared_ptr<Buckets<T>>>>(this->ysize, std::vector<std::shared_ptr<Buckets<T>>>(this->zsize, nullptr)));

      for (auto i = 0; i < this->data.size(); ++i)
         for (auto j = 0; j < this->data[i].size(); ++j)
            for (auto k = 0; k < this->data[i][j].size(); ++k) {
               auto bounds = getBounds({i, j, k});
               //! この方法で，`data[i][j][k]`を８分割しバケツを作成する．\label{buckets_generateTree}
               buckets[i][j][k] = std::make_shared<Buckets<T>>(bounds, this->dL * 0.5 + 1e-10);
               buckets[i][j][k]->level = this->level + 1;
               for (const auto &p : this->data[i][j][k]) {
                  // if (condition(p))
                  {
                     buckets[i][j][k]->add(ToX(p), p);
                  }
               }
            }

      this->has_tree = true;
   }
   /* -------------------------------------------------------------------------- */
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
   Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : CoordinateBounds(c_bounds) { initialize(this->bounds, dL_IN); };
   Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : CoordinateBounds(boundingboxIN) { initialize(this->bounds, dL_IN); };

   void initialize(const auto &boundingboxIN, const double dL_IN) {
      // Clear existing data
      this->all_stored_objects.clear();

      // Set bounds and dimensions
      CoordinateBounds::setBounds(boundingboxIN);
      this->dL = dL_IN;

      // Calculate size based on bounds and resolution
      this->xsize = std::ceil((std::get<1>(this->xbounds()) - std::get<0>(this->xbounds())) / this->dL);
      this->ysize = std::ceil((std::get<1>(this->ybounds()) - std::get<0>(this->ybounds())) / this->dL);
      this->zsize = std::ceil((std::get<1>(this->zbounds()) - std::get<0>(this->zbounds())) / this->dL);
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
      return {static_cast<ST>(std::floor((std::get<0>(x) - std::get<0>(this->xbounds())) / this->dL)),
              static_cast<ST>(std::floor((std::get<1>(x) - std::get<0>(this->ybounds())) / this->dL)),
              static_cast<ST>(std::floor((std::get<2>(x) - std::get<0>(this->zbounds())) / this->dL))};
   };
   constexpr ST3 indices(const Tddd &x) const {
      /*
      intキャストはゼロ方向へ実数を切り捨てた結果を返すので，static_cast<int>によって正しくセルのインデックスに変換できる．
       0.**    1.**  2.**   3.**
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
      for (const auto &p : P) {
         auto ijk = indices(ToX(p));
         this->map_to_ijk[p] = ijk;
         ret = ret && add(ijk, p);
      }
      return ret;
   };
   /* -------------------------------------------------------------------------- */
   auto getBucket(const Tddd &x) const {
      const auto ijk = this->indices(x);
      return this->data[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
   };

   auto getBucket(const Tddd &x, const double &range) const {
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
         return !std::any_of(this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::any_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::any_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  return std::any_of(Bijk.cbegin(), Bijk.cend(), func);
               });
            });
         });
      } else {
         return !std::any_of(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::any_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::any_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  return std::any_of(Bijk.cbegin(), Bijk.cend(), func);
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
         return std::all_of(this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::all_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::all_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  return std::all_of(Bijk.cbegin(), Bijk.cend(), func);
               });
            });
         });
      } else {
         return std::all_of(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::all_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::all_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  return std::all_of(Bijk.cbegin(), Bijk.cend(), func);
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
         return std::any_of(this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::any_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::any_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  return std::any_of(Bijk.cbegin(), Bijk.cend(), func);
               });
            });
         });
      } else {
         return std::any_of(this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&](const auto &Bi) {
            return std::any_of(Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&](const auto &Bij) {
               return std::any_of(Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&](const auto &Bijk) {
                  return std::any_of(Bijk.cbegin(), Bijk.cend(), func);
               });
            });
         });
      }
   }

   //! apply
   void apply(const Tddd &x, const double d, const std::function<void(const T &)> &func) const {
      if (this->data.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");

      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);

      if (!this->vector_is_set) {
         std::for_each(std::execution::unseq, this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
            std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
               std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
                  for (const auto &p : Bijk) func(p);
               });
            });
         });
      } else {
         std::for_each(std::execution::unseq, this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
            std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
               std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
                  for (const auto &p : Bijk) func(p);
               });
            });
         });
      }
   };

   //! apply
   void apply(const Tddd &x, const double d, const std::function<void(const int, const int, const int)> &func) const {
      if (this->data.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");
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
   void apply(const Tddd &x, const Tddd dx_dy_dz,
              const std::function<void(const T &)> &func) const {
      if (this->data.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D data is empty");

      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, dx_dy_dz);

      if (!this->vector_is_set) {
         std::for_each(std::execution::unseq, this->data.cbegin() + i_min, this->data.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
            std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
               std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
                  for (const auto &p : Bijk) func(p);
               });
            });
         });
      } else {
         std::for_each(std::execution::unseq, this->data_vector.cbegin() + i_min, this->data_vector.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
            std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
               std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
                  for (const auto &p : Bijk) func(p);
               });
            });
         });
      }
   };
};

// template <typename T>
// struct Buckets : public BaseBuckets<T> {
//    Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<T>(c_bounds, dL_IN){};
//    Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<T>(boundingboxIN, dL_IN){};
// };

#endif