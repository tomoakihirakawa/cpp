#ifndef lib_spatial_partitioning_H
#define lib_spatial_partitioning_H

// b$ ------------------------------------------------------ */
// b$         近傍のオブジェクト探査のための，空間分割バケツ         */
// b$ ------------------------------------------------------ */
//  テンプレートどんなオブジェクトでもできるはず
//   #define debug_BaseBuckets

template <typename T>
struct BaseBuckets {
   using sizeType = size_t;
   using ST = sizeType;
   using ST2 = std::array<sizeType, 2>;
   using ST3 = std::array<sizeType, 3>;
   using ST6 = std::array<sizeType, 6>;
   //! getX()でxyz座標を取得できるオブジェクトTのための，バケツ
   Tdd xbounds, ybounds, zbounds;
   ST xsize, ysize, zsize;
   T3Tdd bounds;
   Tddd center;
   ST3 dn;
   // bool hashing_done;
   std::vector</*x*/ std::vector</*y*/ std::vector</*z*/ std::unordered_set<T>>>> buckets;
   // std::unordered_map<Tiii, std::unordered_set<T> *> hashed_buckets;
   std::unordered_set<T> all_stored_objects;
   std::unordered_map<T, ST3> map_to_ijk;
   /*
   bucketsベクトルは，unordered_setで重複を認めないが，他のベクトルには，同じオブジェクトが入っていてもよい．
   例えば，鏡像関係を作りたい場合，他のバケット位置に，同じオブジェクトを入れておけばいい．
   */
   double dL;
   double bucketVolume() const { return std::pow(this->dL, 3.); };
   double bucketVolume(const ST depth) const {
      //    depth  1  2  3  4
      // edge len  1  3  5  7   = 2*d-1
      //  buckets  1  9  25  49 = (2*d-1)^2
      auto vol = std::pow(this->dL, 3.);
      return pow(2. * depth - 1., 2.) * vol;
   };
   BaseBuckets(const CoordinateBounds &c_bounds, const double dL_IN) : bounds(c_bounds.bounds), dL(dL_IN), dn(ST3{0, 0, 0}) {
      initialize(this->bounds, dL_IN);
   };
   BaseBuckets(const T3Tdd &boundingboxIN, const double dL_IN) : bounds(boundingboxIN), dL(dL_IN), dn(ST3{0, 0, 0}) {
      initialize(this->bounds, dL_IN);
   };

   void initialize(const T3Tdd &boundingboxIN, const double dL_IN) {
      this->buckets.clear();
      this->all_stored_objects.clear();
      //
      this->bounds = boundingboxIN;
      this->dL = dL_IN;
      this->xbounds = std::get<0>(this->bounds);
      this->ybounds = std::get<1>(this->bounds);
      this->zbounds = std::get<2>(this->bounds);
      this->center = {Mean(this->xbounds), Mean(this->ybounds), Mean(this->zbounds)};
      this->xsize = std::ceil((std::get<1>(this->xbounds) - std::get<0>(this->xbounds)) / this->dL);
      this->ysize = std::ceil((std::get<1>(this->ybounds) - std::get<0>(this->ybounds)) / this->dL);
      this->zsize = std::ceil((std::get<1>(this->zbounds) - std::get<0>(this->zbounds)) / this->dL);
      this->dn = {xsize, ysize, zsize};
      this->buckets.resize(xsize, std::vector</*y*/ std::vector</*z*/ std::unordered_set<T>>>(ysize, std::vector</*z*/ std::unordered_set<T>>(zsize, std::unordered_set<T>{})));

      std::cout << "１バケット幅 = " << this->dL << std::endl;
      std::cout << "バウンド = " << this->bounds << std::endl;
      std::cout << "サイズ = " << this->dn << std::endl;
   }
   const std::unordered_set<T> &getAll() const { return this->all_stored_objects; };

   //@ -------------------------------- インデックス変換 -------------------------------- */
   // x座標を内包するバケツのインデックスを返す
   Tddd itox(const auto i, const auto j, const auto k) const {
      return {this->dL * 0.5 + this->dL * i + std::get<0>(this->xbounds),
              this->dL * 0.5 + this->dL * j + std::get<0>(this->ybounds),
              this->dL * 0.5 + this->dL * k + std::get<0>(this->zbounds)};
   };
   Tddd itox(const ST3 &ijk) const { return itox(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk)); };
   ST3 indices(const Tddd &x) const {
      // uint32_tのキャストはゼロ方向へ実数を切り捨てた結果を返す
      /*
       0.**    1.**  2.**   3.**  <= X-minX というわけで，static_cast<uint32_t>によって正しくセルのインデックスに変換できる
      <-dL-> <-dL-> <-dL-> <-dL->
      *-----*------*------*------*
      |  0  |   1  |   2  |   3  |
      *-----*------*------*------*
      */
      return {std::clamp(static_cast<ST>((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL), static_cast<ST>(0), this->xsize - 1),
              std::clamp(static_cast<ST>((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL), static_cast<ST>(0), this->ysize - 1),
              std::clamp(static_cast<ST>((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL), static_cast<ST>(0), this->zsize - 1)};
   };
   ST6 indices_range(const Tddd &x, const double d) const {
      auto [i_min, j_min, k_min] = indices(x - d);
      auto [i_max, j_max, k_max] = indices(x + d);
      return {i_min, i_max, j_min, j_max, k_min, k_max};
   };
   // void indices(const Tddd &x, ST &i, ST &j, ST &k) const {
   //    i = static_cast<ST>((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL);
   //    j = static_cast<ST>((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL);
   //    k = static_cast<ST>((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL);
   // };
   //@ ------------------------ インデックスがboundsに収まっているかどうか ------------------------ */
   bool isInside(const ST i, const ST j, const ST k) const {
      return (i >= 0 && j >= 0 && k >= 0 && i < this->xsize && j < this->ysize && k < this->zsize);
   };

   bool isInside(const ST3 &ijk) const {
      auto [i, j, k] = ijk;
      return isInside(i, j, k);
   };

   bool isInside(const Tddd &x) const {
      return isInside(indices(x));
   };

   void erase(T const p) {
      auto it = this->map_to_ijk.find(p);
      if (it != this->map_to_ijk.end()) {
         auto [i, j, k] = it->second;
         this->all_stored_objects.erase(p);
         this->buckets[i][j][k].erase(p);
         this->map_to_ijk.erase(it);
      }
   };
   //@ -------------------------------------------------------------------------- */
   //! explicitly specify coordinates
   bool add(const ST3 &ijk, const T p) {
      // hashing_done = false;
      // auto [i, j, k] = ijk;
      // try {
      //    if (!(i < 0 || j < 0 || k < 0 || i >= this->xsize || j >= this->ysize || k >= this->zsize)) {
      //       this->map_to_ijk[p] = ijk;
      //       return (this->buckets[i][j][k].emplace(p)).second && (this->all_stored_objects.emplace(p)).second;
      //    }
      // } catch (std::exception &e) {
      //    std::cerr << e.what() << colorOff << std::endl;
      //    std::stringstream ss;
      //    ss << "[i,j,k] = [" << i << "," << j << "," << k << "]" << std::endl;
      //    ss << "[xsize,ysize,zsize] = [" << xsize << "," << ysize << "," << zsize << "]" << std::endl;
      //    throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      // };
      // return false;
      // hashing_done = false;
      auto [i, j, k] = ijk;
      if (isInside(i, j, k)) {
         this->map_to_ijk[p] = ijk;
         bool bucket_inserted = this->buckets[i][j][k].emplace(p).second;
         bool all_objects_inserted = this->all_stored_objects.emplace(p).second;
         return bucket_inserted && all_objects_inserted;
      }
      return false;
   };
   bool add(const Tddd &x, const T p) { return add(indices(x), p); };

   //! automatically specify coordinates
   bool add(const std::unordered_set<T> &P) {
      // hashing_done = false;
      bool ret = true;
      auto func0 = [&]() {
               ST3 ijk;
               for (const auto &p : P) {
                  ijk = indices(ToX(p));
                  auto [i, j, k] = ijk;
                  this->map_to_ijk[p] = ijk;
               }; };
      auto func1 = [&]() {
               for (const auto &p : P) {
               auto [i, j, k] = indices(ToX(p));
               ret = ret && (this->buckets[i][j][k].emplace(p)).second;
               ret = ret && (this->all_stored_objects.emplace(p)).second;
            } };
#pragma omp parallel sections
      {
#pragma omp section
         func0();
#pragma omp section
         func1();
      };
      return ret;
   };
   bool add(const std::vector<T> &p) {
      bool ret = true;
      for (const auto &q : p)
         ret = ret && add(ToX(q), q);
      return ret;
   };
   /* -------------------------------------------------------------------------- */

   // 座標を入力し，バケツを指定する．
   std::vector<T> getAllBuckets() const {
      return this->buckets;
   };
   std::vector<T4Tddd> getT4Tddd() const {
      std::vector<T4Tddd> ret;
      for (auto k = 0; k < this->zsize; ++k)
         for (auto j = 0; j < this->ysize; ++j)
            for (auto i = 0; i < this->xsize; ++i) {
               Tddd X0 = {(i * this->dL + std::get<0>(this->xbounds)),
                          (j * this->dL + std::get<0>(this->ybounds)),
                          (k * this->dL + std::get<0>(this->zbounds))};
               Tddd X1 = {((i + 1) * this->dL + std::get<0>(this->xbounds)),
                          (j * this->dL + std::get<0>(this->ybounds)),
                          (k * this->dL + std::get<0>(this->zbounds))};
               Tddd X2 = {((i + 1) * this->dL + std::get<0>(this->xbounds)),
                          ((j + 1) * this->dL + std::get<0>(this->ybounds)),
                          (k * this->dL + std::get<0>(this->zbounds))};
               Tddd X3 = {(i * this->dL + std::get<0>(this->xbounds)),
                          ((j + 1) * this->dL + std::get<0>(this->ybounds)),
                          (k * this->dL + std::get<0>(this->zbounds))};
               Tddd X4 = {(i * this->dL + std::get<0>(this->xbounds)),
                          (j * this->dL + std::get<0>(this->ybounds)),
                          ((k + 1) * this->dL + std::get<0>(this->zbounds))};
               Tddd X5 = {((i + 1) * this->dL + std::get<0>(this->xbounds)),
                          (j * this->dL + std::get<0>(this->ybounds)),
                          ((k + 1) * this->dL + std::get<0>(this->zbounds))};
               Tddd X6 = {((i + 1) * this->dL + std::get<0>(this->xbounds)),
                          ((j + 1) * this->dL + std::get<0>(this->ybounds)),
                          ((k + 1) * this->dL + std::get<0>(this->zbounds))};
               Tddd X7 = {(i * this->dL + std::get<0>(this->xbounds)),
                          ((j + 1) * this->dL + std::get<0>(this->ybounds)),
                          ((k + 1) * this->dL + std::get<0>(this->zbounds))};
               ret.push_back({X0, X1, X2, X3});
               ret.push_back({X0, X1, X5, X4});
               ret.push_back({X1, X2, X6, X5});
               ret.push_back({X0, X4, X7, X3});
               ret.push_back({X2, X3, X7, X6});
               ret.push_back({X4, X5, X6, X7});
            }
      return ret;
   };
   std::vector<T> getBucket(const Tddd &x) const {
      Tiii i = this->indices(x);
      return this->buckets[std::get<0>(i)][std::get<1>(i)][std::get<2>(i)];
   };
   /* ------------------------------------------------------ */
   void clear() {
      // hashing_done = false;
      this->buckets.clear();
      this->all_stored_objects.clear();
      // for (auto &vvu_b : this->buckets)
      // 	for (auto &vu_b : vvu_b)
      // 		for (auto &u_b : vu_b)
      // 			u_b.clear();
   };
   //@ 2021年11月9日さらに早くするために，
   //@ std::vector<std::vector<T >> getObjects(const Tddd &x, const int limit_depth /*limit depth*/, const int limit_number = 100000) const
   //@ を改良
   /* -------------- デフォルトのバケツは，深さ毎に粒子を保存していく -------------- */
   // i=3,d=2
   // index = | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
   // depth =     | -2| -1| 0 | 1 | 2 |
   // int bind_range(const int i, const Tii &minmax) const { return (i <= std::get<0>(minmax) ? std::get<0>(minmax) : (i >= std::get<1>(minmax) ? std::get<1>(minmax) : i)); };
   // ST min_index_x(const ST i, const ST d) const { return std::clamp(i - d, static_cast<ST>(0), this->xsize - 1); }
   // ST max_index_x(const ST i, const ST d) const { return std::clamp(i + d, static_cast<ST>(0), this->xsize - 1); }
   // ST min_index_y(const ST i, const ST d) const { return std::clamp(i - d, static_cast<ST>(0), this->ysize - 1); }
   // ST max_index_y(const ST i, const ST d) const { return std::clamp(i + d, static_cast<ST>(0), this->ysize - 1); }
   // ST min_index_z(const ST i, const ST d) const { return std::clamp(i - d, static_cast<ST>(0), this->zsize - 1); }
   // ST max_index_z(const ST i, const ST d) const { return std::clamp(i + d, static_cast<ST>(0), this->zsize - 1); }

   // int min_index_x(const int i, const int d) const { return ((i - d) <= 0 ? 0 : ((i - d) >= (this->xsize - 1) ? this->xsize - 1 : i - d)); };
   // int max_index_x(const int i, const int d) const { return ((i + d) >= (this->xsize - 1) ? this->xsize - 1 : i + d); };
   // int min_index_y(const int i, const int d) const { return ((i - d) <= 0 ? 0 : i - d); };
   // int max_index_y(const int i, const int d) const { return ((i + d) >= (this->ysize - 1) ? this->ysize - 1 : i + d); };
   // int min_index_z(const int i, const int d) const { return ((i - d) <= 0 ? 0 : i - d); };
   // int max_index_z(const int i, const int d) const { return ((i + d) >= (this->zsize - 1) ? this->zsize - 1 : i + d); };

   // std::unordered_set<T> getObjects_unorderedset(const Tddd &x, const ST limit_depth /*limit depth*/, const ST limit_number = 100000) const {
   //    std::unordered_set<T> ret;
   //    // ret.reserve(limit_depth);
   //    // std::vector<T > ret_last(0);
   //    auto [i0, j0, k0] = this->indices(x);
   //    //* ------------------------ depth=0 ------------------------ */
   //    // auto &r = *ret.rbegin();
   //    if (isInside(i0, j0, k0))
   //       ret.insert(this->buckets[i0][j0][k0].begin(), this->buckets[i0][j0][k0].end());
   //    if (ret.size() >= limit_number)
   //       return ret;
   //    /* ------------------------------------------------------ */
   //    ST j_min, j_max, k_min, k_max, i_min_1, i_max_1, j_min_1, j_max_1, i0_m_d, i0_p_d, j0_m_d, j0_p_d, k0_m_d, k0_p_d;
   //    for (auto d = 1; d <= limit_depth; ++d) {
   //       j_min = min_index_y(j0, d);
   //       j_max = max_index_y(j0, d);
   //       k_min = min_index_z(k0, d);
   //       k_max = max_index_z(k0, d);
   //       j_min_1 = min_index_y(j0, d - 1);
   //       j_max_1 = max_index_y(j0, d - 1);
   //       i0_m_d = i0 - d;
   //       i0_p_d = i0 + d;
   //       for (auto j = j_min; j <= j_max /*need equal (=)*/; ++j)
   //          for (auto k = k_min; k <= k_max; ++k) {
   //             if (i0_m_d >= 0 && i0_m_d < this->xsize)
   //                ret.insert(std::begin(this->buckets[i0_m_d][j][k]), std::end(this->buckets[i0_m_d][j][k]));
   //             if (i0_p_d >= 0 && i0_p_d < this->xsize)
   //                ret.insert(std::begin(this->buckets[i0_p_d][j][k]), std::end(this->buckets[i0_p_d][j][k]));
   //          }
   //       i_min_1 = min_index_x(i0, d - 1);
   //       i_max_1 = max_index_x(i0, d - 1);
   //       //
   //       j0_m_d = j0 - d;
   //       j0_p_d = j0 + d;
   //       //
   //       k0_m_d = k0 - d;
   //       k0_p_d = k0 + d;
   //       //
   //       for (auto i = i_min_1; i <= i_max_1; ++i) {
   //          for (auto k = k_min; k <= k_max; ++k) {
   //             if (j0_m_d >= 0 && j0_m_d < this->ysize)
   //                ret.insert(std::begin(this->buckets[i][j0_m_d][k]), std::end(this->buckets[i][j0_m_d][k]));
   //             if (j0_p_d >= 0 && j0_p_d < this->ysize)
   //                ret.insert(std::begin(this->buckets[i][j0_p_d][k]), std::end(this->buckets[i][j0_p_d][k]));
   //          }
   //          for (auto j = j_min_1; j <= j_max_1; ++j) {
   //             if (k0_m_d >= 0 && k0_m_d < this->zsize)
   //                ret.insert(std::begin(this->buckets[i][j][k0_m_d]), std::end(this->buckets[i][j][k0_m_d]));
   //             if (k0_p_d >= 0 && k0_p_d < this->zsize)
   //                ret.insert(std::begin(this->buckets[i][j][k0_p_d]), std::end(this->buckets[i][j][k0_p_d]));
   //          }
   //       }
   //       if (ret.size() >= limit_number)
   //          return ret;
   //    }
   //    return ret;
   // };
   // std::unordered_set<T> getObjects_unorderedset(const Tddd &x, const double smoothing_length, const ST limit_number = 100000) const {
   //    ST depth = std::ceil(smoothing_length / this->dL);
   //    return getObjects_unorderedset(x, depth, limit_number);
   // };
   // b@ -------------------------------------------------------------------------- */
   // b@                             STL like functions                             */
   // b@ -------------------------------------------------------------------------- */
   //! count_if
   ST count_if(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      ST ret = 0;
      if (!this->buckets.empty()) {
         const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_range(x, d);
         auto it = this->buckets.begin();
         auto jt = it->begin();
         auto kt = jt->begin();
         for (it = std::next(this->buckets.begin(), i_min); it != std::next(this->buckets.begin(), i_max); ++it)
            for (jt = std::next(it->begin(), j_min); jt != std::next(it->begin(), j_max); ++jt)
               for (kt = std::next(jt->begin(), k_min); kt != std::next(jt->begin(), k_max); ++kt)
                  for (const auto &p : *kt)
                     if (func(p))
                        ret++;
      }
      return ret;
   };
   ST count_if(const std::function<bool(const T &)> &func) const {
      ST ret = 0;
      for (const auto &p : this->all_stored_objects)
         if (func(p))
            ret++;
      return ret;
   };
   //! none_of
   bool none_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (!this->buckets.empty()) {
         const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_range(x, d);
         /* ------------------------------ イテレータを使ったループ ------------------------------ */
         // std::cout << "i, j, k = " << i << ", " << j << ", " << k << std::endl;
         // std::cout << "d = " << d << std::endl;
         // std::cout << "i = " << i_min << ", " << i_max << ", " << this->xsize << std::endl;
         // std::cout << "j = " << j_min << ", " << j_max << ", " << this->ysize << std::endl;
         // std::cout << "k = " << k_min << ", " << k_max << ", " << this->zsize << std::endl;
         // for (auto i = i_min; i < i_max; ++i)
         //    for (auto j = j_min; j < j_max; ++j)
         //       for (auto k = k_min; k < k_max; ++k)
         //          for (const auto &p : this->buckets[i][j][k])
         //             if (func(p))
         //                return false;
         auto it = this->buckets.begin();
         auto jt = it->begin();
         auto kt = jt->begin();
         for (it = std::next(this->buckets.begin(), i_min); it != std::next(this->buckets.begin(), i_max); ++it)
            for (jt = std::next(it->begin(), j_min); jt != std::next(it->begin(), j_max); ++jt)
               for (kt = std::next(jt->begin(), k_min); kt != std::next(jt->begin(), k_max); ++kt)
                  for (const auto &p : *kt)
                     if (func(p))
                        return false;
      }
      return true;
   };
   bool none_of(const std::function<bool(const T &)> &func) const {
      for (const auto &p : this->all_stored_objects)
         if (func(p))
            return false;
      return true;
   };
   //! all_of
   bool all_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (!this->buckets.empty()) {
         const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_range(x, d);
         /* ------------------------------ イテレータを使ったループ ------------------------------ */
         auto it = this->buckets.begin();
         auto jt = it->begin();
         auto kt = jt->begin();
         for (it = std::next(this->buckets.begin(), i_min); it != std::next(this->buckets.begin(), i_max); ++it)
            for (jt = std::next(it->begin(), j_min); jt != std::next(it->begin(), j_max); ++jt)
               for (kt = std::next(jt->begin(), k_min); kt != std::next(jt->begin(), k_max); ++kt)
                  for (const auto &p : *kt)
                     if (!func(p))
                        return false;
      }
      return true;
   };
   //! any_of
   bool any_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (!this->buckets.empty()) {
         const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_range(x, d);
         /* -------------------------------------------------------------------------- */
         auto it = this->buckets.begin();
         auto jt = it->begin();
         auto kt = jt->begin();
         for (it = std::next(this->buckets.begin(), i_min); it != std::next(this->buckets.begin(), i_max); ++it)
            for (jt = std::next(it->begin(), j_min); jt != std::next(it->begin(), j_max); ++jt)
               for (kt = std::next(jt->begin(), k_min); kt != std::next(jt->begin(), k_max); ++kt)
                  for (const auto &p : *kt)
                     if (func(p) /*一つでも見つかったらtrue*/)
                        return true;
      }
      return false;
   };
   //! apply
   void apply(const Tddd &x, const double d, const std::function<void(const T &)> &func) const {
      if (this->buckets.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D buckets is empty");
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_range(x, 2. * d);
      /* -------------------------------- simple traditional for loop -------------------------------- */
      // for (auto i = i_min; i <= i_max; ++i)
      //    for (auto j = j_min; j <= j_max; ++j)
      //       for (auto k = k_min; k <= k_max; ++k)
      //          for (const auto &p : this->buckets[i][j][k])
      //             func(p);
      /* ------------------------- for loop with iterator ------------------------- */
      // auto it = this->buckets.begin();
      // auto jt = it->begin();
      // auto kt = jt->begin();
      // for (it = std::next(this->buckets.begin(), i_min); it != std::next(this->buckets.begin(), i_max); ++it) {
      //    for (jt = std::next(it->begin(), j_min); jt != std::next(it->begin(), j_max); ++jt) {
      //       for (kt = std::next(jt->begin(), k_min); kt != std::next(jt->begin(), k_max); ++kt) {
      //          for (const auto &p : *kt) func(p);
      //       }
      //    }
      // }
      /* ----------------------------- for_each with next---------------------------- */
      // std::for_each(std::execution::unseq, std::next(this->buckets.cbegin(), i_min), std::next(this->buckets.cbegin(), i_max), [&](const auto &Bi) {
      //    std::for_each(std::execution::unseq, std::next(Bi.cbegin(), j_min), std::next(Bi.cbegin(), j_max), [&](const auto &Bij) {
      //       std::for_each(std::execution::unseq, std::next(Bij.cbegin(), k_min), std::next(Bij.cbegin(), k_max), [&](const auto &Bijk) {
      //          for (const auto &p : Bijk) func(p);
      //       });
      //    });
      // });
      /* ----------------------------- for_each with addition ---------------------------- */
      std::for_each(std::execution::unseq, this->buckets.cbegin() + i_min, this->buckets.cbegin() + i_max + 1, [&func, &j_min, &j_max, &k_min, &k_max](const auto &Bi) {
         std::for_each(std::execution::unseq, Bi.cbegin() + j_min, Bi.cbegin() + j_max + 1, [&func, &k_min, &k_max](const auto &Bij) {
            std::for_each(std::execution::unseq, Bij.cbegin() + k_min, Bij.cbegin() + k_max + 1, [&func](const auto &Bijk) {
               for (const auto &p : Bijk) func(p);
            });
         });
      });
   };

   void apply(const std::function<void(const T &)> &func) const {
      for (const auto &p : this->all_stored_objects) func(p);
   };
   //! ------------------- バケツ内のオブジェクトではなく，バケツを引数として受け取る関数を渡す ------------------- */
   void Apply(const Tddd &x, const ST d, const std::function<void(const std::unordered_set<T> &)> &func) const {
      if (this->buckets.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D buckets is empty");
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_range(x, d);
      /* ----------------------------- for_eachを使ったループ ---------------------------- */
      for_each(std::execution::unseq, std::next(this->buckets.cbegin(), i_min), std::next(this->buckets.cbegin(), i_max), [&](const auto &Bi) {
         for_each(std::execution::unseq, std::next(Bi.cbegin(), j_min), std::next(Bi.cbegin(), j_max), [&](const auto &Bij) {
            for_each(std::execution::unseq, std::next(Bij.cbegin(), k_min), std::next(Bij.cbegin(), k_max), [&](const auto &Bijk) { func(Bijk); });
         });
      });
   };
   void Apply(const Tddd &x, const double d, const std::function<void(const std::unordered_set<T> &)> &func) const {
      Apply(x, static_cast<ST>(std::ceil(d / this->dL)), func);
   };
   /* -------------------------------------------------------------------------- */
   std::vector<T> copy_if(const Tddd &x, const ST d, const std::function<bool(const T &)> &func) const {
      std::vector<T> ret;
      ret.reserve(1000);
      this->apply(x, d, [&](const auto &p) {
         if (func(p)) ret.emplace_back(p);
      });
      return ret;
   };

   std::vector<T> get(const Tddd &x, const ST d) const {
      std::vector<T> ret;
      ret.reserve(1000);
      this->Apply(x, d, [&](const auto &Bijk) { ret.insert(ret.end(), Bijk.begin(), Bijk.end()); });
      return ret;
   };
   std::vector<T> get(const Tddd &x, const double d) const { return get(x, static_cast<ST>(std::ceil(d / this->dL))); };

   std::vector<std::vector<T>> get(const Tddd &x, const Tii &depth_minmax) const {
      std::vector<std::vector<T>> ret(1 + std::get<1>(depth_minmax) - std::get<0>(depth_minmax));
      for (auto i = std::get<0>(depth_minmax); i <= std::get<1>(depth_minmax); ++i)
         ret[i] = this->get(x, i);
      return ret;
   };
   std::vector<std::vector<T>> get(const Tddd &x, const Tdd &depth_minmax) const {
      return get(x, Tii{(ST)std::ceil(std::get<0>(depth_minmax) / this->dL),
                        (ST)std::ceil(std::get<1>(depth_minmax) / this->dL)});
   };

   void getNearest(const std::unordered_set<T> &list, const Tddd &x, T &ret) const {
      if (!list.empty()) {
         double distance = 1E+20;
         for (const auto &l : list) {
            if (distance > Norm(ToX(l) - x)) {
               distance = Norm(ToX(l) - x);
               ret = l;
            }
         }
      }
   };

   void getNearest(const std::unordered_set<T> &list, const T &p, T &ret) const {
      if (!list.empty()) {
         // リスト内で最も近いものを返す
         auto x = ToX(p);
         double distance = 1E+20;
         if (ret)
            distance = Norm(ToX(ret) - x);
         double tmp;
         for (const auto &l : list)
            if (p != l && p && l) {
               tmp = Norm(ToX(l) - x);
               if (distance > tmp) {
                  distance = tmp;
                  ret = l;
               }
            }
      }
   };

   void getNearest(const T &p, const ST d, T &r) const {
      // あるレベルdにおいて最も近い物を返す
      ST i0, j0, k0;
      auto x = ToX(p);
      this->indices(x, i0, j0, k0);
      if (d == 0) {
         if (isInside(i0, j0, k0))
            getNearest(this->buckets[i0][j0][k0], p, r);
      } else {
         ST i, j, k;
         for (j = (j0 - d < 0 ? 0 : j0 - d); j <= j0 + d && j < this->ysize; ++j)
            for (k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; ++k) {
               i = i0 + d;
               if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
                  getNearest(this->buckets[i][j][k], p, r);
               i = i0 - d;
               if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
                  getNearest(this->buckets[i][j][k], p, r);
            }
         for (i = (i0 - d + 1 < 0 ? 0 : i0 - d + 1); i < i0 + d && i < this->xsize; ++i) {
            for (k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; ++k) {
               j = j0 + d;
               if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
                  getNearest(this->buckets[i][j][k], p, r);
               j = j0 - d;
               if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
                  getNearest(this->buckets[i][j][k], p, r);
            }
            for (j = (j0 - d + 1 < 0 ? 0 : j0 - d + 1); j < j0 + d && j < this->ysize; ++j) {
               k = k0 + d;
               if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
                  getNearest(this->buckets[i][j][k], p, r);
               k = k0 - d;
               if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
                  getNearest(this->buckets[i][j][k], p, r);
            }
         }
      }
   };

   void getNearest(const T &p, const Tii &minmaxdepth, T &r) const {
      double current_nearest_distance = 1E+10;
      auto x = ToX(p);
      for (auto d = std::get<0>(minmaxdepth); d < std::get<1>(minmaxdepth); ++d) {
         if (d == 0 || d == 1) {
            getNearest(p, d, r);
            if (r)
               current_nearest_distance = Norm(ToX(r) - x);
         } else {
            if (dL * (d - 1) > current_nearest_distance)
               break;
            else {
               getNearest(p, d, r);
               if (r)
                  current_nearest_distance = Norm(ToX(r) - x);
            }
         }
      }
   };
   // void getNearest(const T &x, const ST d, T &r) const { getNearest(ToX(x), d, r); };

   std::vector<std::vector<T>> getObjects(const Tddd &x, ST limit_depth /*limit depth*/, ST limit_number = 100000) const {
      std::vector<std::vector<T>> ret(limit_depth);
      // ret.reserve(limit_depth);
      // std::vector<T > ret_last(0);
      auto [i0, j0, k0] = this->indices(x);
      //* ------------------------ depth=0 ------------------------ */
      if (isInside(i0, j0, k0))
         ret[0].insert(ret[0].end(), this->buckets[i0][j0][k0].begin(), this->buckets[i0][j0][k0].end());
      ST tot = ret[0].size();
      if (tot >= limit_number)
         return ret;
      /* ---------------------------------------------------------- */
      ST i, j, k;
      for (auto d = 1; d < limit_depth; d++) {
         auto &r = ret[d];
         for (auto j = (j0 - d < 0 ? 0 : j0 - d); j <= j0 + d && j < this->ysize; j++)
            for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; k++) {
               i = i0 + d;
               if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
                  ret[d].insert(ret[d].begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
               i = i0 - d;
               if (!(i < 0 || i >= this->xsize || this->buckets[i][j][k].empty()))
                  ret[d].insert(ret[d].begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
            }
         for (auto i = (i0 - d + 1 < 0 ? 0 : i0 - d + 1); i < i0 + d && i < this->xsize; i++) {
            for (auto k = (k0 - d < 0 ? 0 : k0 - d); k <= k0 + d && k < this->zsize; k++) {
               j = j0 + d;
               if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
                  ret[d].insert(ret[d].begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
               j = j0 - d;
               if (!(j < 0 || j >= this->ysize || this->buckets[i][j][k].empty()))
                  ret[d].insert(ret[d].begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
            }
            for (auto j = (j0 - d + 1 < 0 ? 0 : j0 - d + 1); j < j0 + d && j < this->ysize; j++) {
               k = k0 + d;
               if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
                  ret[d].insert(ret[d].begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
               k = k0 - d;
               if (!(k < 0 || k >= this->zsize || this->buckets[i][j][k].empty()))
                  ret[d].insert(ret[d].begin(), this->buckets[i][j][k].begin(), this->buckets[i][j][k].end());
            }
         }
         tot += ret[d].size();
         if (tot >= limit_number)
            return ret;
      }
      return ret;
   };
   std::vector<std::vector<T>> getObjects(const Tddd &x, const double smoothing_length, const ST limit_number = 100000) const {
      ST depth = std::ceil(smoothing_length / this->dL);
      return getObjects(x, depth, limit_number);
   };
   /* ------------------------------------------------------ */
   std::vector<T> getObjectsFlattened(const Tddd &x, const ST d /*limit depth*/) const {
      std::vector<T> ret;
      Tii i, j, k;
      indices(x, d, i, j, k);
      auto [i_beg, i_end] = i;
      auto [j_beg, j_end] = j;
      auto [k_beg, k_end] = k;
      for (auto it = std::next(this->buckets.begin(), i_beg); it != std::next(this->buckets.begin(), i_end + 1); ++it)
         for (auto jt = std::next(it->begin(), j_beg); jt != std::next(it->begin(), j_end + 1); ++jt)
            for (auto kt = std::next(jt->begin(), k_beg); kt != std::next(jt->begin(), k_end + 1); ++kt)
               ret.insert(ret.end(), kt->begin(), kt->end());
      return ret;
   };
   std::vector<T> getObjectsFlattened(const Tddd &x, const double radius) const {
      std::vector<T> ret;
      Tii i, j, k;
      indices(x, radius, i, j, k);
      auto [i_beg, i_end] = i;
      auto [j_beg, j_end] = j;
      auto [k_beg, k_end] = k;
      for (auto it = std::next(this->buckets.begin(), i_beg); it != std::next(this->buckets.begin(), i_end + 1); ++it)
         for (auto jt = std::next(it->begin(), j_beg); jt != std::next(it->begin(), j_end + 1); ++jt)
            for (auto kt = std::next(jt->begin(), k_beg); kt != std::next(jt->begin(), k_end + 1); ++kt)
               ret.insert(ret.end(), kt->begin(), kt->end());
      return ret;
   };
   /* -------------------------------------------------------------------------- */
   Tdd xrange(const ST i, const double &x) const {
      auto [v0, v1] = (this->dL * i + this->xbounds) - x;  // 注意！v0,v1どちらが最大かわからない
      return {Norm(v0), Norm(v1)};
   };
   Tdd yrange(const ST i, const double &x) const {
      auto [v0, v1] = (this->dL * i + this->ybounds) - x;  // 注意！v0,v1どちらが最大かわからない
      return {Norm(v0), Norm(v1)};
   };
   Tdd zrange(const ST i, const double &x) const {
      auto [v0, v1] = (this->dL * i + this->zbounds) - x;  // 注意！v0,v1どちらが最大かわからない
      return {Norm(v0), Norm(v1)};
   };
   template <typename U>
   bool anyInside(const Tddd &center, const double radius, const U &exceptions) const {
      if (this->buckets.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D buckets is empty");
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_range(center, static_cast<ST>(std::ceil(radius / this->dL)));
      /* ----------------------------- for_eachを使ったループ ---------------------------- */
      bool found = false;
      for_each(std::next(this->buckets.cbegin(), i_min), std::next(this->buckets.cbegin(), i_max), [&](const auto &Bi) {
         if (!found)
            for_each(std::next(Bi.cbegin(), j_min), std::next(Bi.cbegin(), j_max), [&](const auto &Bij) {
               if (!found)
                  for_each(std::next(Bij.cbegin(), k_min), std::next(Bij.cbegin(), k_max), [&](const auto &Bijk) {
                     if (!found)
                        for (const auto &p : Bijk) {
                           if (radius >= Norm(ToX(p) - center) && !MemberQ(exceptions, p)) {
                              found = true;
                              break;
                           }
                        }
                  });
            });
      });
      return found;
   };

   template <typename U>
   bool anyInside(const Tddd &center, const double radius) const {
      if (this->buckets.empty())
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "'s 3D buckets is empty");
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_range(center, static_cast<ST>(std::ceil(radius / this->dL)));
      /* ----------------------------- for_eachを使ったループ ---------------------------- */
      bool found = false;
      for_each(std::next(this->buckets.cbegin(), i_min), std::next(this->buckets.cbegin(), i_max), [&](const auto &Bi) {
         if (!found)
            for_each(std::next(Bi.cbegin(), j_min), std::next(Bi.cbegin(), j_max), [&](const auto &Bij) {
               if (!found)
                  for_each(std::next(Bij.cbegin(), k_min), std::next(Bij.cbegin(), k_max), [&](const auto &Bijk) {
                     if (!found)
                        for (const auto &p : Bijk) {
                           if (radius >= Norm(ToX(p) - center)) {
                              found = true;
                              break;
                           }
                        }
                  });
            });
      });
      return found;
   };
   /* ------------------------------------------------------ */
   //    void apply(const std::function<void(T)> &fun, const Tddd &x, const int d /*limit depth*/, const int limit_number = 100000) {
   //       Tii i, j, k;
   //       indices(x, d, i, j, k);
   //       auto [i_beg, i_end] = i;
   //       auto [j_beg, j_end] = j;
   //       auto [k_beg, k_end] = k;
   //       for (auto it = this->buckets.begin() + i_beg; it != this->buckets.begin() + i_end; ++it)
   //          for (auto jt = it->begin() + j_beg; jt != it->begin() + j_end; ++jt)
   //             for (auto kt = jt->begin() + k_beg; kt != jt->begin() + k_end; ++kt)
   //                for (auto objT = kt->begin(); objT != kt->end(); ++objT)
   //                   fun(*objT);
   //    };
};

template <typename T>
struct Buckets : public BaseBuckets<T> {
   Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<T>(c_bounds, dL_IN){};
   Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<T>(boundingboxIN, dL_IN){};
};

#endif