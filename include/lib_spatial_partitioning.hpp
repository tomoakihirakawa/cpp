#ifndef lib_spatial_partitioning_H
#define lib_spatial_partitioning_H

template <typename T>
struct BaseBuckets {
   using sizeType = int;
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
   std::vector<std::vector<std::vector<std::unordered_set<T>>>> buckets;
   std::unordered_set<T> all_stored_objects;
   std::unordered_map<T, ST3> map_to_ijk;
   double dL;
   double bucketVolume() const { return std::pow(this->dL, 3.); };
   BaseBuckets(const CoordinateBounds &c_bounds, const double dL_IN) : bounds(c_bounds.bounds), dL(dL_IN), dn(ST3{0, 0, 0}) { initialize(this->bounds, dL_IN); };
   BaseBuckets(const T3Tdd &boundingboxIN, const double dL_IN) : bounds(boundingboxIN), dL(dL_IN), dn(ST3{0, 0, 0}) { initialize(this->bounds, dL_IN); };
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
   Tddd itox(const ST i, const ST j, const ST k) const {
      return {this->dL * 0.5 + this->dL * i + std::get<0>(this->xbounds),
              this->dL * 0.5 + this->dL * j + std::get<0>(this->ybounds),
              this->dL * 0.5 + this->dL * k + std::get<0>(this->zbounds)};
   };
   Tddd itox(const ST3 &ijk) const { return itox(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk)); };
   ST3 indices_no_clamp(const Tddd &x) const {
      return {static_cast<ST>((std::get<0>(x) - std::get<0>(this->xbounds)) / this->dL),
              static_cast<ST>((std::get<1>(x) - std::get<0>(this->ybounds)) / this->dL),
              static_cast<ST>((std::get<2>(x) - std::get<0>(this->zbounds)) / this->dL)};
   };
   ST3 indices(const Tddd &x) const {
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
   ST6 indices_ranges(const Tddd &x, const double d) const {
      auto [i_min, j_min, k_min] = indices(x - d);
      auto [i_max, j_max, k_max] = indices(x + d);
      return {i_min, i_max, j_min, j_max, k_min, k_max};
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
         this->buckets[i][j][k].erase(p);
         this->map_to_ijk.erase(it);
      }
   };
   //@ -------------------------------------------------------------------------- */
   //! explicitly specify coordinates
   bool add_force(const ST3 &ijk, const T p) {
      const auto [i, j, k] = ijk;
      this->map_to_ijk[p] = ijk;
      const bool bucket_inserted = this->buckets[i][j][k].emplace(p).second;
      const bool all_objects_inserted = this->all_stored_objects.emplace(p).second;
      return bucket_inserted && all_objects_inserted;
   };
   bool add(const Tddd &x, const T p) { return add_force(indices(x), p); };
   //! automatically specify coordinates
   bool add(const std::unordered_set<T> &P) {
      bool ret = true;
      auto add_ijk = [&]() {
         for (const auto &p : P)
            this->map_to_ijk[p] = indices(ToX(p));
      };
      auto add_object = [&]() {
               for (const auto &p : P) {
                  ret = ret && add_force(indices(ToX(p)),p);
            } };
#pragma omp parallel sections
      {
#pragma omp section
         add_ijk();
#pragma omp section
         add_object();
      };
      return ret;
   };
   /* -------------------------------------------------------------------------- */
   auto getBucket(const Tddd &x) const {
      const auto ijk = this->indices(x);
      return this->buckets[std::get<0>(ijk)][std::get<1>(ijk)][std::get<2>(ijk)];
   };
   void clear() {
      this->buckets.clear();
      this->all_stored_objects.clear();
   };
   // b@ -------------------------------------------------------------------------- */
   // b@                             STL like functions                             */
   // b@ -------------------------------------------------------------------------- */
   //! count_if
   ST count_if(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      ST ret = 0;
      if (!this->buckets.empty()) {
         const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
         for (auto it = std::next(this->buckets.begin(), i_min); it != std::next(this->buckets.begin(), i_max); ++it)
            for (auto jt = std::next(it->begin(), j_min); jt != std::next(it->begin(), j_max); ++jt)
               for (auto kt = std::next(jt->begin(), k_min); kt != std::next(jt->begin(), k_max); ++kt)
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
         const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
         for (auto it = std::next(this->buckets.begin(), i_min); it != std::next(this->buckets.begin(), i_max); ++it)
            for (auto jt = std::next(it->begin(), j_min); jt != std::next(it->begin(), j_max); ++jt)
               for (auto kt = std::next(jt->begin(), k_min); kt != std::next(jt->begin(), k_max); ++kt)
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
         const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
         for (auto it = std::next(this->buckets.begin(), i_min); it != std::next(this->buckets.begin(), i_max); ++it)
            for (auto jt = std::next(it->begin(), j_min); jt != std::next(it->begin(), j_max); ++jt)
               for (auto kt = std::next(jt->begin(), k_min); kt != std::next(jt->begin(), k_max); ++kt)
                  for (const auto &p : *kt)
                     if (!func(p))
                        return false;
      }
      return true;
   };
   //! any_of
   bool any_of(const Tddd &x, const double d, const std::function<bool(const T &)> &func) const {
      if (!this->buckets.empty()) {
         const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, d);
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
      const auto [i_min, i_max, j_min, j_max, k_min, k_max] = indices_ranges(x, 2. * d);
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
};

template <typename T>
struct Buckets : public BaseBuckets<T> {
   Buckets(const CoordinateBounds &c_bounds, const double dL_IN) : BaseBuckets<T>(c_bounds, dL_IN){};
   Buckets(const T3Tdd &boundingboxIN, const double dL_IN) : BaseBuckets<T>(boundingboxIN, dL_IN){};
};

#endif