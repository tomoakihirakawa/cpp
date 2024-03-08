#include <array>
#include "basic_constants.hpp"
#include "lib_multipole_expansion.hpp"

/*DOC_EXTRACT 0_2_1_translation_of_a_multipole_expansion

## `multipole_expansion`クラスのチェック

FMMアルゴリズムでは，展開中心から遠くにある遠方原点の値は，モーメントを計算した後に渡される．
ここでチェックするのは，その計算過程を行うクラス`multipole_expansion`が問題なく動作するかどうかである．

NOTE:境界要素法におけるモーメントは，極そのものではなく，極の面積分（３D）である．

NOTE:多重極モーメントを計算するために，極の値を与えられなければならない．
\cite{Liu_2009}

NOTE:効率化するために要求されるオペレーションは，極の値が変化した際に，できるだけ少ない計算でモーメントを更新することである．

1. モーメントの計算（近傍にある複数の極を変数分離し足し合わせる）
2. 遠方の原点を決めて渡し，計算しておいたモーメントと積和を計算する
3. この計算結果と，展開しない計算結果との差をプロット

一つ前の例では，展開位置を変えることで，多重極展開の精度がどのように変化するかを調べた．
原点位置の移動による展開精度の変化は，展開中心の移動による展開精度の変化と同じである．
展開精度は，（多分）相対距離を規格化した上での，展開中心と極と原点との相対的位置関係で決まっているからである．

## 展開中心の移動（M2M）

多数の極を空間的にグループ分けして，
グループの中心位置を展開中心として多重極展開したとする．

次に，そのグループをさらにまとめて新たな多重極展開を行うことを考える．
この操作は，１ステップ目で得られた各グループの多重極展開係数を利用することで効率的に行うことができる．
各極に対する多重極展開は計算せずに済むからである．

変更されるのは，多重極係数ではなく，球面調和関数自体と，少しの係数のみである．

ここでは，始めに，１ステップ目として座標原点を中心とした多重極展開を行い，
次に，様々な場所での多重極展開を行って，前回同様に精度を検証する．

もし，２ステップ目において，展開中心が１ステップ目同様に原点であれば，
前回と同じ結果が得られるはずである．


```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=test_translation_of_a_multipole_expansion.cpp
make
./test_translation_of_a_multipole_expansion
```

*/

struct FMM_pole {
   double value_FMM = 1.;
   std::array<double, 3> X_FMM = {0., 0., 0.};

   FMM_pole(double value, std::array<double, 3> X) : value_FMM(value), X_FMM(X) {}
};

struct FMM_triangle_face {
   std::array<FMM_pole*, 3> poles = {nullptr, nullptr, nullptr};

   FMM_triangle_face(FMM_pole* p0, FMM_pole* p1, FMM_pole* p2) : poles({p0, p1, p2}) {}

   std::array<std::array<double, 3>, 3> getVerticesPosition() { return {poles[0]->X_FMM, poles[1]->X_FMM, poles[2]->X_FMM}; }
};

template <std::size_t max_n>
struct FMM_multipole_moment {

   std::array<double, 3> center;
   std::vector<FMM_triangle_face*> triangles;
   std::vector<std::tuple<FMM_pole*, std::array<std::array<std::complex<double>, max_n + 1>, max_n + 1>>> pole_values;
   std::array<std::tuple<double, double, double, Tddd>, __array_GW6xGW6__.size()> t0_t1_ww_N012;

   FMM_multipole_moment(const std::vector<FMM_triangle_face*>& triangles, const std::array<double, 3>& center)
       : center(center), triangles(triangles) {
      for (int i = 0; const auto& [t0, t1, ww] : __array_GW6xGW6__)
         t0_t1_ww_N012[i++] = {t0, t1, ww, ModTriShape<3>(t0, t1)};

      // std::array<std::complex<double>, (max_n + 1) * (max_n + 1)> M;
      std::array<std::array<std::complex<double>, max_n + 1>, max_n + 1> M;
      M.fill(0);
      double weight, c;
      int m, n;
      for (const auto& triangle : triangles) {
         auto X012 = triangle->getVerticesPosition();
         for (const auto& pole : triangle->poles) {
            M.fill(0);
            //! make M
            for (const auto& [t0, t1, ww, N01] : this->t0_t1_ww_N012) {

               weight = ww * (1. - t0);
               auto [r_near, a_near, b_near] = ToSphericalCoordinates(Dot(N01, X012) - this->center);

               for (n = 0; n <= max_n; ++n) {
                  c = weight * std::pow(r_near, n);
                  for (m = -n; m <= n; ++m)
                     M[n][max_n + m] += c * Y(n, -m, a_near, b_near);
               }
               
            }
            //! emplace_back M
            this->pole_values.emplace_back(pole, M);
         }
      }
      //! pole_values is ready
   }

   double G_FMM(const std::array<double, 3>& X_far) const {
      int m, n;
      const auto [r_far, a_far, b_far] = ToSphericalCoordinates(X_far - this->center);
      double inv_r_far = 1. / r_far, c;

      //! make S
      std::array<std::complex<double>, (max_n + 1) * (max_n + 1)> S;
      for (n = 0; n <= max_n; ++n)
         for (m = -n; m <= n; ++m)
            S[this->getIndex(n, m)] += inv_r_far * Y(n, m, a_far, b_far);

      //! retrive pole_values
      int i = 0;
      std::complex<double> accum = 0;
      for (const auto& [pole, M] : this->pole_values)
         for (n = 0; n <= max_n; ++n) {
            for (m = -n; m <= n; ++m) {
               i = this->getIndex(n, m);
               accum += (M[i] * inv_r_far) * S[i];
            }
         }
      return accum.real();
   }

   double G_original(std::array<double, 3> X_far) const {
      double accum = 0, inv_r;
      for (const auto& triangle : this->triangles) {
         auto X012 = triangle->getVerticesPosition();
         for (const auto& pole : triangle->poles) {
            for (const auto& [t0, t1, ww, N012] : this->t0_t1_ww_N012)
               accum += ww * (1. - t0) / Norm(Dot(N012, X012) - X_far);
         }
      }
      return accum;
   }

   // const std::complex<double>& getM(const int n, const int m) { return M[(n * n + 1 /*一つ小さいnのサイズ*/) + n + m]; };
   int getIndex(const int n, const int m) const { return (n * n + 1 /*一つ小さいnのサイズ*/) + n + m; };
};

auto p0 = new FMM_pole(1., {0., 0., 0.});
auto p1 = new FMM_pole(1., {1., 0., 0.});
auto p2 = new FMM_pole(1., {0., 1., 0.});
auto p3 = new FMM_pole(1., {0., 0., 1.});

auto face0 = new FMM_triangle_face(p0, p1, p2);
auto face1 = new FMM_triangle_face(p0, p1, p3);

const std::array<double, 3> center = {{0., 0., 0.}};

int main() {

   const double range = 20.0;  // Range of grid from the center
   const double dx = 0.5;      // Grid spacing

   const FMM_multipole_moment<8> multipole_moment({face0, face1}, center);

   // Output file
   std::string name = "./output/multipole_expansion_comparison.dat";
   std::ofstream ofs(name);
   if (!ofs.is_open()) {
      std::cerr << "Failed to open " << name << std::endl;
      return -1;
   }

   std::cout << "X Y Z G_FMM G_Original Error(log10)" << std::endl;

   for (double x = -range; x <= range; x += dx) {
      for (double y = -range; y <= range; y += dx) {
         double z = 0;
         std::array<double, 3> X_far = {x, y, z};
         double G_FMM = multipole_moment.G_FMM(X_far);
         double G_Original = multipole_moment.G_original(X_far);
         double error = std::log10(std::abs(G_FMM - G_Original));

         ofs << x << " " << y << " " << z << " " << G_FMM << " " << G_Original << " " << error << std::endl;
      }
   }

   std::cout << "Data output to " << name << std::endl;
}