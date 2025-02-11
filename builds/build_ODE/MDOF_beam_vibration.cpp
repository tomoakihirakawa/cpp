#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "basic_linear_systems.hpp"
#include "integrationOfODE.hpp"

/*DOC_EXTRACT beam_vibration

# 梁の振動

水平面で固定され垂直に立てられた梁の振動を考える．
梁の上端部にステップ荷重を加えた後，梁はどのように振動するだろうか．

梁の運動方程式は，次のように表される．慣性力と曲げ応力と圧縮応力の和がゼロとなるという原理に基づいている．
慣性力を表す１階時間偏微分と，曲げ応力と圧縮応力を表す４階空間偏微分を含む偏微分方程式である．

```math
\begin{align*}
\frac{\partial^2 u}{\partial t^2} &= \frac{E I}{\rho A} \frac{\partial^4 u}{\partial x^4} - \frac{E I}{\rho A} \frac{\partial^2}{\partial x^2} \left( \frac{\partial^2 u}{\partial x^2} \right) + \frac{E I}{\rho A} \frac{\partial^2}{\partial x^2} \left( \frac{\partial^2 u}{\partial x^2} \right)\\
&= \frac{E I}{\rho A} \frac{\partial^4 u}{\partial x^4}
\end{align*}
```

$`L = T - U`$をラグランジアンとすると，ラグランジアンのオイラー・ラグランジュ方程式は次のようになる．

```math
\begin{align*}
\frac{d}{dt} \left( \frac{\partial L}{\partial \dot{u}} \right) - \frac{\partial L}{\partial u} &= 0\\
\frac{d}{dt} \left( \frac{\partial T}{\partial \dot{u}} \right) - \frac{\partial T}{\partial u} &= 0\\
\frac{d}{dt} \left( \frac{\partial}{\partial \dot{u}} \left( \frac{1}{2} \rho A \dot{u}^2 \right) \right) - \frac{\partial}{\partial u} \left( \frac{1}{2} \rho A \dot{u}^2 \right) &= 0\\
\rho A \ddot{u} &= \frac{\partial}{\partial x} \left( E I \frac{\partial^2 u}{\partial x^2} \right)
\end{align*}
```

これを離散化し，モード展開法と，動的解析法を用いた方法を行なってみよう．
そして，最終的に，ステップ荷重をケーブル荷重に置き換えた場合，どのような違いが生じるかを調べてみよう．

---

## **1. 元の偏微分方程式 (PDE)**
梁の振動を記述する式は、空間的な**４階微分項**と時間的な**２階微分項**を含む偏微分方程式です。これは、**オイラー・ベルヌーイの梁理論**に基づいて次のように書かれます：

\[
\rho A \frac{\partial^2 u}{\partial t^2} = -E I \frac{\partial^4 u}{\partial x^4} + f(x, t)
\]

### **項の意味**:
1. \( \rho A \frac{\partial^2 u}{\partial t^2} \): 梁の単位長さあたりの質量に基づく慣性力（時間方向の加速度に比例）。
2. \( -E I \frac{\partial^4 u}{\partial x^4} \): 曲げ剛性による復元力（空間方向の４階微分）。
3. \( f(x, t) \): 外力（位置 \( x \) と時間 \( t \) の関数）。

この式は**連続体モデル**を基にしており、無限の自由度を持つモデルとして梁の挙動を完全に記述します。

## **4. PDEとODEの関係をまとめる**
1. **PDEの連続性**:
   - \( \rho A \frac{\partial^2 u}{\partial t^2} = -E I \frac{\partial^4 u}{\partial x^4} + f(x, t) \) は無限の自由度を持つ連続体のモデル。

2. **空間離散化**:
   - 有限要素法やモード解析により、空間の自由度を有限に縮約。
   - ４階微分項が剛性マトリクス \( K \) に、質量項が質量マトリクス \( M \) に対応。

3. **ODEによる記述**:
   - PDEが節点ごとの運動方程式 \( M \ddot{u} + K u = F(t) \) に変換される。

4. **モード解析**:
   - ODEをモード座標に変換し、固有振動数と振動形状で簡略化して解く。

## モード展開法と動的解析法の比較

動的解析とモード展開法の結果をフーリエ空間で比較してみると，
モード展開で得られた固有振動数は，動的解析法で得られた固有振動数と一致することがわかる．
当然，分割数に伴うズレも含めて両者の結果が一致することが確認できる．

<img src="sample_modal_analysis.png" width="600px">

解いている離散化された運動方程式が同じであるから当然の結果である．

比較とは別に，節点の増加に伴う階の収束性を知ることができる．
SDOFの結果は，節点数が多いMDOFと大きく異なる．

## ステップ荷重は，ケーブル荷重の代わりになるか？

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=MDOF_beam_vibration.cpp
make
./MDOF_beam_vibration <numNodes>
```

*/

// 質量マトリクスと剛性マトリクスの生成
struct BeamVibration {

   std::vector<std::vector<double>> massMatrix;
   std::vector<std::vector<double>> stiffnessMatrix;
   std::vector<std::vector<double>> dampingMatrix;
   std::vector<double> freq;
   std::vector<std::vector<double>> mode_shape;

   std::vector<std::vector<double>> M, K, C;

   const double E = 205 * 1e9;  // ヤング率 205GPa
   // double h = 0.02;
   // double b = 0.035;
   const double b = 0.01;
   const double h = 0.005;
   const double rho = 7.668 / 1e3 / std::pow(0.01, 3);  // 密度 7.8g/cm^3
   const double A = b * h;                              // 断面積
   const double L = 0.9;                                // 長さ
   double dx;
   const double I = b * std::pow(h, 3) / 12;  // 断面二次モーメント
   const double k = E * I / L / L / L;
   double zeta = 0.;  // 減衰比
   double mass_of_a_node;
   double c;
   // double c = 2.0 * zeta * std::sqrt(mass_of_a_node * k);

   BeamVibration(int numNodes) : dx(L / numNodes) { update(numNodes); }

   // オイラーベルヌーイ梁の質量マトリクス，剛性マトリクス，減衰マトリクスの更新
   // ただし，固定の始めに二点は式に含めず常に変異ゼロとする
   void update(int numNodes) {
      double EI_dx4 = E * I / pow(dx, 4);
      double total_mass = rho * A * L;
      this->mass_of_a_node = total_mass / numNodes;
      c = 2.0 * zeta * std::sqrt(mass_of_a_node * k);

      massMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
      stiffnessMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
      dampingMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
      massMatrix.assign(numNodes, std::vector<double>(numNodes, 0.0));
      stiffnessMatrix.assign(numNodes, std::vector<double>(numNodes, 0.0));
      dampingMatrix.assign(numNodes, std::vector<double>(numNodes, 0.0));

      // 剛性マトリクスの内部ノード (4階中心差分)
      for (int i = 0; i < numNodes - 2; ++i) {
         if (i == 0) {
            // do nothing
         } else if (i == 1) {
            stiffnessMatrix[i][i - 1] = -4 * EI_dx4;
         } else {
            stiffnessMatrix[i][i - 2] = EI_dx4;
            stiffnessMatrix[i][i - 1] = -4 * EI_dx4;
         }
         stiffnessMatrix[i][i] = 6 * EI_dx4;
         stiffnessMatrix[i][i + 1] = -4 * EI_dx4;
         stiffnessMatrix[i][i + 2] = EI_dx4;
      }

      for (int i = 0; i < numNodes; ++i)
         massMatrix[i][i] = rho * A;

      bool based_on_the_reference = false;

      if (based_on_the_reference) {
         // 自由端（せん断力ゼロ）
         stiffnessMatrix[numNodes - 2][numNodes - 4] = EI_dx4;
         stiffnessMatrix[numNodes - 2][numNodes - 3] = -3 * EI_dx4;
         stiffnessMatrix[numNodes - 2][numNodes - 2] = 3 * EI_dx4;
         stiffnessMatrix[numNodes - 2][numNodes - 1] = -EI_dx4;

         // 自由端（曲げモーメントゼロ）
         stiffnessMatrix[numNodes - 1][numNodes - 3] = EI_dx4 / 2;
         stiffnessMatrix[numNodes - 1][numNodes - 2] = -2 * EI_dx4 / 2;
         stiffnessMatrix[numNodes - 1][numNodes - 1] = EI_dx4 / 2;

         // 2 (wi - 2 wm1 + wm2)
         massMatrix[numNodes - 1][numNodes - 1] = rho * A;
         massMatrix[numNodes - 2][numNodes - 2] = rho * A;
         massMatrix[numNodes - 2][numNodes - 1] = rho * A;
      } else {
         stiffnessMatrix[numNodes - 2][numNodes - 4] = EI_dx4;
         stiffnessMatrix[numNodes - 2][numNodes - 3] = -4 * EI_dx4;
         stiffnessMatrix[numNodes - 2][numNodes - 2] = 5 * EI_dx4;
         stiffnessMatrix[numNodes - 2][numNodes - 1] = -2 * EI_dx4;

         stiffnessMatrix[numNodes - 1][numNodes - 3] = 2 * EI_dx4;
         stiffnessMatrix[numNodes - 1][numNodes - 2] = -4 * EI_dx4;
         stiffnessMatrix[numNodes - 1][numNodes - 1] = 2 * EI_dx4;
      }

      massMatrix[numNodes - 1][numNodes - 1] /= 2;

      auto eigen = Eigensystem(stiffnessMatrix, massMatrix);
      this->freq = eigen.first;
      this->mode_shape = eigen.second;
      auto tr = Transpose(this->mode_shape);
      for (auto& shape : tr)
         shape = Normalize(shape);  // モード毎に規格化する
      this->mode_shape = Transpose(tr);

      // ** 減衰マトリクス（単純な対角行列）**
      // ** レイリー減衰行列の適用 **
      double alpha = 0.05;  // 質量比例減衰
      double beta = 0.001;  // 剛性比例減衰

      // double alpha = 0.;  // 質量比例減衰
      // double beta = 0.;   // 剛性比例減衰

      for (int i = 0; i < numNodes; ++i) {
         for (int j = 0; j < numNodes; ++j) {
            dampingMatrix[i][j] = alpha * massMatrix[i][j] + beta * stiffnessMatrix[i][j];
         }
      }

      for (auto& m : massMatrix)
         for (auto& mm : m)
            mm *= dx;

      for (auto& s : stiffnessMatrix)
         for (auto& ss : s)
            ss *= dx;

      for (auto& d : dampingMatrix)
         for (auto& dd : d)
            dd *= dx;

      /*
       モード座標系の質量マトリクス，剛性マトリクス，減衰マトリクスを計算
      */
      auto mode_shape_tr = Transpose(this->mode_shape);
      this->M = Dot(mode_shape_tr, Dot(this->massMatrix, this->mode_shape));
      this->K = Dot(mode_shape_tr, Dot(this->stiffnessMatrix, this->mode_shape));
      this->C = Dot(mode_shape_tr, Dot(this->dampingMatrix, this->mode_shape));
   }
};

// モーダル座標系の加速度計算
std::vector<double> modalAcceleration(const std::vector<double>& q,
                                      const std::vector<double>& q_dot,
                                      const BeamVibration& beam,
                                      const std::vector<double>& f_q /*[N]*/) {
   std::vector<double> ans(q.size(), 0.0), accel(q.size(), 0.0);
   lapack_lu lu(beam.M);
   lu.solve(Dot(beam.C, q_dot), ans);
   accel -= ans;
   lu.solve(Dot(beam.K, q), ans);
   accel -= ans;
   auto _f = Dot(Transpose(beam.mode_shape), f_q);
   lu.solve(_f, ans);
   accel += ans;
   return accel;
}

// 加速度の計算
std::vector<double> acceleration(const std::vector<double>& x,
                                 const std::vector<double>& v,
                                 const BeamVibration& beam,
                                 const std::vector<double>& externalForceVector /*[N]*/) {
   std::vector<double> ans(x.size(), 0.0), accel(x.size(), 0.0);
   lapack_lu lu(beam.massMatrix);
   lu.solve(Dot(beam.dampingMatrix, v), ans);
   accel -= ans;
   lu.solve(Dot(beam.stiffnessMatrix, x), ans);
   accel -= ans;
   lu.solve(externalForceVector, ans);
   accel += ans;
   return accel;
}

// get numNode fron argument
int main(int argc, char* argv[]) {

   if (argc < 3) {
      std::cerr << "Usage: " << argv[0] << " <numNodes>" << " <index the force is applied>" << std::endl;
      return 1;
   }

   const int numNodes = std::stoi(argv[1]);
   const int index_force_exerted = std::stoi(argv[2]);

   // int numNodes = 10;  // 節点の数

   const double end_time = 5.0;  // 終了時間
   double dt = 0.00001;          // 時間ステップ
   const double t0 = 0.0;        // 初期時間

   std::vector<std::vector<double>> massMatrix, stiffnessMatrix, dampingMatrix;

   auto beam = BeamVibration(numNodes);
   // beam.zeta = 1.;
   beam.update(numNodes);

   std::cout << Red;
   std::cout << "質量マトリクス" << std::endl;
   std::cout << beam.massMatrix << std::endl;
   std::cout << Magenta;
   std::cout << "剛性マトリクス" << std::endl;
   std::cout << beam.stiffnessMatrix << std::endl;
   std::cout << Yellow;
   std::cout << "減衰マトリクス" << std::endl;
   std::cout << beam.dampingMatrix << std::endl;
   std::cout << Blue;
   std::cout << "時間依存部分，振動数" << std::endl;
   std::cout << beam.freq << std::endl;
   std::cout << Cyan << "空間依存部分，モード形状" << std::endl;
   std::cout << beam.mode_shape << std::endl;
   std::cout << colorReset;
   std::cout << "質量マトリクス in モード座標系" << std::endl;
   std::cout << beam.M << std::endl;
   std::cout << "剛性マトリクス in モード座標系" << std::endl;
   std::cout << beam.K << std::endl;
   std::cout << "減衰マトリクス in モード座標系" << std::endl;
   std::cout << beam.C << std::endl;

   const std::vector<double> initialDisplacement(numNodes, 0), initialVelocity(numNodes, 0);

   /* -------------------------------------------------------------------------- */

   //% 外力の定義（例：節点2に対して時間tに依存する外力を与える）

   auto externalForce = [&numNodes, &index_force_exerted](int node_index, double t) -> double {
      // if (node_index == numNodes - 2 && 0. < t)
      //    return 1;  // 10.0 * sin(t);  // 例として正弦波の外力
      if (node_index == index_force_exerted)
         return 1;  // 10.0 * sin(t);  // 例として正弦波の外力

      // if (node_index == int((double)numNodes / 2.) - 1 && 0.5 < t)
      //    return 0.1;  // 10.0 * sin(t);  // 例として正弦波の外力
      return 0.0;
   };

   auto externalForceVector = [&](double t) {
      std::vector<double> f(numNodes, 0.0);
      for (size_t i = 0; i < numNodes; ++i)
         f[i] = externalForce(i, t);
      return f;
   };

   /* -------------------------------------------------------------------------- */

   // Runge-Kutta法による解の計算
   std::vector<double> time, time_modal;
   std::vector<std::vector<double>> displacement, displacement_modal;
   std::vector<std::vector<double>> velocity, velocity_modal;
   double t = t0;
   std::vector<double> x = initialDisplacement;
   std::vector<double> v = initialVelocity;

   double output_start_time = t0, output_dt = 0.001;
   x[0] = 0;
   v[0] = 0;
   bool updated = false;
   std::vector<double> accel;

   std::cout << Red << "実際の空間座標系での動的解析 ... numNodes = " << numNodes << colorReset << std::endl;
   /* ---------------------------------- 実際の空間座標系での動的解析 ---------------------------------- */
   double dt_tmp = dt;
   for (auto j = 0; j < 1000000000; ++j) {
      if (j < 100)
         dt = 1E-7;
      else
         dt = dt_tmp;
      RungeKutta<std::vector<double>> RK_x(dt, t, x, 4), RK_v(dt, t, v, 4);
      do {
         x = RK_x.get_x();
         v = RK_v.get_x();
         RK_x.push(v);
         RK_v.push(accel = acceleration(x, v, beam, externalForceVector(t)));
      } while (!RK_x.finished);
      t = RK_x.get_t();
      x = RK_x.get_x();
      v = RK_v.get_x();
      if (output_start_time < t) {
         if (isFinite(x[0]) == false)
            break;
         std::vector<double> row_x = {}, row_v = {};
         row_x.insert(row_x.end(), x.begin(), x.end());
         row_v.insert(row_v.end(), v.begin(), v.end());
         output_start_time += output_dt;
         //
         time.emplace_back(t);
         displacement.emplace_back(row_x);
         velocity.emplace_back(row_v);
      }

      if (end_time < t)
         break;
   }

   /* ---------------------------------- モード座標系での動的解析 ---------------------------------- */
   t = t0;
   output_start_time = t0;

   // 時間積分のための変数
   std::vector<double> q = Dot(Transpose(beam.mode_shape), initialDisplacement);
   std::vector<double> q_dot = Dot(Transpose(beam.mode_shape), initialVelocity);
   std::cout << Green << "モード座標系での動的解析 ... numNodes = " << numNodes << colorReset << std::endl;
   for (auto j = 0; j < 1000000000; ++j) {
      RungeKutta<std::vector<double>> RK_q(dt, t, q, 4), RK_q_dot(dt, t, q_dot, 4);
      do {
         q = RK_q.get_x();
         q_dot = RK_q_dot.get_x();
         RK_q.push(q_dot);
         RK_q_dot.push(accel = modalAcceleration(q, q_dot, beam, externalForceVector(t)));
      } while (!RK_q.finished);
      t = RK_q.get_t();
      q = RK_q.get_x();
      q_dot = RK_q_dot.get_x();
      if (output_start_time < t) {
         if (isFinite(x[0]) == false)
            break;
         std::vector<double> row_q = {}, row_q_dot = {}, row_x = {};
         row_q.insert(row_q.end(), q.begin(), q.end());
         std::vector<double> X = Dot(beam.mode_shape, q);
         row_x.insert(row_x.end(), X.begin(), X.end());
         row_q_dot.insert(row_q_dot.end(), q_dot.begin(), q_dot.end());
         output_start_time += output_dt;
         //
         time_modal.emplace_back(t);
         displacement_modal.emplace_back(row_q);
         velocity_modal.emplace_back(row_q_dot);
      }
      if (end_time < t)
         break;
   }

   // 結果の保存

   JSONoutput json;

   json.set("time", time);
   json.set("displacement", displacement);
   json.set("velocity", velocity);
   json.set("index_force_exerted", std::vector<int>{index_force_exerted});
   json.set("x_force_exerted", std::vector<double>{(index_force_exerted + 1.) * beam.L / numNodes});
   json.set("time_modal", time_modal);
   json.set("displacement_modal", displacement_modal);
   json.set("velocity_modal", velocity_modal);
   json.set("mode_shape", beam.mode_shape);  // モード形状を追加
   json.set("freq", beam.freq);              // 振動数を追加
   json.set("M", beam.massMatrix);           // モード座標系の質量マトリクスを追加
   json.set("K", beam.stiffnessMatrix);      // モード座標系の剛性マトリクスを追加

   std::string output_filename = "./output/result_" + std::to_string(numNodes) + "_" + std::to_string(index_force_exerted) + ".json";
   std::ofstream os(output_filename);
   json.output(os);
   os.close();

   std::cout << "結果が" + output_filename + "に保存されました。" << std::endl;

   return 0;
}