#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "integrationOfODE.hpp"

// 質量マトリクスと剛性マトリクスの生成
void generateMatrices(int numNodes, std::vector<std::vector<double>>& massMatrix, std::vector<std::vector<double>>& stiffnessMatrix, std::vector<std::vector<double>>& dampingMatrix) {
   massMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
   stiffnessMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));
   dampingMatrix.resize(numNodes, std::vector<double>(numNodes, 0.0));

   double E = 200 * 1e8;      // ヤング率
   double A = 0.035 * 0.002;  // 断面積
   double L = 0.9;            // 長さ
   double k = E * A / L;

   for (int i = 0; i < numNodes; ++i) {
      massMatrix[i][i] = 1.0;     // 単純な質量マトリクス
      dampingMatrix[i][i] = 1.0;  // 単純な減衰マトリクス
      if (i > 0) {
         stiffnessMatrix[i][i] = 2.0 * k;
         stiffnessMatrix[i][i - 1] = -1.0 * k;
         stiffnessMatrix[i - 1][i] = -1.0 * k;
      }
   }
   stiffnessMatrix[0][1] = -1.0 * k;  // 固定端の剛性
   stiffnessMatrix[0][0] = 100. * k;  // 固定端の剛性
}

// 初期条件
void initializeConditions(int numNodes, std::vector<double>& initialDisplacement, std::vector<double>& initialVelocity) {
   initialDisplacement.resize(numNodes, 0.0);
   initialVelocity.resize(numNodes, 0.0);
   //    initialDisplacement[numNodes - 1] = 0.1;  // 初期変位
}

// 加速度の計算
std::vector<double> acceleration(const std::vector<double>& x, const std::vector<double>& v, const std::vector<std::vector<double>>& massMatrix, const std::vector<std::vector<double>>& stiffnessMatrix, const std::vector<std::vector<double>>& dampingMatrix, const std::function<double(int, double)>& externalForce, double t) {
   std::vector<double> a(x.size(), 0.0);
   for (size_t i = 0; i < x.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
         a[i] += -(dampingMatrix[i][j] * v[j] + stiffnessMatrix[i][j] * x[j]) / massMatrix[i][i];
      }
      a[i] += externalForce(i, t) / massMatrix[i][i];  // 外力を加える
   }
   return a;
}

int main() {
   int numNodes = 10;       // 節点の数
   double end_time = 20.0;  // 終了時間
   double dt = 0.000001;    // 時間ステップ
   double t0 = 0.0;         // 初期時間

   std::vector<std::vector<double>> massMatrix, stiffnessMatrix, dampingMatrix;
   generateMatrices(numNodes, massMatrix, stiffnessMatrix, dampingMatrix);

   std::vector<double> initialDisplacement, initialVelocity;
   initializeConditions(numNodes, initialDisplacement, initialVelocity);

   // 外力の定義（例：節点2に対して時間tに依存する外力を与える）
   auto externalForce = [&numNodes](int node, double t) -> double {
      if (node == numNodes - 1 && 1 < t && t < 10)
         return 50;  // 10.0 * sin(t);  // 例として正弦波の外力
      return 0.0;
   };

   // Runge-Kutta法による解の計算
   std::vector<std::vector<double>> RK_t_x({{t0}});
   RK_t_x[0].insert(RK_t_x[0].end(), initialDisplacement.begin(), initialDisplacement.end());
   std::vector<std::vector<double>> RK_t_v({{t0}});
   RK_t_v[0].insert(RK_t_v[0].end(), initialVelocity.begin(), initialVelocity.end());
   double t = t0;
   std::vector<double> x = initialDisplacement;
   std::vector<double> v = initialVelocity;
   double t_to_output = 0.0;
   for (auto j = 0; j < 1000000000; ++j) {
      x[0] = 0;
      v[0] = 0;
      RungeKutta<std::vector<double>> RK_x(dt, t, x, 4);
      RungeKutta<std::vector<double>> RK_v(dt, t, v, 4);
      do {
         x = RK_x.get_x();
         v = RK_v.get_x();
         RK_x.push(v);
         RK_v.push(acceleration(x, v, massMatrix, stiffnessMatrix, dampingMatrix, externalForce, t));
      } while (!RK_x.finished);
      t = RK_x.get_t();
      x = RK_x.get_x();
      v = RK_v.get_x();
      std::vector<double> row_x = {t};
      std::vector<double> row_v = {t};
      if (t_to_output < t) {
         row_x.insert(row_x.end(), x.begin(), x.end());
         RK_t_x.push_back(row_x);
         row_v.insert(row_v.end(), v.begin(), v.end());
         RK_t_v.push_back(row_v);
         t_to_output += 0.1;
      }
      if (end_time < t) {
         break;
      }
   }
   // 結果の保存
   auto result = RK_t_x;
   std::ofstream outFile("result.csv");
   outFile << "Time";
   for (int i = 0; i < numNodes; ++i) {
      outFile << ",Displacement" << i + 1;
   }
   outFile << std::endl;
   for (const auto& row : result) {
      for (size_t i = 0; i < row.size(); ++i) {
         outFile << row[i];
         if (i < row.size() - 1) {
            outFile << ",";
         }
      }
      outFile << std::endl;
   }
   outFile.close();

   std::cout << "結果がresult.csvに保存されました。" << std::endl;

   return 0;
}