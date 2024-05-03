#include <chrono>
#include <iostream>
#include <tuple>
#include <vector>

/*

```cpp
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=DataStructurePerformanceComparison.cpp
make
./DataStructurePerformanceComparison
```

*/

struct MyStruct {
   std::size_t i;
   int x;
   double y;
};

int main() {
   const int N = 10000000;
   std::vector<MyStruct> vecStruct;
   std::vector<std::tuple<std::size_t, int, double>> vecTuple;
   std::vector<std::size_t> vecI(N);
   std::vector<int> vecX(N);
   std::vector<double> vecY(N);

   // データの初期化
   for (int i = 0; i < N; i++) {
      vecStruct.push_back(MyStruct{i, i, i * 1.1});
      vecTuple.push_back(std::make_tuple(i, i, i * 1.1));
      vecI[i] = i;
      vecX[i] = i;
      vecY[i] = i * 1.1;
   }

   std::size_t sum_i = 0;
   int sum_x = 0;
   double sum_y = 0;

   auto start = std::chrono::high_resolution_clock::now();
   // MyStruct のループ
   for (const auto& item : vecStruct) {
      sum_i += item.i;
      sum_x += item.x;
      sum_y += item.y;
   }
   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = end - start;
   std::cout << "MyStruct Loop Time: " << elapsed.count() << " seconds\n";

   start = std::chrono::high_resolution_clock::now();
   // tuple のループ
   for (const auto& item : vecTuple) {
      sum_i += std::get<0>(item);
      sum_x += std::get<1>(item);
      sum_y += std::get<2>(item);
   }
   end = std::chrono::high_resolution_clock::now();
   elapsed = end - start;
   std::cout << "Tuple Loop Time: " << elapsed.count() << " seconds\n";

   start = std::chrono::high_resolution_clock::now();
   // SoA のループ
   for (int i = 0; i < N; i++) {
      sum_i += vecI[i];
      sum_x += vecX[i];
      sum_y += vecY[i];
   }
   end = std::chrono::high_resolution_clock::now();
   elapsed = end - start;
   std::cout << "SoA Loop Time: " << elapsed.count() << " seconds\n";

   // 結果を出力して最適化を防ぐ
   std::cout << "Sums: " << sum_i << ", " << sum_x << ", " << sum_y << "\n";

   return 0;
}
