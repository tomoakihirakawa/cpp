
#include "Network.hpp"
#include "NetworkUtility.hpp"

/*DOC_EXTRACT REMESH

# Fusion360を使って計算用objファイルを生成

| step1 | step2 | step3 | step4 |
|:-----:|:-----:|:-----:|:-----:|
|![./sample_fusion360_step1.png](./sample_fusion360_step1.png)|![./sample_fusion360_step2.png](./sample_fusion360_step2.png)|![./sample_fusion360_step3.png](./sample_fusion360_step3.png)|![./sample_fusion360_step4.png](./sample_fusion360_step4.png)|

# 計算用にメッシュの細分化

## メッシュの細分化の方法

\insert{flip}

\insert{smoothing_vector}

## 実行ファイルの作成方法（`remesh.cpp`のコンパイル方法）

1. `sh clean`で古いファイルを削除する．
2. `cmake`を使って，`CMakeLists.txt`から`Makefile`を生成する．（リリースビルドタイプを指定し，ソースファイルを`remesh.cpp`に設定）
3. `make`コマンドで，`Makefile`に基づいて実行ファイルをコンパイルする.

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp
make
```

`remesh`という実行ファイルができる．`n`回の細分化を行う場合は，以下のように実行する．

```
./remesh input_file output_dir output_name n
```

![./sample2.gif](./sample2.gif)

出力は，`output_dir/output_name*.vtu`と`output_dir/output_name*.obj`．


![./sample.gif](./sample.gif)

次は，入力ファイルを生成し，計算をする．

*/

const std::array<std::string, 4> icons = {Red + '|' + colorReset, Magenta + '/' + colorReset, Blue + '-' + colorReset, Green + '\\' + colorReset};
const double small = M_PI / 180 * 0.01;

JSONoutput output_json;

int main(int arg, char **argv) {
   Timer time;
   std::string input_file{argv[1]};        // input
   std::string output_directory{argv[2]};  // output dir
   std::string output_name{argv[3]};       // output name
   int remesh = std::atoi(argv[4]);
   Network net(input_file, output_name);
   mk_vtu(output_directory + "/" + output_name + ".vtu", {net.getFaces()});
   net.displayStates();

   PVDWriter PVD(output_directory + output_name + ".pvd");

   V_netLp lines;
   Histogram Histo;

   uomap_P_d P_isFlat;
   auto getData = [&]() -> VV_VarForOutput {
      P_isFlat.clear();
      for (const auto &p : net.getPoints())
         P_isFlat[p] = isFlat(p);
      return {{"isFlat", P_isFlat}};
   };

   V_d hist_bin_width;
   V_d hist_where_divide;
   std::vector<V_i> hist_count;  // 角瓶の中にあるデータの数
   std::vector<V_d> hist_cumulative;
   std::vector<V_i> hist_diff;
   std::vector<V_i> hist_diffdiff;
   std::vector<V_d> hist_interval;
   std::vector<V_d> hist_mid_interval;

   int output_count = 0;
   for (auto count = 0; count <= remesh; ++count) {

      std::cout << "\r" << icons[(count + 1) % icons.size()] << " time:" << time() << std::flush;
      int index = 0;
      /* ------------------------*/
      //* | 0 | 1 | 2 |[3]| 4 |     diff, mid_interval
      //!   0   1   2  [3]  4   5   interval, cumulative_count
      /* ------------------------*/
      Histo.set(extLength(net.getLines()));
      auto it_cumulative = std::ranges::find_if(Histo.cumulative, [](double value) { return value >= 0.9; });
      if (it_cumulative != Histo.cumulative.end())
         index = std::distance(Histo.cumulative.begin(), it_cumulative);
      auto it_count = std::max_element(Histo.count.begin() + index, Histo.count.end());
      if (it_count != Histo.count.end())
         index = std::distance(Histo.count.begin(), it_count);

      /* -------------------------------------------------- */
      // Output to JSON file
      hist_where_divide.push_back(Histo.interval[index]);
      hist_bin_width.push_back(Histo.bin_width);
      hist_count.push_back(Histo.count);
      hist_cumulative.push_back(Histo.cumulative);
      hist_diff.push_back(Histo.diff);
      hist_diffdiff.push_back(Histo.diffdiff);
      hist_interval.push_back(Histo.interval);
      hist_mid_interval.push_back(Histo.mid_interval);
      std::ofstream file("histogram.json");
      file << "{";
      std::array<std::string, 2> bracket = {"[", "]"};
      file << "\"where_divide\" : " << std::make_tuple(hist_where_divide, bracket) << ", ";
      file << "\"bin_width\" : " << std::make_tuple(hist_bin_width, bracket) << ", ";
      file << "\"count\" : " << std::make_tuple(hist_count, bracket) << ", ";
      file << "\"cumulative\" : " << std::make_tuple(hist_cumulative, bracket) << ", ";
      file << "\"diff\" : " << std::make_tuple(hist_diff, bracket) << ", ";
      file << "\"diffdiff\" : " << std::make_tuple(hist_diffdiff, bracket) << ", ";
      file << "\"interval\" : " << std::make_tuple(hist_interval, bracket) << ", ";
      file << "\"mid_interval\" : " << std::make_tuple(hist_mid_interval, bracket);
      file << "}";
      file.close();
      /* -------------------------------------------------- */

      // find bin that has max bin size, i greater than index
      // for (auto i = index; i < Histo.bins.size(); ++i)
      //    if (Histo.bins[index].size() < Histo.bins[i].size())
      //       index = i;

      auto tmp = net.Lines;
      Divide(tmp, [&](auto l) { return l->length() > Histo.interval[index]; });

      /*
       * flipをしながらスムージングを行う
       */

      auto ps = net.getPoints();
#define remesh_type 4
      for (auto i = 0; i < 5; ++i) {
         for (auto j = 0; j < 50 /*more than about 10*/; ++j) {
            for (const auto &p : ps) {
               SmoothingPreserveShape(p, [&](auto p, auto X) {
#if remesh_type == 1
                  //! cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp -DOUTPUT_NAME=remeshA; make;./remeshA ./remesh_test_cases/random.obj ./remesh_test_cases randomA 10
                  return NeighborAverageSmoothingVector(p, X);
#elif remesh_type == 2
                      //! cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp -DOUTPUT_NAME=remeshB; make;./remeshB ./remesh_test_cases/random.obj ./remesh_test_cases randomB 10
                      return IncenterAverageSmoothingVector(p, X);
#elif remesh_type == 3
                      //! cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp -DOUTPUT_NAME=remeshC; make;./remeshC ./remesh_test_cases/random.obj ./remesh_test_cases randomC 10
                      return EquilateralVertexAveragingVector(p, X);
#elif remesh_type == 4
                      //! cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp -DOUTPUT_NAME=remeshD; make;./remeshD ./remesh_test_cases/random.obj ./remesh_test_cases randomD 10
                      return (0.8 * NeighborAverageSmoothingVector(p, X) + 0.1 * EquilateralVertexAveragingVector2(p, X));

#endif
               });
            }
            for (const auto &l : net.Lines) l->flipIfBetter(M_PI / 180.);
            if (i % 2 == 0 && j % 5 == 0)
               for (const auto &l : net.Lines) l->flipIfTopologicallyBetter(M_PI / 180., M_PI / 180.);
         }
      }
      if (count % 1 == 0) {
         std::string filename = output_directory + "/" + output_name + std::to_string(output_count) + ".vtu";
         mk_vtu(filename, net.getFaces(), getData());
         std::ofstream ofs(output_directory + "/" + output_name + std::to_string(output_count) + ".obj");
         createOBJ(ofs, net);
         ofs.close();
         //! PVD
         PVD.push(filename, count);
         PVD.output();
         output_count++;
      }

      std::cout << "\r \r" << std::flush;
   }
}
