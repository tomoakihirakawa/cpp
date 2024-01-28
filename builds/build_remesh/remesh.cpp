
#include <filesystem>
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
   std::string input_file{argv[1]};                  // input
   std::filesystem::path output_directory{argv[2]};  // output dir
   std::string output_name{argv[3]};                 // output name
   std::filesystem::path output_name_vtu = output_name + ".vtu";
   std::filesystem::path output_name_pvd = output_name + ".pvd";
   std::filesystem::path inner_output_name_pvd = output_name + "_inner.pvd";
   int remesh = std::atoi(argv[4]);
   Network net(input_file, output_name);
   mk_vtu(output_directory / output_name_vtu, {net.getFaces()});
   net.displayStates();

   PVDWriter PVD(output_directory / output_name_pvd);
   PVDWriter PVD_inner(output_directory / inner_output_name_pvd);

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
      // auto it_cumulative = std::ranges::find_if(Histo.cumulative, [](double value) { return value >= 0.9; });
      // if (it_cumulative != Histo.cumulative.end())
      //    index = std::distance(Histo.cumulative.begin(), it_cumulative);
      // auto it_count = std::max_element(Histo.count.begin() + index, Histo.count.end());
      // if (it_count != Histo.count.end())
      //    index = std::distance(Histo.count.begin(), it_count);

      // copy a histogram of the length of the line, then sort it

      double value = 1E+20;
      auto findLargeJamp = [&](const double longer_than) {
         auto temp = Histo.diff;
         std::sort(temp.begin(), temp.end(), [](auto a, auto b) { return a > b; });
         auto it = Histo.count.begin();
         for (auto i = 1; i < temp.size(); ++i) {
            auto IT = std::find(Histo.diff.begin(), Histo.diff.end(), temp[i]);
            auto INDEX = std::distance(Histo.diff.begin(), IT) + 1;
            if (Histo.interval[INDEX] > longer_than) {
               index = INDEX;
               value = Histo.interval[INDEX];
               break;
            }
         }
      };

      findLargeJamp(0.01);

      /* -------------------------------------------------- */
      // Output to JSON file
      hist_where_divide.push_back(Histo.interval[index]);
      std::cout << "count: " << Histo.count << std::endl;
      std::cout << "interval: " << Histo.interval << std::endl;
      std::cout << "index: " << index << std::endl;
      std::cout << "interval"
                << "[" << index << "]"
                << ": " << Histo.interval[index] << std::endl;

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
      if (count > 0)
         Divide(tmp, [&](auto l) {
            // return l->length() > 0.005 && l->length() > Histo.interval[index];
            // return l->length() > Histo.interval[index];
            return l->length() > value;
         });

      /*
       * flipをしながらスムージングを行う
       */

      auto ps = ToVector(net.getPoints());
      int inner_loop_count = 0;
      int c = 0;

      double rad = M_PI / 180;
      // flipIf(net, {10 * rad, 10 * rad}, {10 * rad, 10 * rad}, true);
      // flipIf(net, {10 * rad, 10 * rad}, {10 * rad, 10 * rad}, true);

      //! この値の変化は大きな影響を与える
      const double limit_angle_to_distinguish_flat = 1E-2 * M_PI / 180.;
      for (auto i = 0; i < 5; ++i) {

         for (const auto &l : net.Lines)
            l->flipIfBetter(limit_angle_to_distinguish_flat);

         for (const auto &l : net.Lines)
            l->flipIfBetter(limit_angle_to_distinguish_flat);

         flipIf(net, {0.1 * rad, 0.1 * rad}, true);

         // for (const auto &p : ps) {
         //    if (!isFlat(p)) {
         //       flipIf(p, {0.1 * rad, 0.1 * rad}, true, 6);
         //       for (auto &q : p->getNeighbors())
         //          flipIf(q, {0.1 * rad, 0.1 * rad}, true, 6);
         //    }
         // }

         for (auto j = 0; j < 5 /*more than about 10*/; ++j) {
            net.setGeometricProperties();
#pragma omp parallel
            for (const auto &p : ps)
#pragma omp single nowait
               SmoothingPreserveShape(p);

            for (const auto &p : ps) {
               auto X = p->X;
               bool is_falt_before = isFlat(p);
               // p->setX(0.1 * (p->X_temporary - p->X) + p->X);
               p->setX(p->X_temporary);
               bool is_falt_after = isFlat(p);
               if (is_falt_before != is_falt_after)
                  p->setX(X);
            }

            net.setGeometricProperties();
         }
      }
      if (count % 1 == 0) {
         std::filesystem::path output_name_vtu = output_name + std::to_string(output_count) + ".vtu";
         std::filesystem::path output_name_obj = output_name + std::to_string(output_count) + ".obj";
         mk_vtu(output_directory / output_name_vtu, net.getFaces(), getData());

         std::ofstream ofs(output_directory / output_name_obj);
         createOBJ(ofs, net);
         ofs.close();
         //! PVD
         PVD.push(output_directory / output_name_vtu, count);
         PVD.output();
         output_count++;
      }

      std::cout << "\r \r" << std::flush;
   }
}

//       for (auto i = 0; i < 2; ++i) {

//          for (const auto &l : net.Lines)
//             l->flipIfBetter(1E-5 * M_PI / 180.);

//          // flipIf(net, {0.1 * rad, 0.1 * rad}, true);

//          /*more than about 10*/
//          for (auto j = 0; j < 5; ++j) {

//             // for (const auto &p : ps)
//             //    if (!isFlat(p))
//             //       flipIf(p, {0.1 * rad, 0.1 * rad}, true);

//             net.setGeometricProperties();
//             for (const auto &p : ps)
//                for (const auto &f : p->getFaces())
//                   f->normal_to_be_preserved = f->normal;

//             for (const auto &p : ps) {
//                p->normal_to_be_preserved.fill(0.);
//                for (const auto &f : p->getFaces())
//                   FusedMultiplyIncrement(f->area, f->normal, p->normal_to_be_preserved);
//                p->normal_to_be_preserved = Normalize(p->normal_to_be_preserved);

//                if (std::ranges::all_of(p->getFaces(), [&](auto f) { return isFlat(f->normal, p->normal_to_be_preserved, limit_angle_to_distinguish_flat); }))
//                   for (const auto &f : p->getFaces())
//                      f->normal_to_be_preserved = Normalize(p->normal_to_be_preserved);
//             }
// #pragma omp parallel
//             for (const auto &p : ps)
// #pragma omp single nowait
//                SmoothingPreserveShape(p);

//             for (const auto &p : ps) {
//                auto X = p->X;
//                bool is_falt_before = isFlat(p);
//                p->setX(p->X_temporary);
//                bool is_falt_after = isFlat(p);
//                if (is_falt_before != is_falt_after)
//                   p->setX(X);
//             }

//             net.setGeometricProperties();
//          }
//       }
