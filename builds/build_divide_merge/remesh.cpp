
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
$ sh clean
$ cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp
$ make
```

`remesh`という実行ファイルができる．

## 実行方法

`n`回の細分化を行う．

```
./remesh input_file output_dir output_name n
```

![./sample2.gif](./sample2.gif)

出力は，`output_dir/output_name*.vtu`と`output_dir/output_name*.obj`．


![./sample.gif](./sample.gif)

次は，入力ファイルを生成し，計算をする．

*/

const std::array<std::string, 4> icons = {Red + '|' + colorOff, Magenta + '/' + colorOff, Blue + '-' + colorOff, Green + '\\' + colorOff};
const double small = M_PI / 180 * 0.01;

JSONoutput output_json;

#define remesh_type 5

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

   for (auto count = 0; count <= remesh; ++count) {

      std::cout << "\r" << icons[(count + 1) % icons.size()] << " time:" << time() << std::flush;
      int num = 0;
      /* ------------------------*/
      //* | 0 | 1 | 2 |[3]| 4 |     diff, mid_interval
      //!   0   1   2  [3]  4   5   interval, cumulative_count
      /* ------------------------*/
      Histo.set(extLength(net.getLines()));
      auto it = std::ranges::find_if(Histo.cumulative, [](double value) { return value >= 0.8; });
      if (it != Histo.cumulative.end()) num = std::distance(Histo.cumulative.begin(), it);

      // find bin that has max bin size, i greater than num
      // for (auto i = num; i < Histo.bins.size(); ++i)
      //    if (Histo.bins[num].size() < Histo.bins[i].size())
      //       num = i;

      auto tmp = net.Lines;
      Divide(tmp, [&](auto l) { return l->length() > Histo.interval[num]; });

      for (auto i = 0; i < 3; ++i) {
#if remesh_type == 1
         LaplacianSmoothingPreserveShape(net.getPoints(), 2);
         for (const auto &l : net.Lines) l->flipIfTopologicallyBetter(M_PI / 180., M_PI / 180.);
         LaplacianSmoothingPreserveShape(net.getPoints(), 2);
         for (const auto &l : net.Lines) l->flipIfBetter(M_PI / 180.);
         LaplacianSmoothingPreserveShape(net.getPoints(), 2);
#elif remesh_type == 2
         AreaWeightedSmoothingPreserveShape(net.getPoints(), 2);
         for (const auto &l : net.Lines) l->flipIfTopologicallyBetter(M_PI / 180., M_PI / 180.);
         AreaWeightedSmoothingPreserveShape(net.getPoints(), 2);
         for (const auto &l : net.Lines) l->flipIfBetter(M_PI / 180.);
         AreaWeightedSmoothingPreserveShape(net.getPoints(), 2);
#elif remesh_type == 3
         DistorsionMeasureWeightedSmoothingPreserveShape(net.getPoints(), 2);
         for (const auto &l : net.Lines) l->flipIfTopologicallyBetter(M_PI / 180., M_PI / 180.);
         DistorsionMeasureWeightedSmoothingPreserveShape(net.getPoints(), 2);
         for (const auto &l : net.Lines) l->flipIfBetter(M_PI / 180.);
         DistorsionMeasureWeightedSmoothingPreserveShape(net.getPoints(), 2);
#elif remesh_type == 4
         DistorsionMeasureWeightedSmoothingPreserveShape(net.getPoints(), 3);
         for (const auto &l : net.Lines) l->flipIfTopologicallyBetter(M_PI / 180., M_PI / 180.);
         DistorsionMeasureWeightedSmoothingPreserveShape(net.getPoints(), 3);
         for (const auto &l : net.Lines) l->flipIfBetter(M_PI / 180.);
         LaplacianSmoothingPreserveShape(net.getPoints(), 3);
#elif remesh_type == 5
         LaplacianSmoothingPreserveShape(net.getPoints(), 5);
         if (count < 50) {
            for (const auto &l : net.Lines) l->flipIfBetter(M_PI / 180.);
            LaplacianSmoothingPreserveShape(net.getPoints(), 5);
         } else if (count % 3 == 0) {
            for (const auto &l : net.Lines) l->flipIfTopologicallyBetter(M_PI / 180., M_PI / 180.);
            AreaWeightedSmoothingPreserveShape(net.getPoints(), 5);
            for (const auto &l : net.Lines) l->flipIfBetter(M_PI / 180.);
            DistorsionMeasureWeightedSmoothingPreserveShape(net.getPoints(), 5);
         } else {
            for (const auto &l : net.Lines) l->flipIfBetter(M_PI / 180.);
            DistorsionMeasureWeightedSmoothingPreserveShape(net.getPoints(), 5);
         }
#endif
      }
      if (count % 10 == 0) {
         std::string step = std::to_string(count);
         std::string filename = output_directory + "/" + output_name + step + ".vtu";
         mk_vtu(filename, net.getFaces(), getData());
         std::ofstream ofs(output_directory + "/" + output_name + step + ".obj");
         createOBJ(ofs, net);
         ofs.close();
         //! PVD
         PVD.push(filename, count);
         PVD.output();
         //
         //! insert histogram data into json
         {
            Histogram Histo;
            Histo.set(extLength(net.getLines()));
            for (const auto v : Histo.bins)
               output_json.push("size" + step, (int)v.size());
            for (const auto v : Histo.mid_interval)
               output_json.push("mid_interval" + step, v);
            std::ofstream os("output_line" + std::to_string(remesh_type) + ".json");
            output_json.output(os);
            os.close();
         }
         {
            Histogram Histo;
            std::vector<double> tmp;
            for (const auto &f : net.getFaces())
               tmp.push_back(CircumradiusToInradius(ToX(f)));
            Histo.set(tmp);
            for (const auto v : Histo.bins)
               output_json.push("size" + step, (int)v.size());
            for (const auto v : Histo.mid_interval)
               output_json.push("mid_interval" + step, v);
            std::ofstream os("output_face" + std::to_string(remesh_type) + ".json");
            output_json.output(os);
            os.close();
         }
      }

      std::cout << "\r \r" << std::flush;
   }
}
