#include "Network.hpp"

/*DOC_EXTRACT REMESH

# Fusion360を使って計算用objファイルを生成

<img src="sample_fusion360_step1.png" width="600px">

<img src="sample_fusion360_step2.png" width="600px">

<img src="sample_fusion360_step3.png" width="600px">

<img src="sample_fusion360_step4.png" width="600px">

# 計算用にメッシュの細分化

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

int main(int arg, char **argv) {
   Timer time;
   std::string input_file{argv[1]};        // input
   std::string output_directory{argv[2]};  // output dir
   std::string output_name{argv[3]};       // output name
   int remesh = std::atoi(argv[4]);
   Network net(input_file, output_name);
   mk_vtu(output_directory + "/" + output_name + ".vtu", {net.getFaces()});
   net.displayStates();

   V_netLp lines;
   Histogram Histo;

   for (auto count = 0; count <= remesh; ++count) {
      if (count % 10 == 0) {
         mk_vtu(output_directory + "/" + output_name + std::to_string(count) + ".vtu", net.getFaces());
         std::ofstream ofs(output_directory + "/" + output_name + std::to_string(count) + ".obj");
         creteOBJ(ofs, net);
         ofs.close();
      }
      std::cout << "\r" << icons[(count + 1) % icons.size()] << " time:" << time() << std::flush;
      int num = 0;
      /* ------------------------*/
      //* | 0 | 1 | 2 |[3]| 4 |     diff, mid_interval
      //!   0   1   2  [3]  4   5   interval, cumulative_count
      /* ------------------------*/
      Histo.set(extLength(net.getLines()));
      for (auto j = 0; j < Histo.cumulative.size(); ++j)
         if (Histo.cumulative[j] >= 0.8) {
            num = j;
            break;
         }
      auto tmp = net.Lines;
      Divide(tmp, [&](auto l) { return l->length() > Histo.interval[num]; });

      for (auto i = 0; i < 5; i++) {
         AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         for (const auto &l : net.Lines)
            l->flipIfTopologicallyBetter(0.5 * M_PI / 180., 0.5 * M_PI / 180.);
         AreaWeightedSmoothingPreserveShape(net.getPoints(), small);
         for (const auto &l : net.Lines)
            l->flipIfBetter(M_PI / 180.);
      }
      LaplacianSmoothingPreserveShape(net.getPoints(), small);
      std::cout << "\r \r" << std::flush;
   }
}
