
/*DOC_EXTRACT 0_2_space_partitioning

## 階層のある空間分割（木構造）

 シンプルな空間分割クラスを拡張し，木構造による空間分割を試みる．

`has_tree`が`true`の場合，`buckets`には`Bucket`クラスのポインタが格納される．
`buckets[i][j][k]`には，上のレベルの`data[i][j][k]`のデータが引き継がれている．
つまり，`buckets[i][j][k]`は，`data[i][j][k]`のデータをさらに分割したものである．
デフォルトでは，`buckets[i][j][k]`は内部に８つの`data`を持つ:

`data[0][0][0]`，`data[0][0][1]`，`data[0][1][0]`，`data[0][1][1]`，`data[1][0][0]`，`data[1][0][1]`，`data[1][1][0]`，`data[1][1][1]`．

\ref{buckets_generateTree}{このツリー生成方法}は，
バウンディングボックスを範囲と，それを分割する幅を指定する．
分割数を指定するよりも，この方法のように分割幅を指定する方が，自分はわかりやすい．

```cpp
buckets[i][j][k] = std::make_shared<Buckets<T>>(bounds, this->dL * 0.5 + 1e-10);
```

![example2_tree_faster.gif](example2_tree_faster.gif)

*/

#include "Network.hpp"
#include "vtkWriter.hpp"

auto obj = std::make_unique<Network>("./bunny.obj");

// Function to write polygons to a file
void writePolygonToFile(const std::string& name, auto polygon) {
   std::ofstream ofs(name);
   vtkPolygonWrite(ofs, polygon);
   ofs.close();
}

// Function to process deeper levels of the tree
void processDeeperLevels(const Tiii ijk, const int l, int& global_l) {
   auto [i, j, k] = ijk;
   if (obj->BucketPoints.has_tree && obj->BucketPoints.buckets[i][j][k] != nullptr) {
      int L = 0;
      auto& bucket = obj->BucketPoints.buckets[i][j][k];
      auto& data = bucket->data;
      int level = bucket->level;

      // Further subdivision
      for (int i = 0; i < data.size(); i++) {
         for (int j = 0; j < data[i].size(); j++) {
            for (int k = 0; k < data[i][j].size(); k++) {
               std::string suffix = "_level" + std::to_string(level) + "_" + std::to_string(l) + "_" + std::to_string(L) + ".vtp";
               std::string suffix2 = "_all_level" + std::to_string(level) + "_" + std::to_string(global_l) + ".vtp";

               writePolygonToFile("./output/bucket_range" + suffix, (T8Tddd)(bucket->getBounds({i, j, k})));
               writePolygonToFile("./output/bucket_range" + suffix2, (T8Tddd)(bucket->getBounds({i, j, k})));
               writePolygonToFile("./output/points_in_bucket" + suffix, bucket->data[i][j][k]);
               writePolygonToFile("./output/points_in_bucket" + suffix2, bucket->data[i][j][k]);

               global_l++;
               L++;
            }
         }
      }
   }
}

int main() {
   obj->makeBucketPoints(obj->getScale() / 4.);
   obj->makeBucketFaces(obj->getScale() / 4.);

   obj->BucketPoints.generateTree();

   int l = 0, global_l = 0;
   PVDWriter pvd("./bunny_buckets.pvd");
   std::vector<T8Tddd> all_cube;

   auto& p_data_level0 = obj->BucketPoints.data;
   auto& f_data_level0 = obj->BucketFaces.data;

   for (int i = 0; i < p_data_level0.size(); i++) {
      for (int j = 0; j < p_data_level0[i].size(); j++) {
         for (int k = 0; k < p_data_level0[i][j].size(); k++) {
            std::string suffix = "_level0_" + std::to_string(l) + ".vtp";

            auto cube = (T8Tddd)obj->BucketPoints.getBounds({i, j, k});
            all_cube.emplace_back(cube);
            writePolygonToFile("./output/bucket_range" + suffix, cube);
            writePolygonToFile("./output/points_in_bucket" + suffix, p_data_level0[i][j][k]);
            writePolygonToFile("./output/faces_in_bucket" + suffix, f_data_level0[i][j][k]);

            processDeeperLevels({i, j, k}, l, global_l);

            std::cout << "i: " << i << ", j: " << j << ", k: " << k << ", p_data_level0[i][j][k].size(): " << p_data_level0[i][j][k].size() << std::endl;

            l++;
         }
      }
   }

   writePolygonToFile("./output/bucket_all_cube.vtp", all_cube);
}
