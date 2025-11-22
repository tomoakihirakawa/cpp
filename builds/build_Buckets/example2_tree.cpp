
/*DOC_EXTRACT 0_2_0_space_partitioning

## 階層のある空間分割（木構造）

```shell
sh clean;
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example2_tree.cpp
make
./example2_tree
```

### `generateTree`

ツリー生成`generateTree`によって`children[i][j][k]`が生成されると，`data[i][j][k]`の内容が`children[i][j][k]->data1D`にコピーされる．

以下二つは全く同じものを保存している:

* `b->data[i][j][k]`
* `b->children[i][j][k]->data1D`　


<div class="slide-container">
<img src="example2_level0.png" width="500px" />
<img src="example2_level1.png" width="500px" />
<img src="example2_level2.png" width="500px" />
<img src="example2_level3.png" width="500px" />
</div>

levelを０から３まで変化させて，各レベルのBucketsクラス`children`にオブジェクトを再帰的に振り分けていく．

*/

#include "Network.hpp"
#include "vtkWriter.hpp"

auto obj = std::make_unique<Network>("./bunny.obj");

// HOME_DIR
std::string home_dir = getenv("HOME");

// Function to write polygons to a file
void writePolygonToFile(const std::string& name, auto polygon) {
   std::ofstream ofs(name);
   vtkPolygonWrite(ofs, polygon);
   ofs.close();
   std::cout << "Wrote " << name << std::endl;
}

void writePolygonToFile(const std::string& name, auto polygon, std::unordered_map<networkPoint*, double>& map) {
   std::ofstream ofs(name);
   vtkPolygonWrite(ofs, polygon, map);
   ofs.close();
   std::cout << "Wrote " << name << std::endl;
}

int current_level = 0, max_level = 3;
int main() {
   obj->makeBucketPoints(obj->getScale() / 4.);
   obj->makeBucketFaces(obj->getScale() / 4.);
   obj->BucketPoints.setLevel(current_level, max_level);
   obj->BucketPoints.generateTree();

   auto get = [&](int level) {
      std::vector<T8Tddd> all_cube;
      std::vector<networkPoint*> all_points;
      std::unordered_map<networkPoint*, double> map;
      int l = 0;
      obj->BucketPoints.forEachAtLevel(level, [&](auto bucket) {
         std::string suffix = "_level" + std::to_string(level) + "_" + std::to_string(l) + ".vtp";
         T8Tddd cube = bucket->getVertices();
         all_cube.emplace_back(cube);
         for (auto p : bucket->data1D) {
            all_points.emplace_back(p);
            map[p] = (double)(l % 10);
         }
         l++;
      });
      writePolygonToFile(home_dir + "/output/all_buckets_level" + std::to_string(level) + ".vtp", all_cube);
      writePolygonToFile(home_dir + "/output/all_points_level" + std::to_string(level) + ".vtp", all_points, map);
      std::cout << "Level " << level << " complete" << std::endl;
   };

   get(0);
   get(1);
   get(2);
   get(3);
}