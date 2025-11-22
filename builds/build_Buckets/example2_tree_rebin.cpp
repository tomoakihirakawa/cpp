/*DOC_EXTRACT 0_2_1_space_partitioning

## 階層のある空間分割（木構造）

```shell
sh clean; cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example2_tree_rebin.cpp
make
./example2_tree_rebin
```

### `rebin`

バケツに保存したオブジェクトの位置が変化した場合，
ツリー構造全体を再構築するのではなく，
バケツの外に出てしまったオブジェクトだけを再配置する方が効率的である．

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

PVDWriter pvdPoints0(home_dir + "/output/points_repeate_LV0.pvd");
PVDWriter pvdPoints1(home_dir + "/output/points_repeate_LV1.pvd");
PVDWriter pvdPoints2(home_dir + "/output/points_repeate_LV2.pvd");
PVDWriter pvdPoints3(home_dir + "/output/points_repeate_LV3.pvd");

PVDWriter pvdBuckets0(home_dir + "/output/buckets_repeate_LV0.pvd");
PVDWriter pvdBuckets1(home_dir + "/output/buckets_repeate_LV1.pvd");
PVDWriter pvdBuckets2(home_dir + "/output/buckets_repeate_LV2.pvd");
PVDWriter pvdBuckets3(home_dir + "/output/buckets_repeate_LV3.pvd");

std::vector<PVDWriter> pvdPoints = {pvdPoints0, pvdPoints1, pvdPoints2, pvdPoints3};
std::vector<PVDWriter> pvdBuckets = {pvdBuckets0, pvdBuckets1, pvdBuckets2, pvdBuckets3};

int current_level = 0, max_level = 3;
int main() {
   obj->makeBucketPoints(obj->getScale() / 4.);
   obj->makeBucketFaces(obj->getScale() / 4.);
   obj->BucketPoints.setLevel(current_level, max_level);
   obj->BucketPoints.generateTree();
   writePolygonToFile(home_dir + "/output/buckets_level0.vtp", obj->BucketPoints.getVertices());

   int repeat = 1;
   auto get = [&repeat]() {
      for (int level = 0; level <= max_level; level++) {
         std::cout << "Level: " << level << ", Repeat: " << repeat << std::endl;
         std::vector<T8Tddd> all_cube;
         std::vector<networkPoint*> all_points;
         std::unordered_map<networkPoint*, double> map;
         int l = 0;
         obj->BucketPoints.forEachAtLevel(level, [&](auto bucket) {
            T8Tddd cube = bucket->getVertices();
            all_cube.emplace_back(cube);
            for (auto p : bucket->data1D) {
               all_points.emplace_back(p);
               map[p] = (double)(l % 10);
            }
            l++;
         });
         std::string suffix = "_level" + std::to_string(level) + "_repeat" + std::to_string(repeat) + ".vtp";
         //
         auto name = home_dir + "/output/all_buckets" + suffix;
         writePolygonToFile(name, all_cube);
         pvdBuckets[level].push(name, repeat);
         //
         name = home_dir + "/output/all_points" + suffix;
         writePolygonToFile(name, all_points, map);
         pvdPoints[level].push(name, repeat);
      }

      std::cout << "Move points" << std::endl;
      for (auto& p : obj->getPoints())
         p->setXSingle(p->X + 0.02 * std::sin(0.1 * (2. * M_PI) * repeat));
      repeat++;

      std::cout << "rebin" << std::endl;
      auto escaped = obj->BucketPoints.rebin();
      if (!escaped.empty()) {
         obj->makeBucketPoints(obj->getScale() / 4.);
         obj->BucketPoints.setLevel(current_level, max_level);
         obj->BucketPoints.generateTree();
      }
   };

   for (auto i = 0; i < 50; i++)
      get();

   for (auto& pvd : pvdBuckets)
      pvd.output();
   for (auto& pvd : pvdPoints)
      pvd.output();
}