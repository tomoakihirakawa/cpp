#include "Network.hpp"
#include "vtkWriter.hpp"

int main() {

   /*DOC_EXTRACT 0_1_Network_space_partitioning

   # 空間分割（space_partitioning）

   ## 等間隔のシンプルな空間分割

   ```shell
   cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_space_partitioning.cpp
   make
   ./example1_space_partitioning
   ```

   \insert{coordinatebounds}

   \insert{space_partitioning}

   ### 例

   この例では，うさぎの３Dモデルを空間分割する．
   配列させたバケット内に，うさぎの点または面が含まれるかを判定し，バケットに保存する．

   ただ，面は広がりがあるので，複数のバケットに含まれることがある．
   面と交わる全バケットを簡単に確実に見つける方法は，現在のところ思いつかない．
   なので，今の所は，面を無数の点に分けて，各点を含むバケットに面を保存することで対応している．

   ![example1_space_partitioning.gif](example1_space_partitioning.gif)

   */

   auto obj = new Network("./bunny.obj");
   obj->makeBucketPoints(obj->getScale() / 4.);
   obj->makeBucketFaces(obj->getScale() / 4.);
   int l = 0;
   PVDWriter pvd("./bunny_buckets.pvd");
   std::vector<T8Tddd> all_faces;
   auto& p_Bucket = obj->BucketPoints.data;
   auto& f_Bucket = obj->BucketFaces.data;
   for (auto i = 0; i < p_Bucket.size(); i++) {
      for (auto j = 0; j < p_Bucket[i].size(); j++) {
         for (auto k = 0; k < p_Bucket[i][j].size(); k++) {
            {
               auto name = "./output/bucket_range" + std::to_string(l) + ".vtp";
               std::ofstream ofs(name);
               T8Tddd face = obj->BucketPoints.getBounds({i, j, k});  // implicit conversion from CoordinateBounds to T8Tddd
               all_faces.emplace_back(face);
               vtkPolygonWrite(ofs, face);
               ofs.close();
            }
            {  // points in p_Bucket
               auto name = "./output/points_in_bucket" + std::to_string(l) + ".vtp";
               std::ofstream ofs(name);
               vtkPolygonWrite(ofs, p_Bucket[i][j][k]);
               ofs.close();
            }
            {  // face in p_Bucket
               auto name = "./output/faces_in_bucket" + std::to_string(l) + ".vtp";
               std::ofstream ofs(name);
               vtkPolygonWrite(ofs, f_Bucket[i][j][k]);
               ofs.close();
            }
            std::cout << "i: " << i << ", j: " << j << ", k: " << k << ", p_Bucket[i][j][k].size(): " << p_Bucket[i][j][k].size() << std::endl;
            l++;
         }
      }
   }
   std::ofstream ofs("./output/bucket_all.vtp");
   vtkPolygonWrite(ofs, all_faces);
   ofs.close();
};