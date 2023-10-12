#include "Network.hpp"
#include "vtkWriter.hpp"

int main() {

   /*DOC_EXTRACT Network_space_partitioning

   # 木構造による空間分割

    シンプルな空間分割クラスを拡張し，木構造による空間分割を試みる．

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