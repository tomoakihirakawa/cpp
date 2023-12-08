#include "Network.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 0_3_line_face_interaction

# 線と三角形の干渉判定

## 等間隔のシンプルな空間分割

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example3_line_face_interaction.cpp
make
./example3_line_face_interaction
```

*/

int main() {

   auto obj = new Network("./bunny.obj");
   obj->makeBucketPoints(obj->getScale() / 4.);
   obj->makeBucketFaces(obj->getScale() / 4.);
   int l = 0;
   PVDWriter pvd("./bunny_buckets.pvd");
   std::vector<T8Tddd> all_faces;
   auto& p_Bucket = obj->BucketPoints.data;
   auto& f_Bucket = obj->BucketFaces.data;

   std::ofstream ofs("./output/bucket_all.vtp");
   vtkPolygonWrite(ofs, all_faces);
   ofs.close();
};