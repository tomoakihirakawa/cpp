#include "Network.hpp"
#include "vtkWriter.hpp"

/*DOC_EXTRACT 0_4_face_to_face_contact

### 面同士の接触判定

\ref{basic_geometry:IntersectQ}{`IntersectQ`}関数は，交差判定には使えるが，接触判定には使えない．

接触は，ギリギリ交差している状態を指すだろうが，
実際に接触判定を応用する場面では，
交差していなくとも接触していると判定させたい場合が多いだろう．
なので，接触判定条件はより緩く設定されることが多い．

```shell
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example3_line_face_interaction.cpp
make
./example3_line_face_interaction
```

![./example3/anim.gif](./example3/anim_faster.gif)

*/

std::string outdir = "./example4/";
PVDWriter pvd_cube(outdir + "cube.pvd");
PVDWriter pvd_cylinder(outdir + "cylinder.pvd");
PVDWriter pvd_triangle(outdir + "triangle.pvd");
PVDWriter pvd_sphere(outdir + "sphere.pvd");

int main() {

   auto center = new Network("./example4/center3.obj");
   center->makeBucketFaces(center->getScale() / 10.);
   auto cube = new Network("./example4/cube3.obj");
   cube->makeBucketFaces(cube->getScale() / 10.);
   auto cylinder = new Network("./example4/cylinder3.obj");
   cylinder->makeBucketFaces(cylinder->getScale() / 10.);
   auto triangle = new Network("./example4/triangle3.obj");
   triangle->makeBucketFaces(triangle->getScale() / 10.);
   auto sphere = new Network("./example4/sphere3.obj");
   sphere->makeBucketFaces(sphere->getScale() / 10.);

   //! gradually move toward the x-direction
   const int steps = 100;
   const double distance = 0.2;
   for (int i = 0; i < steps; i++) {
      double x = -distance / steps;
      cube->translate({x, 0, 0});
      cylinder->translate({x, 0, 0});
      triangle->translate({x, 0, 0});
      sphere->translate({x, 0, 0});

      {
         std::string name = "cube" + std::to_string(i) + ".vtp";
         std::ofstream ofs(outdir + name);
         vtkPolygonWrite(ofs, cube->getFaces());
         ofs.close();
         pvd_cube.push(name, i);
      }
      {
         std::string name = "cylinder" + std::to_string(i) + ".vtp";
         std::ofstream ofs(outdir + name);
         vtkPolygonWrite(ofs, cylinder->getFaces());
         ofs.close();
         pvd_cylinder.push(name, i);
      }
      {
         std::string name = "triangle" + std::to_string(i) + ".vtp";
         std::ofstream ofs(outdir + name);
         vtkPolygonWrite(ofs, triangle->getFaces());
         ofs.close();
         pvd_triangle.push(name, i);
      }
      {
         std::string name = "sphere" + std::to_string(i) + ".vtp";
         std::ofstream ofs(outdir + name);
         vtkPolygonWrite(ofs, sphere->getFaces());
         ofs.close();
         pvd_sphere.push(name, i);
      }
   }
   pvd_cube.output();
   pvd_cylinder.output();
   pvd_triangle.output();
   pvd_sphere.output();
}