/*DOC_EXTRACT 9_9_CGAL

# CGALを使って四面体を生成する

## CGALを使って四面体を生成する

```shell
brew install CGAL
```

*/

#include "Network.hpp"
#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"
//
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_3 CGALPoint;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;  // Defining Delaunay type

/* -------------------------------------------------------------------------- */

CGALPoint to_CGAL_point(const std::array<double, 3>& arr) {
   return CGALPoint(arr[0], arr[1], arr[2]);
}

CGALPoint to_CGAL_point(const networkPoint* p) {
   return to_CGAL_point(p->X);
}

std::array<double, 3> to_std_array(const CGALPoint& point) {
   return {point.x(), point.y(), point.z()};
}

/* -------------------------------------------------------------------------- */

int main() {
   /* -------------------------------------------------------------------------- */
   auto net = new Network("./bunny.obj", "duck");
   {
      std::ofstream ofs(_HOME_DIR_ + "/output_CGAL/polygon.vtp");
      vtkPolygonWrite(ofs, net->getFaces());
      ofs.close();
   }

   std::vector<CGALPoint> vertices;
   for (const auto& p : net->getPoints())
      vertices.push_back(to_CGAL_point(p));

   /*
   Delaunay_triangulation_3は，頂点が追加されると，自動的に四面体分割を行う．
   */
   Delaunay T;
   for (const auto& p : vertices)
      T.insert(p);

   /* -------------------------------------------------------------------------- */
   // std::ifstream input("./bunny.off");
   // Polyhedron poly;
   // input >> poly;
   // input.close();
   // Delaunay T;
   // for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) {
   //    CGALPoint p = v->point();
   //    T.insert(p);
   // }
   /* -------------------------------------------------------------------------- */

   std::cout << "Number of tetrahedra: " << T.number_of_finite_cells() << std::endl;
   std::unordered_map<CGALPoint, networkPoint*> CGALPoint2networkPoint;

   net->makeBucketFaces(0.5 * Mean(extLength(net->getLines())));
   net->BucketFaces.setVector();

   int count = 0;
   auto get_p = [&](const auto& CGALP) {
      auto [it, inserted] = CGALPoint2networkPoint.try_emplace(CGALP, nullptr);
      if (inserted)
         it->second = new networkPoint(net, to_std_array(CGALP));
      return it->second;
   };

   for (auto c = T.finite_cells_begin(); c != T.finite_cells_end(); ++c) {
      auto CGALP0 = c->vertex(0)->point();
      auto CGALP1 = c->vertex(1)->point();
      auto CGALP2 = c->vertex(2)->point();
      auto CGALP3 = c->vertex(3)->point();

      auto X0 = to_std_array(CGALP0);
      auto X1 = to_std_array(CGALP1);
      auto X2 = to_std_array(CGALP2);
      auto X3 = to_std_array(CGALP3);

      // if (net->isInside_MethodBucket(Incenter(X0, X1, X2, X3)))
      if (net->isInside_MethodBucket(Mean(X0, X1, X2, X3)))
         genTetra(net, get_p(CGALP0), get_p(CGALP1), get_p(CGALP2), get_p(CGALP3));
   }

   //
   {
      std::ofstream ofs(_HOME_DIR_ + "/output_CGAL/tetras.vtp");
      vtkPolygonWrite(ofs, net->getTetras());
      ofs.close();
   }

   return 0;
}
