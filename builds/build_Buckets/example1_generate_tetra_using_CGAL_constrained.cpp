/*DOC_EXTRACT 9_9_CGAL

## CGALを使って四面体を生成する

WARNING:　コンパイルできない

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_generate_tetra_using_CGAL_constrained.cpp
make
```

`CGAL::Mesh_polyhedron_3<K>::type` is typically a typedef for a polyhedron data structure that is compatible with CGAL's mesh generation algorithms.
`CGAL::Polyhedron_3<K>` is a standard CGAL polyhedron class.

*/

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <fstream>
#include "CGAL_Network.hpp"
#include "Network.hpp"
#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"
//
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
// Triangulation for Meshing
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Weighted_point_3<CGAL::Epick> CGAL_Weight_Point;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

std::array<double, 3> to_std_array(const auto& point) {
   return {point.x(), point.y(), point.z()};
}

int main() {
   // Load input surface mesh

   // for (const auto& detect_angle : {150, 140, 130, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10})
   {
      // std::cout << "detect_angle = " << detect_angle << std::endl;
      std::string filename = "./bunny.obj";
      std::ifstream input(filename);
      auto net = new Network(filename);
      auto polyhedron = buildMeshPolyhedronFromNetwork(*net);
      std::cout << "polyhedron number of vertices = " << polyhedron.size_of_vertices() << std::endl;
      // Create domain
      Mesh_domain domain(polyhedron);
      // domain.detect_features(detect_angle);  // 例：60度以上のエッジを特徴として扱う
      // std::cout << "domain number of features = " << domain.number_of_features() << std::endl;

      // Mesh criteria
      Mesh_criteria criteria(CGAL::parameters::edge_size = 0.001,
                             CGAL::parameters::facet_angle = 1,
                             CGAL::parameters::facet_size = 0.01,
                             CGAL::parameters::facet_distance = 0.0004,
                             CGAL::parameters::cell_radius_edge_ratio = 3);

      // Mesh generation
      C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

      // Tetrahedral remeshing
      const double target_edge_length = 0.1;                                                                                       // Adjust as needed
      CGAL::tetrahedral_isotropic_remeshing(c3t3.triangulation(), target_edge_length, CGAL::parameters::number_of_iterations(3));  // Adjust iterations as needed

      // Output the remeshed mesh
      // Iterate over cells (tetrahedra) in the complex
      std::unordered_map<CGAL_Weight_Point, networkPoint*> CGALPoint2networkPoint;
      auto get_p = [&](const CGAL_Weight_Point& CGALP) {
         auto [it, inserted] = CGALPoint2networkPoint.try_emplace(CGALP, nullptr);
         if (inserted)
            it->second = new networkPoint(net, to_std_array(CGALP));
         return it->second;
      };

      for (auto cit = c3t3.cells_begin(); cit != c3t3.cells_end(); ++cit) {
         auto CGALP0 = cit->vertex(0)->point();
         auto CGALP1 = cit->vertex(1)->point();
         auto CGALP2 = cit->vertex(2)->point();
         auto CGALP3 = cit->vertex(3)->point();

         auto X0 = to_std_array(CGALP0);
         auto X1 = to_std_array(CGALP1);
         auto X2 = to_std_array(CGALP2);
         auto X3 = to_std_array(CGALP3);
         //   if (net->isInside_MethodBucket(Mean(X0, X1, X2, X3)))
         genTetra(net, get_p(CGALP0), get_p(CGALP1), get_p(CGALP2), get_p(CGALP3));
      }

      {
         // std::ofstream ofs("./tetras" + std::to_string(detect_angle) + ".vtp");
         std::ofstream ofs("./tetras.vtp");
         vtkPolygonWrite(ofs, net->getTetras());
         ofs.close();
         std::cout << "paraview tetras.vtp" << std::endl;
      }
   }
   return EXIT_SUCCESS;
}
