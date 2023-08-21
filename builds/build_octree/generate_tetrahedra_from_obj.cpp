#include "Network.hpp"
#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"
//
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef K::Point_3 CGALPoint;

std::array<double, 3> to_std_array(const CGALPoint& point) {
   return {point.x(), point.y(), point.z()};
}

int main() {
   std::ifstream input("./bunny.off");
   Polyhedron poly;
   input >> poly;
   input.close();

   Delaunay T;
   for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) {
      CGALPoint p = v->point();
      T.insert(p);
   }

   std::cout << "Number of tetrahedra: " << T.number_of_finite_cells() << std::endl;

   // output
   auto net = new Network;
   auto object = new Network("./bunny.off", "bunny");
   object->genOctreeOfFaces({1, 5}, 1);

   for (auto c = T.finite_cells_begin(); c != T.finite_cells_end(); ++c) {

      auto X0 = to_std_array(c->vertex(0)->point());
      auto X1 = to_std_array(c->vertex(1)->point());
      auto X2 = to_std_array(c->vertex(2)->point());
      auto X3 = to_std_array(c->vertex(3)->point());

      auto [isinside, _, __] = object->isInside_MethodOctree(Incenter(X0, X1, X2, X3));

      if (isinside) {
         auto p0 = new networkPoint(net, X0);
         auto p1 = new networkPoint(net, X1);
         auto p2 = new networkPoint(net, X2);
         auto p3 = new networkPoint(net, X3);
         genTetra(net, p0, p1, p2, p3);
      }
   }

   //
   std::ofstream ofs(_HOME_DIR_ + "/output_CGAL/tetras.vtp");
   vtkPolygonWrite(ofs, net->getTetras());
   ofs.close();

   return 0;
}
