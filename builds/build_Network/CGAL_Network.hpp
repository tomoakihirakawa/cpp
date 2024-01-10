#ifndef CGAL_Network_H
#define CGAL_Network_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <unordered_map>
#include "Network.hpp"
#include "Network.hpp"  // Your Network class

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Mesh_Polyhedron;
typedef Mesh_Polyhedron::Point_3 Point_3;
typedef Mesh_Polyhedron::HalfedgeDS HalfedgeDS;

// Builder class to build Mesh_Polyhedron from Network
class BuildMeshPolyhedron : public CGAL::Modifier_base<HalfedgeDS> {
  public:
   BuildMeshPolyhedron(const Network& net) : net_(net) {}

   void operator()(HalfedgeDS& hds) {
      CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(hds, true);
      builder.begin_surface(net_.getPoints().size(), net_.getFaces().size());

      std::unordered_map<networkPoint*, int> vhandles;
      int index = 0;

      // Add vertices to Mesh_Polyhedron
      for (const auto& p : net_.getPoints()) {
         auto [x, y, z] = p->X;
         vhandles[p] = index++;
         builder.add_vertex(Point_3(x, y, z));
      }

      // Add faces to Mesh_Polyhedron
      for (const auto& f : net_.getFaces()) {
         auto [p0, p1, p2] = f->getPoints();
         builder.begin_facet();
         builder.add_vertex_to_facet(vhandles[p0]);
         builder.add_vertex_to_facet(vhandles[p1]);
         builder.add_vertex_to_facet(vhandles[p2]);
         builder.end_facet();
      }

      builder.end_surface();
   }

  private:
   const Network& net_;
};

// Function to build Mesh_Polyhedron from Network
Mesh_Polyhedron buildMeshPolyhedronFromNetwork(const Network& net) {
   Mesh_Polyhedron P;
   std::cout << "Building Mesh_Polyhedron from Network..." << std::endl;
   BuildMeshPolyhedron builder(net);
   std::cout << "Building BuildMeshPolyhedrone" << std::endl;
   P.delegate(builder);
   std::cout << "Building delegate builder... done" << std::endl;
   return P;
}

#endif
