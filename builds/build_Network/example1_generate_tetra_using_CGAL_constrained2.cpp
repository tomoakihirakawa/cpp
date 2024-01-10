#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include "Network.hpp"  // Your custom Network class

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

int main() {
   // Load your polyhedron (representing the given triangle faces)
   Polyhedron polyhedron;
   std::ifstream input("./bunny.off");

   // Create a domain with features
   Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(polyhedron);

   // Define mesh criteria
   Mesh_criteria criteria(
       CGAL::parameters::edge_size = 0.01,
       CGAL::parameters::facet_angle = 5,
       CGAL::parameters::facet_size = 0.01,
       CGAL::parameters::facet_distance = 0.001,
       CGAL::parameters::cell_radius_edge_ratio = 3);

   // Mesh generation
   C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

   // Mesh refinement (optional)
   // CGAL::refine_mesh_3(c3t3, domain, criteria);

   // Exporting mesh data using your Network class
   Network network;
   network.processMesh(c3t3);               // Assume this method processes and stores mesh data
   network.exportToVTP("output_mesh.vtp");  // Export to VTP file
}
