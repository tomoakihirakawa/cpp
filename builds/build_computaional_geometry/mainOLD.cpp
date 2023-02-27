#include "basic_IO.hpp"
#include "basic.hpp"
#include "lib_computationalGeometry.hpp"
namespace gn = geometric_network;
#include "vtkWriter.hpp"
#include "lib_measurement.hpp"

#define case_simple2
#if defined(case_simple2)
//
// Points --> Tetrahedrons --> Edges and Triangles
//
int main()
{
    Timer timer;
    auto net = new gn::Network;
    auto p0 = new gn::Point(net, 0., 0., 0.);
    auto p1 = new gn::Point(net, 1., 0., 0.);
    auto p2 = new gn::Point(net, 1., 1., 0.);
    auto p3 = new gn::Point(net, 1., 1., 1.);
    auto p4 = new gn::Point(net, 1., 1., -1.);
    auto p5 = new gn::Point(net, 1., -1., -1.);
    std::cout << "time:" << timer() << std::endl;
    //
    auto tet0 = new gn::Tetrahedron(net, p0, p1, p2, p3);
    tet0->fill(tet0->getPoints());
    auto tet1 = new gn::Tetrahedron(net, p0, p1, p2, p4);
    tet1->fill(tet1->getPoints());
    //
    auto tet2 = new gn::Tetrahedron(net);
    auto e05 = new gn::Edge(net, p0, p5);
    auto e15 = new gn::Edge(net, p1, p5);
    auto e25 = new gn::Edge(net, p2, p5);
    tet2->fill(tet0->E0, tet0->E1, tet0->E2, e05, e15, e25);
    std::cout << "time:" << timer() << std::endl;
    /* ------------------------------------------------------ */
    vtkPolygonWriter<gn::Point *> vtp;
    vtp.add({p0, p1, p2, p3, p4, p5});
    vtp.addPolygon(tet0->T0->getPoints());
    vtp.addPolygon(tet0->T1->getPoints());
    vtp.addPolygon(tet0->T2->getPoints());
    vtp.addPolygon(tet0->T3->getPoints());
    //
    vtp.addPolygon(tet1->T0->getPoints());
    vtp.addPolygon(tet1->T1->getPoints());
    vtp.addPolygon(tet1->T2->getPoints());
    vtp.addPolygon(tet1->T3->getPoints());
    //
    vtp.addPolygon(tet2->T0->getPoints());
    vtp.addPolygon(tet2->T1->getPoints());
    vtp.addPolygon(tet2->T2->getPoints());
    vtp.addPolygon(tet2->T3->getPoints());
    //
    // vtp.addPolygon(tri1->getPoints());
    // vtp.addPolygon(tri2->getPoints());
    // vtp.addPolygon(tri3->getPoints());
    // std::cout << tri0->getPoints() << std::endl;
    std::ofstream ofs("./output.vtp");
    vtp.write(ofs);
};
#elif defined(case_simple)
int main()
{
    Timer timer;
    auto net = new gn::Network;
    auto p0 = new gn::Point(net, 0., 0., 0.);
    auto p1 = new gn::Point(net, 1., 0., 0.);
    auto p2 = new gn::Point(net, 1., 1., 0.);
    auto p3 = new gn::Point(net, 1., 1., 1.);
    std::cout << "time:" << timer() << std::endl;
    //
    auto e01 = new gn::Edge(net, p0, p1);
    auto e12 = new gn::Edge(net, p1, p2);
    auto e20 = new gn::Edge(net, p2, p0);
    auto e03 = new gn::Edge(net, p0, p3);
    auto e13 = new gn::Edge(net, p1, p3);
    auto e23 = new gn::Edge(net, p2, p3);
    std::cout << "time:" << timer() << std::endl;
    //
    auto tri0 = new gn::Triangle(net, e01, e12, e20);
    auto tri1 = new gn::Triangle(net, e01, e13, e03);
    auto tri2 = new gn::Triangle(net, e20, e03, e23);
    auto tri3 = new gn::Triangle(net, e13, e23, e12);
    std::cout << "time:" << timer() << std::endl;
    //
    auto tet3 = new gn::Tetrahedron(net, tri0, tri1, tri2, tri3);
    std::cout << "time:" << timer() << std::endl;
    /* ------------------------------------------------------ */
    vtkPolygonWriter<gn::Point *> vtp;
    vtp.add({p0, p1, p2, p3});
    vtp.addPolygon(tri0->getPoints());
    vtp.addPolygon(tri1->getPoints());
    vtp.addPolygon(tri2->getPoints());
    vtp.addPolygon(tri3->getPoints());
    std::cout << tri0->getPoints() << std::endl;
    std::ofstream ofs("./output.vtp");
    vtp.write(ofs);
};
#endif