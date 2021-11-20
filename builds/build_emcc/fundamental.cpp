/**
 * emscriptenは，pybind11と違って，std::vectorとの親和性がよくない．
 * doubleやintの値をjavascriptとc++でやりとりにするときに使うことにしよう
 */

#define NOMINMAX
#define _CRT_SECURE_NO_WARNINGS
#include "../../include/fundamental.hpp"

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>

#include <tuple>

#include "../../include/InterpolationRBF.hpp"
#include "../../include/Network.hpp"
#include "../../include/rootFinding.hpp"
/* ------------------------------------------------------ */
// #ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
// #endif
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
// V_d operator*(const Quaternion &A, const Quaternion &B)
// {
//   V_d v = -Dot(A.v, B.v) + A.a * B.v + B.a * A.v + Cross(A.v, B.v);
//   return {(A.a) * (B.a), v[0], v[1], v[2]};
// };
/* ------------------------------------------------------ */
float divide(float x, float y)
{
  return x / y;
};
// #ifdef __EMSCRIPTEN__
// EMSCRIPTEN_BINDINGS(stl_wrappers)
EMSCRIPTEN_BINDINGS(stl_wrappers)
{
  // ベクトルのバインドを登録する必要がある．
  emscripten::register_vector<int>("V_i");
  emscripten::register_vector<double>("V_d");
  emscripten::register_vector<float>("V_f");
  emscripten::register_vector<std::vector<double>>("VV_d");
  emscripten::register_vector<std::vector<std::vector<double>>>("VVV_d");
  emscripten::register_vector<std::vector<int>>("VV_i");
  emscripten::register_vector<std::vector<float>>("VV_f");
}

EMSCRIPTEN_BINDINGS(fundamental)
{
  /* ------------------------------------------------------ */
  // m.doc() = "fundamental module";
  // /* ------------------------------------------------------ */
  // m.def("GaussianQuadratureWeights", &GaussianQuadratureWeights, "");
  // m.def("SingularGaussianQuadratureWeights", &SingularGaussianQuadratureWeights, "");
  // m.def("GQW", &GaussianQuadratureWeights, "");
  // m.def("SGQW", &SingularGaussianQuadratureWeights, "");
  // /* ------------------------------------------------------ */
  // m.def("Inverse", &Inverse, "");
  /* ------------------------------------------------------ */
  emscripten::function("divide", divide);
  emscripten::function("Flatten", emscripten::select_overload<V_d(const VV_d &)>(&Flatten));
  emscripten::function("Flatten", emscripten::select_overload<VV_d(const VVV_d &)>(&Flatten));
  // emscripten::function("Flatten", [](const VV_d &U)
  //                      { return Flatten(U); });
  // emscripten::function("Flatten", [](const VV_d &U)
  //                      { return Flatten(U); });
  // m.def("Flatten", [](const VV_netPp &U)
  //       { return Flatten(U); });
  // m.def("Flatten", [](const VV_netFp &U)
  //       { return Flatten(U); });
  // /* ------------------------------------------------------ */
  // m.def("Join", [](const V_d &U, const V_d &V)
  //       { return Join(U, V); });
  // m.def("Join", [](const VV_d &U, const VV_d &V)
  //       { return Join(U, V); });
  // /* ------------------------------------------------------ */
  // m.def("Dot", [](const V_d &U, const V_d &V)
  //       { return Dot(U, V); });
  // m.def("Dot", [](const V_d &U, const VV_d &V)
  //       { return Dot(U, V); });
  // m.def("Dot", [](const VV_d &U, const V_d &V)
  //       { return Dot(U, V); });
  // m.def("Dot", [](const VV_d &U, const VV_d &V)
  //       { return Dot(U, V); });
  // /* ------------------------------------------------------ */
  // m.def("Max", [](const V_d &U)
  //       { return Max(U); });
  // /* ------------------------------------------------------ */
  // m.def("Min", [](const V_d &U)
  //       { return Min(U); });
  // /* ------------------------------------------------------ */
  // m.def("Total", [](const V_d &U)
  //       { return std::accumulate(U.cbegin(), U.cend(), 0.); });
  // /* ------------------------------------------------------ */
  // m.def("Minus", [](const V_d &U)
  //       { return -U; });
  // /* ------------------------------------------------------ */
  // m.def("Subtract", [](const V_d &U, const V_d &V)
  //       { return U - V; });
  // m.def("Subtract", [](const V_d &U, const double a)
  //       { return U - a; });
  // m.def("Subtract", [](const double a, const V_d &U)
  //       { return a - U; });
  // /* ------------------------------------------------------ */
  // m.def("Add", [](const V_d &U, const V_d &V)
  //       { return U + V; });
  // m.def("Add", [](const V_d &U, const double a)
  //       { return U + a; });
  // m.def("Add", [](const double a, const V_d &U)
  //       { return U + a; });
  // /* ------------------------------------------------------ */
  // m.def("Norm", [](const V_d &V)
  //       { return Norm(V); });
  // m.def("Normalize", [](const V_d &V)
  //       { return Normalize(V); });
  // /* ------------------------------------------------------ */
  // m.def("Transpose", [](const VV_d &M)
  //       { return Transpose(M); });
  // m.def("Transpose", [](const VVV_d &M)
  //       { return Transpose(M); });

  // m.def("RandomReal", [](const V_d &minmax)
  //       { return RandomReal(minmax); });
  // /* ------------------------------------------------------ */
  // pybind11::class_<InterpolationVectorRBF>(m, "InterpolationVectorRBF")
  //     .def(pybind11::init<const VV_d &, const VV_d &>())
  //     .def("__call__", &InterpolationVectorRBF::operator())
  //     .def("div", &InterpolationVectorRBF::div)
  //     .def("N", pybind11::overload_cast<const V_d &>(&InterpolationVectorRBF::N, pybind11::const_))
  //     .def("N", pybind11::overload_cast<const double, const double>(&InterpolationVectorRBF::N, pybind11::const_))
  //     .def("gradN", pybind11::overload_cast<const V_d &>(&InterpolationVectorRBF::gradN, pybind11::const_))
  //     .def("gradN", pybind11::overload_cast<const double, const double>(&InterpolationVectorRBF::gradN, pybind11::const_))
  //     .def("cross", &InterpolationVectorRBF::cross)
  //     .def("J", &InterpolationVectorRBF::J);
  // /* ------------------------------------------------------ */
  // pybind11::class_<Quaternion>(m, "Quaternion")
  //     .def(pybind11::init<const double, const double, const double, const double>())
  //     .def(pybind11::init<const V_d &>())
  //     .def("R", pybind11::overload_cast<>(&Quaternion::R, pybind11::const_))
  //     .def("R", pybind11::overload_cast<const V_d &>(&Quaternion::R, pybind11::const_))
  //     .def("__call__", &Quaternion::operator())
  //     .def("__mul__", [](const Quaternion &A, const Quaternion &B)
  //          { return A * B; })
  //     .def("yaw", &Quaternion::yaw)
  //     .def("pitch", &Quaternion::pitch)
  //     .def("roll", &Quaternion::roll)
  //     .def("YPR", &Quaternion::YPR)
  //     .def("YRP", &Quaternion::YRP)
  //     .def("set", &Quaternion::set);
  // /* ------------------------------------------------------ */
  // pybind11::class_<NewtonRaphson>(m, "NewtonRaphson")
  //     .def(pybind11::init<const V_d &>())
  //     .def("update", &NewtonRaphson::update)
  //     .def_readwrite("X", &NewtonRaphson::X)
  //     .def_readwrite("dX", &NewtonRaphson::dX);
  // /* ------------------------------------------------------ */
  // // https://pybind11.readthedocs.io/en/stable/classes.html#overloaded-methods
  // // pybind11::const_を忘れない
  // pybind11::class_<NetworkObj>(m, "NetworkObj")
  //     .def(pybind11::init<const std::string &>())
  //     .def("getPoints", pybind11::overload_cast<>(&NetworkObj::getPoints, pybind11::const_))
  //     .def("getFaces", &NetworkObj::getFaces)
  //     .def("getPointIndex", pybind11::overload_cast<const VV_netPp &>(&NetworkObj::getPointIndex, pybind11::const_))
  //     .def("getPointIndex", pybind11::overload_cast<const V_netPp &>(&NetworkObj::getPointIndex, pybind11::const_))
  //     .def("getPointIndex", pybind11::overload_cast<const netPp>(&NetworkObj::getPointIndex, pybind11::const_))
  //     .def("getFaceIndex", pybind11::overload_cast<const VV_netFp &>(&NetworkObj::getFaceIndex, pybind11::const_))
  //     .def("getFaceIndex", pybind11::overload_cast<const V_netFp &>(&NetworkObj::getFaceIndex, pybind11::const_))
  //     .def("getFaceIndex", pybind11::overload_cast<const netFp>(&NetworkObj::getFaceIndex, pybind11::const_));

  // pybind11::class_<networkPoint>(m, "networkPoint")
  //     .def(pybind11::init<Network *, Network *, const V_d &>());
  // // .def("getNeighborsPolarAsTuple", &networkPoint::getNeighborsPolarAsTuple)
  // // .def("getNeighborsPolarAsVariant", &networkPoint::getNeighborsPolarAsVariant);

  // pybind11::class_<networkFace>(m, "networkFace")
  //     .def(pybind11::init<Network *, Network *, const V_netLp &>());
  // // .def("getPoints", pybind11::overload_cast<const netPp>(&networkFace::getPoints, pybind11::const_));

  // pybind11::class_<networkLine>(m, "networkLine")
  //     .def(pybind11::init<Network *, netP *, netP *>());

  // m.def("extractX", [](const V_netPp &ps)
  //       { return extractX(ps); });
  // m.def("extractX", [](const V_netFp &fs)
  //       { return extractX(fs); });
  // m.def("extractPoints", [](const V_netFp &fs)
  //       { return extractPoints(fs); });
}

// #endif