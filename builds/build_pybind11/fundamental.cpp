#define NOMINMAX
#define _CRT_SECURE_NO_WARNINGS
#include "../../include/fundamental.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <tuple>

#include "../../include/InterpolationRBF.hpp"
#include "../../include/Network.hpp"
#include "../../include/rootFinding.hpp"
#include "../../include/fusion.hpp"
/* ------------------------------------------------------ */
PYBIND11_MODULE(fundamental, m)
{
	/* ------------------------------------------------------ */
	m.doc() = "fundamental module";
	/* ------------------------------------------------------ */
	m.def("GaussianQuadratureWeights", &GaussianQuadratureWeights, "");
	m.def("SingularGaussianQuadratureWeights", &SingularGaussianQuadratureWeights, "");
	m.def("GQW", &GaussianQuadratureWeights, "");
	m.def("SGQW", &SingularGaussianQuadratureWeights, "");
	/* ------------------------------------------------------ */
	m.def("Inverse", &Inverse, "");
	/* ------------------------------------------------------ */
	m.def("Subdivide", [](double a, double b, int n)
		  { return Subdivide(a, b, n); });
	/* ------------------------------------------------------ */
	m.def("Flatten", [](const VV_d &U)
		  { return Flatten(U); });
	m.def("Flatten", [](const VV_netPp &U)
		  { return Flatten(U); });
	m.def("Flatten", [](const VV_netFp &U)
		  { return Flatten(U); });
	/* ------------------------------------------------------ */
	m.def("Join", [](const V_d &U, const V_d &V)
		  { return Join(U, V); });
	m.def("Join", [](const VV_d &U, const VV_d &V)
		  { return Join(U, V); });
	/* ------------------------------------------------------ */
	m.def("Times", [](const double V, const Tddd &U)
		  { return U * V; });
	m.def("Times", [](const Tddd &U, const double V)
		  { return U * V; });
	m.def("Times", [](const Tddd &U, const Tddd &V)
		  { return U * V; });
	m.def("Times", [](const V_d &U, const V_d &V)
		  { return U * V; });
	m.def("Times", [](const V_d &U, const double d)
		  { return U * d; });
	m.def("Times", [](const double d, const V_d &U)
		  { return U * d; });
	m.def("Times", [](const Quaternion &Q, const double d)
		  { return Q * d; });
	m.def("Times", [](const double d, const Quaternion &Q)
		  { return Q * d; });
	/* ------------------------------------------------------ */
	m.def("__div__", [](const double d, const Tddd &U)
		  { return d / U; });
	m.def("__div__", [](const Tddd &U, const double d)
		  { return U / d; });
	m.def("__div__", [](const Tddd &V, const Tddd &U)
		  { return V / U; });
	m.def("Divide", [](const double d, const Tddd &U)
		  { return d / U; });
	m.def("Divide", [](const T4d &U, const double d)
		  { return U / d; });
	m.def("Divide", [](const T4d &V, const T4d &U)
		  { return V / U; });
	m.def("Divide", [](const Tddd &U, const double d)
		  { return U / d; });
	m.def("Divide", [](const Tddd &V, const Tddd &U)
		  { return V / U; });
	//
	m.def("Divide", [](const VV_d &U, const double d)
		  { return U / d; });
	m.def("Divide", [](const V_d &U, const double d)
		  { return U / d; });
	m.def("Divide", [](const double d, const V_d &U)
		  { return d / U; });
	/* ------------------------------------------------------ */
	m.def("Dot", [](const V_d &U, const V_d &V)
		  { return Dot(U, V); });
	m.def("Dot", [](const V_d &U, const VV_d &V)
		  { return Dot(U, V); });
	m.def("Dot", [](const VV_d &U, const V_d &V)
		  { return Dot(U, V); });
	m.def("Dot", [](const VV_d &U, const VV_d &V)
		  { return Dot(U, V); });
	/* ------------------------------------------------------ */
	m.def("Cross", [](const V_d &U, const V_d &V)
		  { return Cross(U, V); });
	/* ------------------------------------------------------ */
	m.def("Max", [](const V_d &U)
		  { return Max(U); });
	/* ------------------------------------------------------ */
	m.def("Min", [](const V_d &U)
		  { return Min(U); });
	/* ------------------------------------------------------ */
	m.def("Total", [](const V_d &U)
		  { return std::accumulate(U.cbegin(), U.cend(), 0.); });
	/* ------------------------------------------------------ */
	m.def("Minus", [](const V_d &U)
		  { return -U; });
	/* ------------------------------------------------------ */
	m.def("Subtract", [](const Tddd &U, const Tddd &V)
		  { return U - V; });
	//
	m.def("Subtract", [](const V_d &U, const V_d &V)
		  { return U - V; });
	m.def("Subtract", [](const V_d &U, const VV_d &V)
		  { return U - V; });
	m.def("Subtract", [](const VV_d &U, const V_d &V)
		  { return U - V; });
	m.def("Subtract", [](const VV_d &U, const VV_d &V)
		  { return U - V; });
	m.def("Subtract", [](const V_d &U, const double a)
		  { return U - a; });
	m.def("Subtract", [](const double a, const V_d &U)
		  { return a - U; });
	/* ------------------------------------------------------ */
	m.def("Add", [](const V_d &U, const V_d &V)
		  { return U + V; });
	m.def("Add", [](const V_d &U, const VV_d &V)
		  { return U + V; });
	m.def("Add", [](const VV_d &U, const V_d &V)
		  { return U + V; });
	m.def("Add", [](const VV_d &U, const VV_d &V)
		  { return U + V; });
	m.def("Add", [](const V_d &U, const double a)
		  { return U + a; });
	m.def("Add", [](const double a, const V_d &U)
		  { return U + a; });
	m.def("Add", [](const Quaternion &Q, const double d)
		  { return Q + d; });
	m.def("Add", [](const double d, const Quaternion &Q)
		  { return Q + d; });
	m.def("Add", [](const Quaternion &Q, const Quaternion &QQ)
		  { return Q + QQ; });
	m.def("Add", [](const T4d &d, const Quaternion &Q)
		  { return Q + d; });
	m.def("Add", [](const Quaternion &Q, const T4d &d)
		  { return Q + d; });
	/* ------------------------------------------------------ */
	m.def("Norm", [](const Tdd &U)
		  { return Norm(U); });
	m.def("Norm", [](const T4d &U)
		  { return Norm(U); });
	m.def("Norm", [](const Tddd &U)
		  { return Norm(U); });
	m.def("Norm", [](const Quaternion &Q)
		  { return Norm(Q); });
	m.def("Norm", [](const V_d &V)
		  { return Norm(V); });
	m.def("Normalize", [](const V_d &V)
		  { return Normalize(V); });
	m.def("Normalize", [](const Tddd &V)
		  { return Normalize(V); });
	m.def("Normalize", [](const T4d &V)
		  { return Normalize(V); });
	/* ------------------------------------------------------ */
	m.def("Total", [](const std::vector<Tddd> &V)
		  { return Total(V); });
	m.def("Mean", [](const std::vector<Tddd> &V)
		  { return Mean(V); });
	/* ------------------------------------------------------ */
	m.def("Transpose", [](const VV_d &M)
		  { return Transpose(M); });
	m.def("Transpose", [](const VVV_d &M)
		  { return Transpose(M); });
	m.def("RandomReal", [](const Tdd &minmax)
		  { return RandomReal(minmax); });
	/* ------------------------------------------------------ */
	/* ------------------------------------------------------ */
	pybind11::class_<geometry::Triangle>(m, "Triangle")
		.def(pybind11::init<const T3Tddd &, double>());
	pybind11::class_<geometry::Sphere>(m, "Sphere")
		.def(pybind11::init<const Tddd &, const double>());
	pybind11::class_<geometry::CoordinateBounds>(m, "CoordinateBounds")
		.def(pybind11::init<const geometry::Triangle &>())
		.def(pybind11::init<const geometry::Sphere &>())
		.def_readwrite("bounds", &geometry::CoordinateBounds::bounds);
	/* ------------------------------------------------------ */
	pybind11::class_<InterpolationLagrange<double>>(m, "InterpolationLagrange1d")
		.def(pybind11::init<const std::vector<std::tuple<double, double>> &>())
		.def("__call__", &InterpolationLagrange<double>::operator());
	pybind11::class_<InterpolationVectorRBF>(m, "InterpolationVectorRBF")
		.def(pybind11::init<const VV_d &, const VV_d &>())
		.def("__call__", &InterpolationVectorRBF::operator())
		.def("div", &InterpolationVectorRBF::div)
		.def("grad", &InterpolationVectorRBF::grad)
		.def("N", pybind11::overload_cast<const V_d &>(&InterpolationVectorRBF::N, pybind11::const_))
		.def("N", pybind11::overload_cast<const double, const double>(&InterpolationVectorRBF::N, pybind11::const_))
		.def("gradN", pybind11::overload_cast<const V_d &>(&InterpolationVectorRBF::gradN, pybind11::const_))
		.def("gradN", pybind11::overload_cast<const double, const double>(&InterpolationVectorRBF::gradN, pybind11::const_))
		.def("cross", &InterpolationVectorRBF::cross)
		.def("J", &InterpolationVectorRBF::J);
	/* ------------------------------------------------------ */
	pybind11::class_<Quaternion>(m, "Quaternion")
		.def(pybind11::init<const double, const double, const double, const double>())
		.def(pybind11::init<const T4d &>())
		.def(pybind11::init<const Tddd &, const double>())
		.def("Rs", pybind11::overload_cast<>(&Quaternion::Rs, pybind11::const_))
		.def("Rs", pybind11::overload_cast<const Tddd &>(&Quaternion::Rs, pybind11::const_))
		.def("Rv", pybind11::overload_cast<>(&Quaternion::Rv, pybind11::const_))
		.def("Rv", pybind11::overload_cast<const Tddd &>(&Quaternion::Rv, pybind11::const_))
		.def("R", pybind11::overload_cast<>(&Quaternion::R, pybind11::const_))
		.def("R", pybind11::overload_cast<const Tddd &>(&Quaternion::R, pybind11::const_))
		.def("__call__", &Quaternion::operator())
		.def("__mul__", [](const Quaternion &A, const Quaternion &B)
			 { return A * B; })
		.def("yaw", &Quaternion::yaw)
		.def("pitch", &Quaternion::pitch)
		.def("roll", &Quaternion::roll)
		.def("YPR", &Quaternion::YPR)
		.def("Ryaw", &Quaternion::Ryaw)
		.def("Rpitch", &Quaternion::Rpitch)
		.def("Rroll", &Quaternion::Rroll)
		.def("set", pybind11::overload_cast<const Quaternion &>(&Quaternion::set))
		.def("set", pybind11::overload_cast<const T4d &>(&Quaternion::set))
		.def("d_dt", &Quaternion::d_dt);
	/* ------------------------------------------------------ */
	pybind11::class_<Histogram>(m, "Histogram")
		.def(pybind11::init<const V_d &>())
		.def_readwrite("data", &Histogram::data)
		.def_readwrite("bins", &Histogram::bins)
		.def_readwrite("interval", &Histogram::interval)
		.def_readwrite("mid_interval", &Histogram::mid_interval)
		.def_readwrite("count", &Histogram::count)
		.def_readwrite("cumulative_count", &Histogram::cumulative_count)
		.def_readwrite("bin_width", &Histogram::bin_width)
		.def_readwrite("diff", &Histogram::diff);
	/* ------------------------------------------------------ */
	pybind11::class_<NewtonRaphson<V_d>>(m, "NewtonRaphson")
		.def(pybind11::init<const V_d &>())
		.def("update", &NewtonRaphson<V_d>::update)
		.def_readwrite("X", &NewtonRaphson<V_d>::X)
		.def_readwrite("dX", &NewtonRaphson<V_d>::dX);
	pybind11::class_<Fusion>(m, "Fusion")
		.def(pybind11::init<const Tddd &, const Tddd &>())
		.def("solveForQuaternion", pybind11::overload_cast<const Tddd &, const Tddd &, const Tddd &, const double>(&Fusion::solveForQuaternion))
		.def("solveForQuaternionModified", pybind11::overload_cast<const Tddd &, const Tddd &, const Tddd &, const double>(&Fusion::solveForQuaternionModified))
		.def("history", &Fusion::history)
		.def("historyQ", &Fusion::historyQ)
		.def("historyA", &Fusion::historyA)
		.def("historyM", &Fusion::historyM)
		.def("historyW", &Fusion::historyW)
		.def("historyAbody", &Fusion::historyAbody)
		.def("interpA", &Fusion::interpA)
		.def("interpM", &Fusion::interpM)
		.def("interpW", &Fusion::interpW)
		.def("interpYPR", &Fusion::interpYPR)
		.def("interpAbody", &Fusion::interpAbody)
		.def("interp", &Fusion::interp)
		.def("calculateOffsetM", &Fusion::calculateOffsetM)
		.def("setMagTransMat", &Fusion::setMagTransMat)
		.def("setOffsetM", &Fusion::setOffsetM)
		.def("setScaleM", &Fusion::setScaleM);

	/* ------------------------------------------------------ */
	// https://pybind11.readthedocs.io/en/stable/classes.html#overloaded-methods
	// pybind11::const_を忘れない
	// pybind11::class_<Network>(m, "Network")
	// 	.def(pybind11::init<const std::string &>())
	// 	.def("getPoints", pybind11::overload_cast<>(&Network::getPoints, pybind11::const_))
	// 	.def("getFaces", &Network::getFaces)
	// 	.def("getPointIndex", pybind11::overload_cast<const VV_netPp &>(&Network::getPointIndex, pybind11::const_))
	// 	.def("getPointIndex", pybind11::overload_cast<const V_netPp &>(&Network::getPointIndex, pybind11::const_))
	// 	.def("getPointIndex", pybind11::overload_cast<const netPp>(&Network::getPointIndex, pybind11::const_))
	// 	.def("getFaceIndex", pybind11::overload_cast<const VV_netFp &>(&Network::getFaceIndex, pybind11::const_))
	// 	.def("getFaceIndex", pybind11::overload_cast<const V_netFp &>(&Network::getFaceIndex, pybind11::const_))
	// 	.def("getFaceIndex", pybind11::overload_cast<const netFp>(&Network::getFaceIndex, pybind11::const_));

	pybind11::class_<networkPoint>(m, "networkPoint")
		.def(pybind11::init<Network *, Network *, const Tddd &>());
	// .def("getNeighborsPolarAsTuple", &networkPoint::getNeighborsPolarAsTuple)
	// .def("getNeighborsPolarAsVariant", &networkPoint::getNeighborsPolarAsVariant);

	pybind11::class_<networkFace>(m, "networkFace")
		.def(pybind11::init<Network *, Network *, const V_netLp &>());
	// .def("getPoints", pybind11::overload_cast<const netPp>(&networkFace::getPoints, pybind11::const_));

	pybind11::class_<networkLine>(m, "networkLine")
		.def(pybind11::init<Network *, netP *, netP *>());

	m.def("extractX", [](const V_netPp &ps)
		  { return extractX(ps); });
	m.def("extractX", [](const V_netFp &fs)
		  { return extractX(fs); });
	m.def("extractPoints", [](const V_netFp &fs)
		  { return extractPoints(fs); });
}
