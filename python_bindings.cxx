//
// python_bindings.cc -- python bindings for Dirac++ user classes
//                       using the Boost.Python C++ interface.
//
// author: richard.t.jones at uconn.edu
// version: (still under construction)

#include <boost/python.hpp>

#include <Double.h>
#include <Complex.h>
#include <TCrossSection.h>
#include <TDiracMatrix.h>
#include <TDiracSpinor.h>
#include <TFourVectorComplex.h>
#include <TFourVectorReal.h>
#include <TLepton.h>
#include <TLorentzBoost.h>
#include <TLorentzTransform.h>
#include <TPauliMatrix.h>
#include <TPauliSpinor.h>
#include <TPhoton.h>
#include <TThreeRotation.h>
#include <TThreeVectorComplex.h>
#include <TThreeVectorReal.h>
#include <constants.h>

Complex_t Complex_abs(const Complex_t val) {
   return std::abs(val);
}

Complex_t Complex_arg(const Complex_t val) {
   return std::arg(val);
}

Complex_t Complex_norm(const Complex_t val) {
   return std::norm(val);
}

// Wrap overloaded methods using functions with unique names.

void (TThreeVectorReal::*TThreeVectorReal_GetCartesian3)(LDouble_t &x, LDouble_t &y, LDouble_t &z) const =
     &TThreeVectorReal::GetCartesian;
void (TThreeVectorReal::*TThreeVectorReal_GetCartesian1)(LDouble_t *array) const =
     &TThreeVectorReal::GetCartesian;

LDouble_t (TThreeVectorReal::*TThreeVectorReal_DistanceTo3)(LDouble_t x, LDouble_t y, LDouble_t z) const =
     &TThreeVectorReal::DistanceTo;
LDouble_t (TThreeVectorReal::*TThreeVectorReal_DistanceTo1)(const LDouble_t *array) const =
     &TThreeVectorReal::DistanceTo;
LDouble_t (TThreeVectorReal::*TThreeVectorReal_DistanceTo)(const TThreeVectorReal &vec2) const =
     &TThreeVectorReal::DistanceTo;

TThreeVectorReal &(TThreeVectorReal::*TThreeVectorReal_Rotate1)(const TThreeRotation &rotOp) =
     &TThreeVectorReal::Rotate;
TThreeVectorReal &(TThreeVectorReal::*TThreeVectorReal_Rotate2)(const TUnitVector &ahat, LDouble_t angle) =
     &TThreeVectorReal::Rotate;
TThreeVectorReal &(TThreeVectorReal::*TThreeVectorReal_Rotate3)(LDouble_t phi, LDouble_t theta, LDouble_t psi) =
     &TThreeVectorReal::Rotate;

TThreeVectorReal &(TThreeVectorReal::*TThreeVectorReal_Cross)(const TThreeVectorReal &other) =
     &TThreeVectorReal::Cross;
TThreeVectorReal &(TThreeVectorReal::*TThreeVectorReal_Cross2)(const TThreeVectorReal &va, const TThreeVectorReal &vb) =
     &TThreeVectorReal::Cross;

void TThreeVectorReal_Print(TThreeVectorReal &obj) {
   obj.Print();
}

LDouble_t TThreeVectorReal_getitem(const TThreeVectorReal &obj, Int_t index) {
   return obj[index];
}

// Create a python module containing all of the user classes
// that are needed to interact with Dirac++ objects from python.
// Here it is named libdiracxx (happens to also be the name of
// the shared library) but when it is loaded from python, the
// name of the module is diracxx.

BOOST_PYTHON_MODULE(libDirac)
{
   boost::python::enum_<EPauliIndex>("EPauliIndex")
      .value("kPauliOne", kPauliOne)
      .value("kPauliSigma1", kPauliSigma1)
      .value("kPauliSigma2", kPauliSigma2)
      .value("kPauliSigma3", kPauliSigma3)
   ;
   boost::python::enum_<EDiracIndex>("EDiracIndex")
      .value("kDiracOne", kDiracOne)
      .value("kDiracGamma0", kDiracGamma0)
      .value("kDiracGamma1", kDiracGamma1)
      .value("kDiracGamma2", kDiracGamma2)
      .value("kDiracGamma3", kDiracGamma3)
      .value("kDiracGamma4", kDiracGamma4)
      .value("kDiracGamma5", kDiracGamma5)
      .value("kDiracSigma1", kDiracSigma1)
      .value("kDiracSigma2", kDiracSigma2)
      .value("kDiracSigma3", kDiracSigma3)
      .value("kDiracKappa1", kDiracKappa1)
      .value("kDiracKappa2", kDiracKappa2)
      .value("kDiracKappa3", kDiracKappa3)
   ;

   boost::python::class_<TThreeVectorReal, TThreeVectorReal*>
         ("TThreeVectorReal",
          "three vector with real components")
      .def(boost::python::init<const LDouble_t, const LDouble_t, const LDouble_t>())
      .def(boost::python::init<const Float_t *>())
      .def(boost::python::init<const LDouble_t *>())
      .def(boost::python::init<const TThreeVectorReal &>())
      .def("__getitem__", &TThreeVectorReal_getitem)
      .def("SetResolution", &TThreeVectorReal::SetResolution)
      .def("Resolution", &TThreeVectorReal::Resolution)
      .def("Length", &TThreeVectorReal::Length)
      .def("LengthSqr", &TThreeVectorReal::LengthSqr)
      .def("Rho", &TThreeVectorReal::Rho)
      .def("RhoSqr", &TThreeVectorReal::RhoSqr)
      .def("Theta", &TThreeVectorReal::Theta)
      .def("CosTheta", &TThreeVectorReal::CosTheta)
      .def("Phi", &TThreeVectorReal::Phi)
      .def("GetPolar", &TThreeVectorReal::GetPolar)
      .def("GetCartesian", TThreeVectorReal_GetCartesian3)
      .def("GetCartesian", TThreeVectorReal_GetCartesian1)
      .def("DistanceTo", TThreeVectorReal_DistanceTo3)
      .def("DistanceTo", TThreeVectorReal_DistanceTo1)
      .def("DistanceTo", TThreeVectorReal_DistanceTo)
      .def(boost::python::self_ns::self += TThreeVectorReal())
      .def(boost::python::self_ns::self -= TThreeVectorReal())
      .def(boost::python::self_ns::self *= LDouble_t())
      .def(boost::python::self_ns::self /= LDouble_t())
      .def(boost::python::self_ns::self + TThreeVectorReal())
      .def(boost::python::self_ns::self - TThreeVectorReal())
      .def(boost::python::self_ns::self * LDouble_t())
      .def(boost::python::self_ns::self / LDouble_t())
      .def(TThreeVectorReal() + boost::python::self_ns::self)
      .def(TThreeVectorReal() - boost::python::self_ns::self)
      .def(LDouble_t() * boost::python::self_ns::self)
      .def("__eq__", &TThreeVectorReal::operator==)
      .def("__ne__", &TThreeVectorReal::operator!=)
      .def("Zero", &TThreeVectorReal::Zero,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SpaceInv", &TThreeVectorReal::SpaceInv,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Normalize", &TThreeVectorReal::Normalize,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetPolar", &TThreeVectorReal::SetPolar,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TThreeVectorReal_Rotate1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TThreeVectorReal_Rotate3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TThreeVectorReal_Rotate2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Cross", TThreeVectorReal_Cross,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Cross", TThreeVectorReal_Cross2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Dot", &TThreeVectorReal::Dot)
      .def("__neg__", &TThreeVectorReal::operator-)
      .def("Print", &TThreeVectorReal::Print)
      .def("Print", &TThreeVectorReal_Print)
   ;
}
