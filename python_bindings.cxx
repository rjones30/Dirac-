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


void (TThreeVectorComplex::*TThreeVectorComplex_GetCartesian3)(Complex_t &x, Complex_t &y, Complex_t &z) const =
     &TThreeVectorComplex::GetCartesian;
void (TThreeVectorComplex::*TThreeVectorComplex_GetCartesian1)(Complex_t *array) const =
     &TThreeVectorComplex::GetCartesian;

LDouble_t (TThreeVectorComplex::*TThreeVectorComplex_DistanceTo3)(Complex_t x, Complex_t y, Complex_t z) const =
     &TThreeVectorComplex::DistanceTo;
LDouble_t (TThreeVectorComplex::*TThreeVectorComplex_DistanceTo1)(const Complex_t *array) const =
     &TThreeVectorComplex::DistanceTo;
LDouble_t (TThreeVectorComplex::*TThreeVectorComplex_DistanceTo)(const TThreeVectorComplex &vec2) const =
     &TThreeVectorComplex::DistanceTo;

TThreeVectorComplex &(TThreeVectorComplex::*TThreeVectorComplex_Rotate1)(const TThreeRotation &rotOp) =
     &TThreeVectorComplex::Rotate;
TThreeVectorComplex &(TThreeVectorComplex::*TThreeVectorComplex_Rotate2)(const TUnitVector &ahat, LDouble_t angle) =
     &TThreeVectorComplex::Rotate;
TThreeVectorComplex &(TThreeVectorComplex::*TThreeVectorComplex_Rotate3)(LDouble_t phi, LDouble_t theta, LDouble_t psi) =
     &TThreeVectorComplex::Rotate;

TThreeVectorComplex &(TThreeVectorComplex::*TThreeVectorComplex_Cross)(const TThreeVectorComplex &other) =
     &TThreeVectorComplex::Cross;
TThreeVectorComplex &(TThreeVectorComplex::*TThreeVectorComplex_Cross2)(const TThreeVectorComplex &va, const TThreeVectorComplex &vb) =
     &TThreeVectorComplex::Cross;

void TThreeVectorComplex_Print(TThreeVectorComplex &obj) {
   obj.Print();
}

Complex_t TThreeVectorComplex_getitem(const TThreeVectorComplex &obj, Int_t index) {
   return obj[index];
}


void (TFourVectorReal::*TFourVectorReal_GetCoord4)(LDouble_t &t, LDouble_t &x, LDouble_t &y, LDouble_t &z) const =
     &TFourVectorReal::GetCoord;
void (TFourVectorReal::*TFourVectorReal_GetCoord1)(LDouble_t *array) const =
     &TFourVectorReal::GetCoord;

LDouble_t (TFourVectorReal::*TFourVectorReal_DistanceTo4)(LDouble_t t, LDouble_t x, LDouble_t y, LDouble_t z) const =
     &TFourVectorReal::DistanceTo;
LDouble_t (TFourVectorReal::*TFourVectorReal_DistanceTo1)(const LDouble_t *array) const =
     &TFourVectorReal::DistanceTo;
LDouble_t (TFourVectorReal::*TFourVectorReal_DistanceTo)(const TFourVectorReal &vec2) const =
     &TFourVectorReal::DistanceTo;

TFourVectorReal &(TFourVectorReal::*TFourVectorReal_Boost)(const TLorentzBoost &boostOp) =
     &TFourVectorReal::Boost;
TFourVectorReal &(TFourVectorReal::*TFourVectorReal_Boost4)(LDouble_t betaX, LDouble_t betaY, LDouble_t betaZ) =
     &TFourVectorReal::Boost;
TFourVectorReal &(TFourVectorReal::*TFourVectorReal_Boost1)(const LDouble_t *beta) =
     &TFourVectorReal::Boost;
TFourVectorReal &(TFourVectorReal::*TFourVectorReal_Boost3)(const TThreeVectorReal &beta) =
     &TFourVectorReal::Boost;
TFourVectorReal &(TFourVectorReal::*TFourVectorReal_Boost2)(const TUnitVector &bhat, LDouble_t beta) =
     &TFourVectorReal::Boost;

void TFourVectorReal_Print(TFourVectorReal &obj) {
   obj.Print();
}

LDouble_t TFourVectorReal_getitem(const TFourVectorReal &obj, Int_t index) {
   return obj[index];
}


void (TFourVectorComplex::*TFourVectorComplex_GetCoord4)(Complex_t &t, Complex_t &x, Complex_t &y, Complex_t &z) const =
     &TFourVectorComplex::GetCoord;
void (TFourVectorComplex::*TFourVectorComplex_GetCoord1)(Complex_t *array) const =
     &TFourVectorComplex::GetCoord;

LDouble_t (TFourVectorComplex::*TFourVectorComplex_DistanceTo4)(Complex_t t, Complex_t x, Complex_t y, Complex_t z) const =
     &TFourVectorComplex::DistanceTo;
LDouble_t (TFourVectorComplex::*TFourVectorComplex_DistanceTo1)(const Complex_t *array) const =
     &TFourVectorComplex::DistanceTo;
LDouble_t (TFourVectorComplex::*TFourVectorComplex_DistanceTo2)(const LDouble_t *array) const =
     &TFourVectorComplex::DistanceTo;
LDouble_t (TFourVectorComplex::*TFourVectorComplex_DistanceTo3)(const Float_t *array) const =
     &TFourVectorComplex::DistanceTo;
LDouble_t (TFourVectorComplex::*TFourVectorComplex_DistanceTo)(const TFourVectorComplex &vec2) const =
     &TFourVectorComplex::DistanceTo;

TFourVectorComplex &(TFourVectorComplex::*TFourVectorComplex_Boost)(const TLorentzBoost &boostOp) =
     &TFourVectorComplex::Boost;
TFourVectorComplex &(TFourVectorComplex::*TFourVectorComplex_Boost4)(LDouble_t betaX, LDouble_t betaY, LDouble_t betaZ) =
     &TFourVectorComplex::Boost;
TFourVectorComplex &(TFourVectorComplex::*TFourVectorComplex_Boost1)(const LDouble_t *beta) =
     &TFourVectorComplex::Boost;
TFourVectorComplex &(TFourVectorComplex::*TFourVectorComplex_Boost2)(const TUnitVector &bhat, LDouble_t beta) =
     &TFourVectorComplex::Boost;
TFourVectorComplex &(TFourVectorComplex::*TFourVectorComplex_Boost3)(const TThreeVectorReal &beta) =
     &TFourVectorComplex::Boost;

Complex_t (TFourVectorComplex::*TFourVectorComplex_ScalarProd)(const TFourVectorComplex &other) =
     &TFourVectorComplex::ScalarProd;
Complex_t (TFourVectorComplex::*TFourVectorComplex_ScalarProd2)(const TFourVectorComplex &v1, const TFourVectorComplex &v2) =
     &TFourVectorComplex::ScalarProd;

void TFourVectorComplex_Print(TFourVectorComplex &obj) {
   obj.Print();
}

Complex_t TFourVectorComplex_getitem(const TFourVectorComplex &obj, Int_t index) {
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

   boost::python::class_<TThreeVectorComplex, TThreeVectorComplex*>
         ("TThreeVectorComplex",
          "three vector with complex components")
      .def(boost::python::init<const Complex_t, const Complex_t, const Complex_t>())
      .def(boost::python::init<const Float_t *>())
      .def(boost::python::init<const LDouble_t *>())
      .def(boost::python::init<const Complex_t *>())
      .def(boost::python::init<const TThreeVectorReal &>())
      .def(boost::python::init<const TThreeVectorComplex &>())
      .def("__getitem__", &TThreeVectorComplex_getitem)
      .def("SetResolution", &TThreeVectorComplex::SetResolution)
      .def("Resolution", &TThreeVectorComplex::Resolution)
      .def("Length", &TThreeVectorComplex::Length)
      .def("LengthSqr", &TThreeVectorComplex::LengthSqr)
      .def("RealPart", &TThreeVectorComplex::RealPart)
      .def("ImagPart", &TThreeVectorComplex::ImagPart)
      .def("GetCartesian", TThreeVectorComplex_GetCartesian3)
      .def("GetCartesian", TThreeVectorComplex_GetCartesian1)
      .def("DistanceTo", TThreeVectorComplex_DistanceTo3)
      .def("DistanceTo", TThreeVectorComplex_DistanceTo1)
      .def("DistanceTo", TThreeVectorComplex_DistanceTo)
      .def(boost::python::self_ns::self += TThreeVectorComplex())
      .def(boost::python::self_ns::self -= TThreeVectorComplex())
      .def(boost::python::self_ns::self *= Complex_t())
      .def(boost::python::self_ns::self *= LDouble_t())
      .def(boost::python::self_ns::self /= Complex_t())
      .def(boost::python::self_ns::self /= LDouble_t())
      .def(boost::python::self_ns::self + TThreeVectorComplex())
      .def(boost::python::self_ns::self - TThreeVectorComplex())
      .def(boost::python::self_ns::self * Complex_t())
      .def(boost::python::self_ns::self / Complex_t())
      .def(TThreeVectorComplex() + boost::python::self_ns::self)
      .def(TThreeVectorComplex() - boost::python::self_ns::self)
      .def(Complex_t() * boost::python::self_ns::self)
      .def("__eq__", &TThreeVectorComplex::operator==)
      .def("__ne__", &TThreeVectorComplex::operator!=)
      .def("Zero", &TThreeVectorComplex::Zero,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Conj", &TThreeVectorComplex::Conj,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SpaceInv", &TThreeVectorComplex::SpaceInv,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Normalize", &TThreeVectorComplex::Normalize,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TThreeVectorComplex_Rotate1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TThreeVectorComplex_Rotate3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TThreeVectorComplex_Rotate2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Cross", TThreeVectorComplex_Cross,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Cross", TThreeVectorComplex_Cross2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Dot", &TThreeVectorComplex::Dot)
      .def("__neg__", &TThreeVectorComplex::operator-)
      .def("Print", &TThreeVectorComplex::Print)
      .def("Print", &TThreeVectorComplex_Print)
   ;

   boost::python::class_<TFourVectorReal, TFourVectorReal*,
          boost::python::bases<TThreeVectorReal> >
         ("TFourVectorReal",
          "four vector with real components")
      .def(boost::python::init<const LDouble_t, const LDouble_t, const LDouble_t, const LDouble_t>())
      .def(boost::python::init<const Float_t *>())
      .def(boost::python::init<const LDouble_t *>())
      .def(boost::python::init<const LDouble_t, const TThreeVectorReal &>())
      .def(boost::python::init<const TFourVectorReal &>())
      .def("__getitem__", &TFourVectorReal_getitem)
      .def("Resolution", &TFourVectorReal::Resolution)
      .def("Invariant", &TFourVectorReal::Invariant)
      .def("InvariantSqr", &TFourVectorReal::InvariantSqr)
      .def("GetCoord", TFourVectorReal_GetCoord4)
      .def("GetCoord", TFourVectorReal_GetCoord1)
      .def("DistanceTo", TFourVectorReal_DistanceTo4)
      .def("DistanceTo", TFourVectorReal_DistanceTo1)
      .def("DistanceTo", TFourVectorReal_DistanceTo)
      .def(boost::python::self_ns::self += TFourVectorReal())
      .def(boost::python::self_ns::self -= TFourVectorReal())
      .def(boost::python::self_ns::self *= LDouble_t())
      .def(boost::python::self_ns::self /= LDouble_t())
      .def(boost::python::self_ns::self + TFourVectorReal())
      .def(boost::python::self_ns::self - TFourVectorReal())
      .def(boost::python::self_ns::self * LDouble_t())
      .def(boost::python::self_ns::self / LDouble_t())
      .def(TFourVectorReal() + boost::python::self_ns::self)
      .def(TFourVectorReal() - boost::python::self_ns::self)
      .def(LDouble_t() * boost::python::self_ns::self)
      .def("__eq__", &TFourVectorReal::operator==)
      .def("__ne__", &TFourVectorReal::operator!=)
      .def("Zero", &TFourVectorReal::Zero,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Transform", &TFourVectorReal::Transform,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorReal_Boost,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorReal_Boost1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorReal_Boost2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorReal_Boost3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorReal_Boost4,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("BoostToRest", &TFourVectorReal::BoostToRest,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("BoostFromRest", &TFourVectorReal::BoostFromRest,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("ScalarProd", &TFourVectorReal::ScalarProd)
      .def("__neg__", &TFourVectorReal::operator-)
      .def("Print", &TFourVectorReal::Print)
      .def("Print", &TFourVectorReal_Print)
   ;

   boost::python::class_<TFourVectorComplex, TFourVectorComplex*>
         ("TFourVectorComplex",
          "four vector with complex components")
      .def(boost::python::init<const Complex_t &, const Complex_t &, const Complex_t &, const Complex_t &>())
      .def(boost::python::init<const Float_t *>())
      .def(boost::python::init<const LDouble_t *>())
      .def(boost::python::init<const Complex_t *>())
      .def(boost::python::init<const Complex_t &, const TThreeVectorComplex &>())
      .def(boost::python::init<const TFourVectorReal &>())
      .def(boost::python::init<const TFourVectorComplex &>())
      .def("__getitem__", &TFourVectorComplex_getitem)
      .def("SetResolution", &TFourVectorComplex::SetResolution)
      .def("Resolution", &TFourVectorComplex::Resolution)
      .def("Invariant", &TFourVectorComplex::Invariant)
      .def("InvariantSqr", &TFourVectorComplex::InvariantSqr)
      .def("RealPart", &TFourVectorComplex::RealPart)
      .def("ImagPart", &TFourVectorComplex::ImagPart)
      .def("GetCoord", TFourVectorComplex_GetCoord4)
      .def("GetCoord", TFourVectorComplex_GetCoord1)
      .def("DistanceTo", TFourVectorComplex_DistanceTo4)
      .def("DistanceTo", TFourVectorComplex_DistanceTo3)
      .def("DistanceTo", TFourVectorComplex_DistanceTo2)
      .def("DistanceTo", TFourVectorComplex_DistanceTo1)
      .def("DistanceTo", TFourVectorComplex_DistanceTo)
      .def(boost::python::self_ns::self += TFourVectorComplex())
      .def(boost::python::self_ns::self += TFourVectorReal())
      .def(boost::python::self_ns::self -= TFourVectorComplex())
      .def(boost::python::self_ns::self -= TFourVectorReal())
      .def(boost::python::self_ns::self *= Complex_t())
      .def(boost::python::self_ns::self /= Complex_t())
      .def(boost::python::self_ns::self + TFourVectorComplex())
      .def(boost::python::self_ns::self - TFourVectorComplex())
      .def(boost::python::self_ns::self * Complex_t())
      .def(boost::python::self_ns::self / Complex_t())
      .def(TFourVectorComplex() + boost::python::self_ns::self)
      .def(TFourVectorComplex() - boost::python::self_ns::self)
      .def(Complex_t() * boost::python::self_ns::self)
      .def("__eq__", &TFourVectorComplex::operator==)
      .def("__ne__", &TFourVectorComplex::operator!=)
      .def("Zero", &TFourVectorComplex::Zero,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Conj", &TFourVectorComplex::Conj,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Transform", &TFourVectorComplex::Transform,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorComplex_Boost1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorComplex_Boost2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorComplex_Boost3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost",TFourVectorComplex_Boost4,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TFourVectorComplex_Boost,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("BoostFromRest", &TFourVectorComplex::BoostFromRest,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("BoostToRest", &TFourVectorComplex::BoostToRest,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("ScalarProd", TFourVectorComplex_ScalarProd)
      .def("ScalarProd", TFourVectorComplex_ScalarProd2)
      .def("__neg__", &TFourVectorComplex::operator-)
      .def("Print", &TFourVectorComplex::Print)
      .def("Print", &TFourVectorComplex_Print)
   ;
}
