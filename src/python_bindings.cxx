//
// python_bindings.cc -- python bindings for Dirac++ user classes
//                       using the Boost.Python C++ interface.
//
// author: richard.t.jones at uconn.edu
// version: (still under construction)

#define BOOST_BIND_GLOBAL_PLACEHOLDERS 1
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


TThreeRotation &(TThreeRotation::*TThreeRotation_SetAxis)(const TThreeVectorReal &axis) =
   &TThreeRotation::SetAxis;
TThreeRotation &(TThreeRotation::*TThreeRotation_SetAxis2)(const TUnitVector &ahat, LDouble_t angle) =
   &TThreeRotation::SetAxis;

void TThreeRotation_Print(TThreeRotation &obj) {
   obj.Print();
}


TLorentzBoost &(TLorentzBoost::*TLorentzBoost_SetBeta)(LDouble_t betaX, LDouble_t betaY, LDouble_t betaZ) =
   &TLorentzBoost::SetBeta;
TLorentzBoost &(TLorentzBoost::*TLorentzBoost_SetBeta1)(const LDouble_t *beta) =
   &TLorentzBoost::SetBeta;
TLorentzBoost &(TLorentzBoost::*TLorentzBoost_SetBeta2)(const TUnitVector &bhat, LDouble_t beta) =
   &TLorentzBoost::SetBeta;
TLorentzBoost &(TLorentzBoost::*TLorentzBoost_SetBeta3)(const TThreeVectorReal &beta) =
   &TLorentzBoost::SetBeta;
TLorentzBoost &(TLorentzBoost::*TLorentzBoost_SetBeta4)(const TFourVectorReal &p) =
   &TLorentzBoost::SetBeta;

void TLorentzBoost_Print(TLorentzBoost &obj) {
   obj.Print();
}


void TLorentzTransform_Print(TLorentzTransform &obj) {
   obj.Print();
}


Complex_t *TPauliMatrix_getitem(const TPauliMatrix &obj, Int_t index) {
   return obj[index];
}

void (TPauliMatrix::*TPauliMatrix_Decompose1)(LDouble_t &a, TThreeVectorReal &b) const =
   &TPauliMatrix::Decompose;
void (TPauliMatrix::*TPauliMatrix_Decompose2)(Complex_t &a, TThreeVectorComplex &b) const =
   &TPauliMatrix::Decompose;

TPauliMatrix &(TPauliMatrix::*TPauliMatrix_Compose1)(LDouble_t a, const TThreeVectorReal &polar) =
   &TPauliMatrix::Compose;
TPauliMatrix &(TPauliMatrix::*TPauliMatrix_Compose2)(const Complex_t &a, const TThreeVectorComplex &polar) =
   &TPauliMatrix::Compose;

TPauliMatrix &(TPauliMatrix::*TPauliMatrix_SetDiagonal1)(const Complex_t &a) =
   &TPauliMatrix::SetDiagonal;
TPauliMatrix &(TPauliMatrix::*TPauliMatrix_SetDiagonal2)(const Complex_t &a11, const Complex_t &a22) =
   &TPauliMatrix::SetDiagonal;

TPauliMatrix &(TPauliMatrix::*TPauliMatrix_SetRotation)(const TThreeRotation &rotOp) =
   &TPauliMatrix::SetRotation;
TPauliMatrix &(TPauliMatrix::*TPauliMatrix_SetRotation1)(const TThreeVectorReal &axis) =
   &TPauliMatrix::SetRotation;
TPauliMatrix &(TPauliMatrix::*TPauliMatrix_SetRotation2)(const TUnitVector &axis, LDouble_t angle) =
   &TPauliMatrix::SetRotation;
TPauliMatrix &(TPauliMatrix::*TPauliMatrix_SetRotation3)(LDouble_t phi, LDouble_t theta, LDouble_t psi) =
   &TPauliMatrix::SetRotation;

void TPauliMatrix_Print(TPauliMatrix &obj) {
   obj.Print();
}

Complex_t *TDiracMatrix_getitem(const TDiracMatrix &obj, int index) {
   return obj[index];
}

Complex_t (TDiracMatrix::*TDiracMatrix_Component1)(const EDiracIndex i) const =
   &TDiracMatrix::Component;
Complex_t (TDiracMatrix::*TDiracMatrix_Component2)(const EDiracIndex i, const EDiracIndex j) const =
   &TDiracMatrix::Component;

TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetDiagonal1)(const Complex_t &a) =
   &TDiracMatrix::SetDiagonal;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetDiagonal4)(const Complex_t &a11, const Complex_t &a22, const Complex_t &a33, const Complex_t &a44) =
   &TDiracMatrix::SetDiagonal;

TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetUUbar1)(const TFourVectorReal &p) =
   &TDiracMatrix::SetUUbar;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetUUbar2)(const TFourVectorReal &p, const Float_t helicity) =
   &TDiracMatrix::SetUUbar;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetUUbar3)(const TFourVectorReal &p, const TDiracSpinor &u) =
   &TDiracMatrix::SetUUbar;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetUUbar4)(const TFourVectorReal &p, const TThreeVectorReal &polar) =
   &TDiracMatrix::SetUUbar;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetUUbar5)(const TFourVectorReal &p, const TPauliMatrix &density) =
   &TDiracMatrix::SetUUbar;

TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetVVbar1)(const TFourVectorReal &p) =
   &TDiracMatrix::SetVVbar;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetVVbar2)(const TFourVectorReal &p, const Float_t helicity) =
   &TDiracMatrix::SetVVbar;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetVVbar3)(const TFourVectorReal &p, const TDiracSpinor &u) =
   &TDiracMatrix::SetVVbar;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetVVbar4)(const TFourVectorReal &p, const TThreeVectorReal &polar) =
   &TDiracMatrix::SetVVbar;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetVVbar5)(const TFourVectorReal &p, const TPauliMatrix &density) =
   &TDiracMatrix::SetVVbar;

TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetRotation)(const TThreeRotation &rotOp) =
   &TDiracMatrix::SetRotation;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetRotation1)(const TThreeVectorReal &axis) =
   &TDiracMatrix::SetRotation;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetRotation2)(LDouble_t phi, LDouble_t theta, LDouble_t psi) =
   &TDiracMatrix::SetRotation;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetRotation3)(LDouble_t betaX, LDouble_t betaY, LDouble_t betaZ) =
   &TDiracMatrix::SetRotation;

TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetBoost)(const TLorentzBoost &boostOp) =
   &TDiracMatrix::SetBoost;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetBoost1)(const TThreeVectorReal &beta) =
   &TDiracMatrix::SetBoost;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetBoost2)(const TUnitVector &bhat, const LDouble_t &beta) =
   &TDiracMatrix::SetBoost;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_SetBoost3)(LDouble_t betaX, LDouble_t betaY, LDouble_t betaZ) =
   &TDiracMatrix::SetBoost;

TDiracMatrix &(TDiracMatrix::*TDiracMatrix_Slash1)(const TFourVectorReal &p) =
   &TDiracMatrix::Slash;
TDiracMatrix &(TDiracMatrix::*TDiracMatrix_Slash2)(const TFourVectorComplex &A) =
   &TDiracMatrix::Slash;

void TDiracMatrix_Print(TDiracMatrix &obj) {
   obj.Print();
}


Complex_t &TPauliSpinor_getitem(const TPauliSpinor &obj, int index) {
   return obj[index];
}

LDouble_t (TPauliSpinor::*TPauliSpinor_DistanceTo1)(const TPauliSpinor &another) const =
   &TPauliSpinor::DistanceTo;
LDouble_t (TPauliSpinor::*TPauliSpinor_DistanceTo2)(const Complex_t *array) const =
   &TPauliSpinor::DistanceTo;

Bool_t (TPauliSpinor::*TPauliSpinor_equals)(const TPauliSpinor &other) const =
   &TPauliSpinor::operator==;
Bool_t (TPauliSpinor::*TPauliSpinor_nequals)(const TPauliSpinor &other) const =
   &TPauliSpinor::operator!=;

TPauliSpinor &(TPauliSpinor::*TPauliSpinor_Normalize)() =
   &TPauliSpinor::Normalize;
TPauliSpinor &(TPauliSpinor::*TPauliSpinor_Normalize1)(const LDouble_t &norm) =
   &TPauliSpinor::Normalize;

TPauliSpinor &(TPauliSpinor::*TPauliSpinor_SetPolar)(const TUnitVector &pol) =
   &TPauliSpinor::SetPolar;
TPauliSpinor &(TPauliSpinor::*TPauliSpinor_SetPolar2)(const LDouble_t &theta, const LDouble_t &phi) =
   &TPauliSpinor::SetPolar;

TPauliSpinor &(TPauliSpinor::*TPauliSpinor_Rotate)(const TThreeRotation &rotOp) =
   &TPauliSpinor::Rotate;
TPauliSpinor &(TPauliSpinor::*TPauliSpinor_Rotate1)(const TThreeVectorReal &axis) =
   &TPauliSpinor::Rotate;
TPauliSpinor &(TPauliSpinor::*TPauliSpinor_Rotate2)(const TUnitVector &axis, LDouble_t angle) =
   &TPauliSpinor::Rotate;
TPauliSpinor &(TPauliSpinor::*TPauliSpinor_Rotate3)(const LDouble_t &phi, const LDouble_t &theta, const LDouble_t &psi) =
   &TPauliSpinor::Rotate;

void TPauliSpinor_Print(TPauliSpinor &obj) {
   obj.Print();
}


Complex_t &TDiracSpinor_getitem(const TDiracSpinor &obj, int index) {
   return obj[index];
}

LDouble_t (TDiracSpinor::*TDiracSpinor_DistanceTo1)(const TDiracSpinor &another) const =
   &TDiracSpinor::DistanceTo;
LDouble_t (TDiracSpinor::*TDiracSpinor_DistanceTo2)(const Complex_t *array) const =
   &TDiracSpinor::DistanceTo;

Bool_t (TDiracSpinor::*TDiracSpinor_equals)(const TDiracSpinor &other) const =
   &TDiracSpinor::operator==;
Bool_t (TDiracSpinor::*TDiracSpinor_nequals)(const TDiracSpinor &other) const =
   &TDiracSpinor::operator!=;

TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Normalize)(const LDouble_t &norm) =
   &TDiracSpinor::Normalize;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Normalize1)(const TFourVectorReal &p) =
   &TDiracSpinor::Normalize;

TDiracSpinor &(TDiracSpinor::*TDiracSpinor_SetStateU)(const TFourVectorReal &p, Float_t helicity) =
   &TDiracSpinor::SetStateU;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_SetStateU2)(const TFourVectorReal &p, const TUnitVector &polar) =
   &TDiracSpinor::SetStateU;

TDiracSpinor &(TDiracSpinor::*TDiracSpinor_SetStateV)(const TFourVectorReal &p, Float_t helicity) =
   &TDiracSpinor::SetStateV;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_SetStateV2)(const TFourVectorReal &p, const TUnitVector &polar) =
   &TDiracSpinor::SetStateV;

TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Rotate)(const TThreeRotation &rotOp) =
   &TDiracSpinor::Rotate;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Rotate1)(const TThreeVectorReal &axis) =
   &TDiracSpinor::Rotate;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Rotate2)(const TUnitVector &axis, LDouble_t angle) =
   &TDiracSpinor::Rotate;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Rotate3)(const LDouble_t &phi, const LDouble_t &theta, const LDouble_t &psi) =
   &TDiracSpinor::Rotate;

TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Boost)(const TLorentzBoost &boostOp) =
   &TDiracSpinor::Boost;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Boost1)(const TThreeVectorReal &beta) =
   &TDiracSpinor::Boost;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Boost2)(const TUnitVector &bhat, LDouble_t beta) =
   &TDiracSpinor::Boost;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Boost3)(LDouble_t betaX, LDouble_t betaY, LDouble_t betaZ) =
   &TDiracSpinor::Boost;
TDiracSpinor &(TDiracSpinor::*TDiracSpinor_Boost4)(const LDouble_t *beta) =
   &TDiracSpinor::Boost;

void TDiracSpinor_Print(TDiracSpinor &obj) {
   obj.Print();
}


TLepton &(TLepton::*TLepton_SetMom)(const TFourVectorReal &p) =
   &TLepton::SetMom;
TLepton &(TLepton::*TLepton_SetMom1)(const TThreeVectorReal &p) =
   &TLepton::SetMom;

void TLepton_Print(TLepton &obj) {
   obj.Print();
}


TPhoton &(TPhoton::*TPhoton_SetMom)(const TFourVectorReal &p) =
   &TPhoton::SetMom;
TPhoton &(TPhoton::*TPhoton_SetMom1)(const TThreeVectorReal &p) =
   &TPhoton::SetMom;

void TPhoton_Print(TPhoton &obj) {
   obj.Print();
}


void TGhoston_Print(TGhoston &obj) {
   obj.Print();
}


void TCrossSection_Print(TCrossSection &obj) {
   obj.Print();
}

///////////////////////////////////////////////////////////
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

   boost::python::class_<TFourVectorComplex, TFourVectorComplex*,
          boost::python::bases<TThreeVectorReal> >
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

   boost::python::class_<TLorentzTransform, TLorentzTransform*>
         ("TLorentzTransform",
          "general Lorentz transform in 4-space")
      .def(boost::python::init<const TLorentzTransform &>())
      .def("SetResolution", &TLorentzTransform::SetResolution)
      .def("Resolution", &TLorentzTransform::Resolution)
      .def("IsNull", &TLorentzTransform::IsNull)
      .def("IsRotation", &TLorentzTransform::IsRotation)
      .def("IsLorentzBoost", &TLorentzTransform::IsLorentzBoost)
      .def("IsOrthogonal", &TLorentzTransform::IsOrthogonal)
      .def("IsIsochronous", &TLorentzTransform::IsIsochronous)
      .def("IsProper", &TLorentzTransform::IsProper)
      .def("Factorize", &TLorentzTransform::Factorize)
      .def("Null", &TLorentzBoost::Null,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("TimeRev", &TLorentzBoost::TimeRev,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SpaceInv", &TLorentzBoost::SpaceInv,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Transpose", &TLorentzBoost::Transpose,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Invert", &TLorentzBoost::Invert,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def(boost::python::self_ns::self *= TLorentzTransform())
      .def(boost::python::self_ns::self * TLorentzTransform())
      .def(boost::python::self_ns::self * TFourVectorComplex())
      .def(boost::python::self_ns::self * TFourVectorReal())
      .def("__eq__", &TLorentzTransform::operator==)
      .def("__ne__", &TLorentzTransform::operator!=)
      .def("Print", &TLorentzBoost::Print)
      .def("Print", &TLorentzBoost_Print)
   ;

   boost::python::class_<TThreeRotation, TThreeRotation*,
          boost::python::bases<TLorentzTransform> >
         ("TThreeRotation",
          "rotation operation in 3-space")
      .def(boost::python::init<const TThreeVectorReal &>())
      .def(boost::python::init<const TUnitVector &, const LDouble_t>())
      .def(boost::python::init<LDouble_t, LDouble_t, LDouble_t>())
      .def(boost::python::init<const TThreeRotation &>())
      .def("Axis", &TThreeRotation::Axis)
      .def("GetAxis", &TThreeRotation::GetAxis)
      .def("GetEuler", &TThreeRotation::GetEuler)
      .def(boost::python::self_ns::self *= TThreeRotation())
      .def(boost::python::self_ns::self * TThreeVectorReal())
      .def(boost::python::self_ns::self * TThreeVectorComplex())
      .def(boost::python::self_ns::self * TThreeRotation())
      .def("GetEuler", &TThreeRotation::GetEuler)
      .def("Null", &TThreeRotation::Null,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Transpose", &TThreeRotation::Transpose,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Invert", &TThreeRotation::Invert,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetAxis", TThreeRotation_SetAxis,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetAxis", TThreeRotation_SetAxis2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetEuler", &TThreeRotation::SetEuler,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Print", &TThreeRotation::Print)
      .def("Print", &TThreeRotation_Print)
   ;

   boost::python::class_<TLorentzBoost, TLorentzBoost*,
          boost::python::bases<TLorentzTransform> >
         ("TLorentzBoost",
          "Lorentz boost in 4-space")
      .def(boost::python::init<const TThreeVectorReal &>())
      .def(boost::python::init<const TFourVectorReal &>())
      .def(boost::python::init<LDouble_t, LDouble_t, LDouble_t>())
      .def(boost::python::init<const LDouble_t*>())
      .def(boost::python::init<const TUnitVector &, LDouble_t>())
      .def(boost::python::init<const TLorentzBoost &>())
      .def("Beta", &TLorentzBoost::Beta)
      .def("Null", &TLorentzBoost::Null,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Transpose", &TLorentzBoost::Transpose,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Invert", &TLorentzBoost::Invert,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBeta", TLorentzBoost_SetBeta,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBeta", TLorentzBoost_SetBeta1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBeta", TLorentzBoost_SetBeta2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBeta", TLorentzBoost_SetBeta3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBeta", TLorentzBoost_SetBeta4,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Print", &TLorentzBoost::Print)
      .def("Print", &TLorentzBoost_Print)
   ;

   boost::python::class_<TPauliMatrix, TPauliMatrix*>
         ("TPauliMatrix",
          "general 3x3 complex matrix")
      .def(boost::python::init<const EPauliIndex>())
      .def(boost::python::init<Float_t>())
      .def(boost::python::init<Double_t>())
      .def(boost::python::init<LDouble_t>())
      .def(boost::python::init<Complex_t>())
      .def(boost::python::init<Complex_t, const TThreeVectorComplex &>())
      .def(boost::python::init<const TPauliMatrix &>())
      .def("__getitem__", &TPauliMatrix_getitem,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetResolution", &TPauliMatrix::SetResolution)
      .def("Resolution", &TPauliMatrix::Resolution)
      .def("IsIdentity", &TPauliMatrix::IsIdentity)
      .def("IsUnitary", &TPauliMatrix::IsUnitary)
      .def("IsDiagonal", &TPauliMatrix::IsDiagonal)
      .def("IsHermetian", &TPauliMatrix::IsHermetian)
      .def("IsIdempotent", &TPauliMatrix::IsIdempotent)
      .def("Trace", &TPauliMatrix::Trace)
      .def("Determ", &TPauliMatrix::Determ)
      .def("Decompose", TPauliMatrix_Decompose1)
      .def("Decompose", TPauliMatrix_Decompose2)
      .def("GetDiagonal", &TPauliMatrix::GetDiagonal)
      .def(boost::python::self_ns::self += TPauliMatrix())
      .def(boost::python::self_ns::self += Complex_t())
      .def(boost::python::self_ns::self += LDouble_t())
      .def(boost::python::self_ns::self -= TPauliMatrix())
      .def(boost::python::self_ns::self -= Complex_t())
      .def(boost::python::self_ns::self -= LDouble_t())
      .def(boost::python::self_ns::self *= TPauliMatrix())
      .def(boost::python::self_ns::self *= Complex_t())
      .def(boost::python::self_ns::self *= LDouble_t())
      .def(boost::python::self_ns::self /= TPauliMatrix())
      .def(boost::python::self_ns::self /= Complex_t())
      .def(boost::python::self_ns::self /= LDouble_t())
      .def(boost::python::self_ns::self + TPauliMatrix())
      .def(boost::python::self_ns::self + Complex_t())
      .def(Complex_t() + boost::python::self_ns::self)
      .def(boost::python::self_ns::self - TPauliMatrix())
      .def(boost::python::self_ns::self - Complex_t())
      .def(Complex_t() - boost::python::self_ns::self)
      .def(boost::python::self_ns::self * TPauliMatrix())
      .def(boost::python::self_ns::self * Complex_t())
      .def(Complex_t() * boost::python::self_ns::self)
      .def(boost::python::self_ns::self / TPauliMatrix())
      .def(boost::python::self_ns::self / Complex_t())
      .def(Complex_t() / boost::python::self_ns::self)
      .def("__eq__", &TPauliMatrix::operator==)
      .def("__ne__", &TPauliMatrix::operator!=)
      .def("Zero", &TPauliMatrix::Zero,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Conj", &TPauliMatrix::Conj,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Invert", &TPauliMatrix::Invert,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Adjoint", &TPauliMatrix::Adjoint,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Transpose", &TPauliMatrix::Transpose,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Compose", TPauliMatrix_Compose1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Compose", TPauliMatrix_Compose2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetDiagonal", TPauliMatrix_SetDiagonal1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetDiagonal", TPauliMatrix_SetDiagonal2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetDensity", &TPauliMatrix::SetDensity,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetRotation", TPauliMatrix_SetRotation,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetRotation", TPauliMatrix_SetRotation1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetRotation", TPauliMatrix_SetRotation2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetRotation", TPauliMatrix_SetRotation3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SimTransform", &TPauliMatrix::SimTransform,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("UniTransform", &TPauliMatrix::UniTransform,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("__neg__", &TPauliMatrix::operator-)
      .def("Print", &TPauliMatrix::Print)
      .def("Print", &TPauliMatrix_Print)
   ;

   boost::python::class_<TDiracMatrix, TDiracMatrix*>
         ("TDiracMatrix",
          "general 4x4 complex matrix")
      .def(boost::python::init<const EDiracIndex>())
      .def(boost::python::init<const EDiracIndex, const EDiracIndex>())
      .def(boost::python::init<LDouble_t>())
      .def(boost::python::init<Complex_t>())
      .def(boost::python::init<const TLorentzTransform &>())
      .def(boost::python::init<const TDiracMatrix &>())
      .def("__getitem__", &TDiracMatrix_getitem,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetResolution", &TDiracMatrix::SetResolution)
      .def("Resolution", &TDiracMatrix::Resolution)
      .def("IsIdentity", &TDiracMatrix::IsIdentity)
      .def("IsUnitary", &TDiracMatrix::IsUnitary)
      .def("IsDiagonal", &TDiracMatrix::IsDiagonal)
      .def("IsAntiUnitary", &TDiracMatrix::IsDiagonal)
      .def("IsHermetian", &TDiracMatrix::IsHermetian)
      .def("IsIdempotent", &TDiracMatrix::IsIdempotent)
      .def("Trace", &TDiracMatrix::Trace)
      .def("Determ", &TDiracMatrix::Determ)
      .def("Component", TDiracMatrix_Component1)
      .def("Component", TDiracMatrix_Component2)
      .def("GetDiagonal", &TDiracMatrix::GetDiagonal)
      .def(boost::python::self_ns::self += TDiracMatrix())
      .def(boost::python::self_ns::self += Complex_t())
      .def(boost::python::self_ns::self += LDouble_t())
      .def(boost::python::self_ns::self -= TDiracMatrix())
      .def(boost::python::self_ns::self -= Complex_t())
      .def(boost::python::self_ns::self -= LDouble_t())
      .def(boost::python::self_ns::self *= TDiracMatrix())
      .def(boost::python::self_ns::self *= Complex_t())
      .def(boost::python::self_ns::self *= LDouble_t())
      .def(boost::python::self_ns::self /= TDiracMatrix())
      .def(boost::python::self_ns::self /= Complex_t())
      .def(boost::python::self_ns::self /= LDouble_t())
      .def(boost::python::self_ns::self + TDiracMatrix())
      .def(boost::python::self_ns::self + Complex_t())
      .def(boost::python::self_ns::self + LDouble_t())
      .def(Complex_t() + boost::python::self_ns::self)
      .def(LDouble_t() + boost::python::self_ns::self)
      .def(boost::python::self_ns::self - TDiracMatrix())
      .def(boost::python::self_ns::self - Complex_t())
      .def(boost::python::self_ns::self - LDouble_t())
      .def(Complex_t() - boost::python::self_ns::self)
      .def(LDouble_t() - boost::python::self_ns::self)
      .def(boost::python::self_ns::self * TDiracMatrix())
      .def(boost::python::self_ns::self * Complex_t())
      .def(boost::python::self_ns::self * LDouble_t())
      .def(Complex_t() * boost::python::self_ns::self)
      .def(LDouble_t() * boost::python::self_ns::self)
      .def(boost::python::self_ns::self / TDiracMatrix())
      .def(boost::python::self_ns::self / Complex_t())
      .def(boost::python::self_ns::self / LDouble_t())
      .def(Complex_t() / boost::python::self_ns::self)
      .def(LDouble_t() / boost::python::self_ns::self)
      .def("__eq__", &TDiracMatrix::operator==)
      .def("__ne__", &TDiracMatrix::operator!=)
      .def("Zero", &TDiracMatrix::Zero,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Conj", &TDiracMatrix::Conj,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Invert", &TDiracMatrix::Invert,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Adjoint", &TDiracMatrix::Adjoint,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Transpose", &TDiracMatrix::Transpose,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetDiagonal", TDiracMatrix_SetDiagonal1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetDiagonal", TDiracMatrix_SetDiagonal4,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetUUbar", TDiracMatrix_SetUUbar1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetUUbar", TDiracMatrix_SetUUbar2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetUUbar", TDiracMatrix_SetUUbar3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetUUbar", TDiracMatrix_SetUUbar4,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetUUbar", TDiracMatrix_SetUUbar5,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetVVbar", TDiracMatrix_SetVVbar1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetVVbar", TDiracMatrix_SetVVbar2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetVVbar", TDiracMatrix_SetVVbar3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetVVbar", TDiracMatrix_SetVVbar4,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetVVbar", TDiracMatrix_SetVVbar5,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetRotation", TDiracMatrix_SetRotation,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetRotation", TDiracMatrix_SetRotation1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetRotation", TDiracMatrix_SetRotation2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetRotation", TDiracMatrix_SetRotation3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBoost", TDiracMatrix_SetBoost,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBoost", TDiracMatrix_SetBoost1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBoost", TDiracMatrix_SetBoost2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetBoost", TDiracMatrix_SetBoost3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetTransform", &TDiracMatrix::SetTransform,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SimTransform", &TDiracMatrix::SimTransform,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("UniTransform", &TDiracMatrix::UniTransform,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Slash", TDiracMatrix_Slash1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Slash", TDiracMatrix_Slash2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("__neg__", &TDiracMatrix::operator-)
      .def("Print", &TDiracMatrix::Print)
      .def("Print", &TDiracMatrix_Print)
   ;

   boost::python::class_<TPauliSpinor, TPauliSpinor*>
         ("TPauliSpinor",
          "general 2x1 complex column matrix")
      .def(boost::python::init<const Complex_t &, const Complex_t &>())
      .def(boost::python::init<const Complex_t *>())
      .def(boost::python::init<LDouble_t, LDouble_t>())
      .def(boost::python::init<const TUnitVector &>())
      .def(boost::python::init<const TPauliSpinor &>())
      .def("__getitem__", &TPauliSpinor_getitem,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetResolution", &TPauliSpinor::SetResolution)
      .def("Resolution", &TPauliSpinor::Resolution)
      .def("Norm", &TPauliSpinor::Norm)
      .def("NormSqr", &TPauliSpinor::NormSqr)
      .def("GetPolar", &TPauliSpinor::GetPolar)
      .def("Polar", &TPauliSpinor::Polar)
      .def("DistanceTo", TPauliSpinor_DistanceTo1)
      .def("DistanceTo", TPauliSpinor_DistanceTo2)
      .def(boost::python::self_ns::self += TPauliSpinor())
      .def(boost::python::self_ns::self -= TPauliSpinor())
      .def(boost::python::self_ns::self *= Complex_t())
      .def(boost::python::self_ns::self /= Complex_t())
      .def(boost::python::self_ns::self + TPauliSpinor())
      .def(boost::python::self_ns::self - TPauliSpinor())
      .def(boost::python::self_ns::self * Complex_t())
      .def(boost::python::self_ns::self * LDouble_t())
      .def(Complex_t() * boost::python::self_ns::self)
      .def(LDouble_t() * boost::python::self_ns::self)
      .def(boost::python::self_ns::self / Complex_t())
      .def("__eq__", TPauliSpinor_equals)
      .def("__ne__", TPauliSpinor_nequals)
      .def("Zero", &TPauliSpinor::Zero,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Conj", &TPauliSpinor::Conj,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Normalize", TPauliSpinor_Normalize,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Normalize", TPauliSpinor_Normalize1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetPolar", TPauliSpinor_SetPolar,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetPolar", TPauliSpinor_SetPolar2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Operate", &TPauliSpinor::Operate,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TPauliSpinor_Rotate,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TPauliSpinor_Rotate1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TPauliSpinor_Rotate2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TPauliSpinor_Rotate3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("InnerProd", &TPauliSpinor::InnerProd)
      .def("ScalarProd", &TPauliSpinor::ScalarProd)
      .def("__neg__", &TPauliSpinor::operator-)
      .def("Print", &TPauliSpinor::Print)
      .def("Print", &TPauliSpinor_Print)
   ;

   boost::python::class_<TDiracSpinor, TDiracSpinor*>
         ("TDiracSpinor",
          "general 4x1 complex column matrix")
      .def(boost::python::init<const Complex_t &, const Complex_t &, const Complex_t &, const Complex_t &>())
      .def(boost::python::init<const Complex_t *>())
      .def(boost::python::init<const TFourVectorReal &, const Float_t>())
      .def(boost::python::init<const TFourVectorReal &, const TUnitVector &>())
      .def(boost::python::init<const TDiracSpinor &>())
      .def("__getitem__", &TDiracSpinor_getitem,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetResolution", &TDiracSpinor::SetResolution)
      .def("Resolution", &TDiracSpinor::Resolution)
      .def("Norm", &TDiracSpinor::Norm)
      .def("NormSqr", &TDiracSpinor::NormSqr)
      .def("Upper", &TDiracSpinor::Upper)
      .def("Lower", &TDiracSpinor::Lower)
      .def("DistanceTo", TDiracSpinor_DistanceTo1)
      .def("DistanceTo", TDiracSpinor_DistanceTo2)
      .def(boost::python::self_ns::self += TDiracSpinor())
      .def(boost::python::self_ns::self -= TDiracSpinor())
      .def(boost::python::self_ns::self *= Complex_t())
      .def(boost::python::self_ns::self /= Complex_t())
      .def(boost::python::self_ns::self + TDiracSpinor())
      .def(boost::python::self_ns::self - TDiracSpinor())
      .def(boost::python::self_ns::self * Complex_t())
      .def(boost::python::self_ns::self * LDouble_t())
      .def(Complex_t() * boost::python::self_ns::self)
      .def(LDouble_t() * boost::python::self_ns::self)
      .def(boost::python::self_ns::self / Complex_t())
      .def("__eq__", TDiracSpinor_equals)
      .def("__ne__", TDiracSpinor_nequals)
      .def("Zero", &TDiracSpinor::Zero,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Conj", &TDiracSpinor::Conj,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Bar", &TDiracSpinor::Bar,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetUpper", &TDiracSpinor::SetUpper)
      .def("SetLower", &TDiracSpinor::SetLower)
      .def("Normalize", TDiracSpinor_Normalize,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Normalize", TDiracSpinor_Normalize1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetStateU", TDiracSpinor_SetStateU,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetStateU", TDiracSpinor_SetStateU2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetStateV", TDiracSpinor_SetStateV,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetStateV", TDiracSpinor_SetStateV2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Operate", &TDiracSpinor::Operate,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TDiracSpinor_Rotate,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TDiracSpinor_Rotate1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TDiracSpinor_Rotate2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Rotate", TDiracSpinor_Rotate3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TDiracSpinor_Boost,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TDiracSpinor_Boost1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TDiracSpinor_Boost2,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TDiracSpinor_Boost3,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Boost", TDiracSpinor_Boost4,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("BoostToRest", &TDiracSpinor::BoostToRest,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("BoostFromRest", &TDiracSpinor::BoostFromRest,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("InnerProd", &TDiracSpinor::InnerProd)
      .def("ScalarProd", &TDiracSpinor::ScalarProd)
      .def("__neg__", &TDiracSpinor::operator-)
      .def("Print", &TDiracSpinor::Print)
      .def("Print", &TDiracSpinor_Print)
   ;

   boost::python::class_<TLepton, TLepton*>
         ("TLepton",
          "object representing a lepton")
      .def(boost::python::init<LDouble_t>())
      .def(boost::python::init<const TFourVectorReal &>())
      .def(boost::python::init<const TFourVectorReal &, LDouble_t>())
      .def(boost::python::init<const TThreeVectorReal &>())
      .def(boost::python::init<const TThreeVectorReal &, LDouble_t>())
      .def(boost::python::init<const TLepton>())
      .def("Mass", &TLepton::Mass)
      .def("Mom", &TLepton::Mom)
      .def("Pol", &TLepton::Pol)
      .def("SDM", &TLepton::SDM,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetMom", TLepton_SetMom,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetMom", TLepton_SetMom1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetPol", &TLepton::SetPol,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("AllPol", &TLepton::AllPol,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("Print", &TLepton::Print)
      .def("Print", &TLepton_Print)
   ;

   boost::python::class_<TPhoton, TPhoton*>
         ("TPhoton",
          "object representing a photon")
      .def(boost::python::init<const TFourVectorReal &>())
      .def(boost::python::init<const TThreeVectorReal &>())
      .def(boost::python::init<const TPhoton &>())
      .def("Mom", &TPhoton::Mom)
      .def("Pol", &TPhoton::Pol)
      .def("Eps", &TPhoton::Eps)
      .def("EpsStar", &TPhoton::EpsStar)
      .def("SDM", &TPhoton::SDM,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetMom", TPhoton_SetMom,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetMom", TPhoton_SetMom1,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetPol", &TPhoton::SetPol,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("AllPol", &TPhoton::AllPol,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetPlanePolarization", &TPhoton::SetPlanePolarization,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("SetEllipticalPolarization", &TPhoton::SetEllipticalPolarization,
           boost::python::return_value_policy<boost::python::reference_existing_object>())
      .def("GetPolarizationPlane", &TPhoton::GetPolarizationPlane)
      .def("Print", &TPhoton::Print)
      .def("Print", &TPhoton_Print)
   ;

   boost::python::class_<TGhoston, TGhoston*,
          boost::python::bases<TPhoton> >
         ("TGhoston",
          "object representing a ghost photon (longitudinal mode)")
      .def(boost::python::init<const TFourVectorReal &>())
      .def(boost::python::init<const TThreeVectorReal &>())
      .def("Eps", &TGhoston::Eps)
      .def("EpsStar", &TGhoston::EpsStar)
      .def("Print", &TGhoston::Print)
      .def("Print", &TGhoston_Print)
   ;

   boost::python::class_<TCrossSection, TCrossSection*>
         ("TCrossSection",
          "methods for computing various QED polarized cross sections")
      .def("Compton", &TCrossSection::Compton)
      .staticmethod("Compton")
      .def("Bremsstrahlung", &TCrossSection::Bremsstrahlung)
      .staticmethod("Bremsstrahlung")
      .def("PairProduction", &TCrossSection::PairProduction)
      .staticmethod("PairProduction")
      .def("TripletProduction", &TCrossSection::TripletProduction)
      .staticmethod("TripletProduction")
      .def("eeBremsstrahlung", &TCrossSection::eeBremsstrahlung)
      .staticmethod("eeBremsstrahlung")
      .def("Print", &TCrossSection::Print)
      .def("Print", &TCrossSection_Print)
   ;
}
