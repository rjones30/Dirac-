//
// TLorentzBoost.cxx
//
// author:  Richard T. Jones  11/16/98
// version:  Dec. 12, 1998  v1.00
//
/*************************************************************************
 * Copyright(c) 1998, University of Connecticut, All rights reserved.    *
 * Author: Richard T. Jones, Asst. Prof. of Physics                      *
 *                                                                       *
 * Permission to use, copy, modify and distribute this software and its  *
 * documentation for non-commercial purposes is hereby granted without   *
 * fee, provided that the above copyright notice appears in all copies   *
 * and that both the copyright notice and this permission notice appear  *
 * in the supporting documentation. The author makes no claims about the *
 * suitability of this software for any purpose.                         *
 * It is provided "as is" without express or implied warranty.           *
 *************************************************************************/
//////////////////////////////////////////////////////////////////////////
//
// Lorentz Algebra Package
//
// The present package implements all the basic algorithms dealing
// with three-vectors and four-vectors, together with their transform-
// ations.  Four-vectors are derived from three-vectors and inherit
// all of their members.  Direct access to the components is provided
// through the subscript operator [i] which covers the range 1...3 for
// three-vectors and 0...3 for four-vectors.  Transformations are
// implemented as a friend class so that they can operate directly on
// the data members of the vector, which are otherwise hidden.  The
// special transformations (rotations and boosts) inherit from the
// general class LorentzTransform.  Products of rotations are other
// rotations, whereas the product of a boost with anything is simply
// a LorentzTransform.  The LorentzTransform objects can be tested
// for the property of being a pure rotation or boost.  They can also
// implement non-isochronous and improper transformations.
//
// Rotations may be specified either by Euler angles or by a rotation
// axis.  All angles are assumed to be in radians.  Vector classes are
// defined for both Real_t and Complex_t generic types.  For complex
// vectors there are several additional member functions to deal with
// operations that are specific to complex numbers.
//
// The classes comprising this package are:
//   TThreeVectorReal is a base class
//   TThreeVectorComplex is a base class
//   TFourVectorReal is a TThreeVectorReal
//   TFourVectorComplex is a TThreeVectorComplex
//   TLorentzTransform is a base class
//   TThreeRotation is a TLorentzTransform
//   TLorentzBoost is a TLorentzTransform
// The following aliases are defined for these classes:
//   TUnitVector is an alias for TThreeVectorReal
//   TThreeVector is an alias for TThreeVectorReal
//   TFourVector is an alias for TFourVectorReal
//
// This package was developed at the University of Connecticut by
// Richard T. Jones
//
//////////////////////////////////////////////////////////////////////////
 
#include "TLorentzBoost.h"
#include "TFourVectorReal.h"
#include "TFourVectorComplex.h"

#include <iostream>
using namespace std;

ClassImp(TLorentzBoost)

 
TLorentzBoost &TLorentzBoost::SetBeta(const TThreeVectorReal &beta)
{
   Double_t betaSqr = beta.LengthSqr();
   if (betaSqr >= 1) {
      Error("TLorentzBoost::SetBeta()","attempt to boost outside light cone");
      return SetBeta(beta,0);
   }
   Double_t gamma = 1/sqrt(1-betaSqr);
   TThreeVectorReal eta = beta*gamma;
   Double_t gammaPlusOne = gamma + 1;
   fMatrix[0][0] = gamma;
   fMatrix[1][1] = 1 + eta.fVector[1]*eta.fVector[1]/gammaPlusOne;
   fMatrix[2][2] = 1 + eta.fVector[2]*eta.fVector[2]/gammaPlusOne;
   fMatrix[3][3] = 1 + eta.fVector[3]*eta.fVector[3]/gammaPlusOne;
   fMatrix[0][1] = fMatrix[1][0] = -gamma*beta.fVector[1];
   fMatrix[0][2] = fMatrix[2][0] = -gamma*beta.fVector[2];
   fMatrix[0][3] = fMatrix[3][0] = -gamma*beta.fVector[3];
   fMatrix[1][2] = fMatrix[2][1] = eta.fVector[1]*eta.fVector[2]/gammaPlusOne;
   fMatrix[1][3] = fMatrix[3][1] = eta.fVector[1]*eta.fVector[3]/gammaPlusOne;
   fMatrix[2][3] = fMatrix[3][2] = eta.fVector[2]*eta.fVector[3]/gammaPlusOne;
   return *this;
}

TLorentzBoost &TLorentzBoost::SetBeta
              (const TUnitVector &bhat, const Double_t beta)
{
   TThreeVectorReal betaV(bhat);
   betaV.Normalize(beta);
   SetBeta(betaV);
   return *this;
}

TLorentzBoost &TLorentzBoost::SetBeta(const TFourVectorReal &p)
{
   TThreeVectorReal betaV = (TThreeVectorReal)p/p[0];
   SetBeta(betaV);
   return *this;
}

void TLorentzBoost::Print(Option_t *option)
{
   // Output an ascii representation for boost matrix.

   cout << "TLorentzBoost matrix" << endl;
   cout << "(" << fMatrix[0][0] << "," << fMatrix[0][1] << ","
               << fMatrix[0][2] << "," << fMatrix[0][3] << ")" << endl;
   cout << "(" << fMatrix[1][0] << "," << fMatrix[1][1] << ","
               << fMatrix[1][2] << "," << fMatrix[1][3] << ")" << endl;
   cout << "(" << fMatrix[2][0] << "," << fMatrix[2][1] << ","
               << fMatrix[2][2] << "," << fMatrix[2][3] << ")" << endl;
   cout << "(" << fMatrix[3][0] << "," << fMatrix[3][1] << ","
               << fMatrix[3][2] << "," << fMatrix[3][3] << ")" << endl;
}

//______________________________________________________________________________


#ifdef R__HPUX

//______________________________________________________________________________
//  These functions should be inline
//______________________________________________________________________________

#endif
