//
// TFourVectorComplex.cxx
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
// defined for both Double_t and Complex_t generic types.  For complex
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
 
#include "TFourVectorComplex.h"
#include "TLorentzBoost.h"

#include <iostream>
using namespace std;

ClassImp(TFourVectorComplex)


TFourVectorComplex &TFourVectorComplex::Transform
                   (const TLorentzTransform &xformOp)
{
   TFourVectorComplex temp = xformOp*(*this);
   return (*this = temp);
}

TFourVectorComplex &TFourVectorComplex::Boost
                   (const TLorentzBoost &boostOp)
{
   TFourVectorComplex temp = boostOp*(*this);
   return (*this = temp);
}

TFourVectorComplex &TFourVectorComplex::Boost(const LDouble_t betaX,
                                              const LDouble_t betaY,
                                              const LDouble_t betaZ)
{
   TLorentzBoost boostOp(betaX,betaY,betaZ);
   return Boost(boostOp);
}

TFourVectorComplex &TFourVectorComplex::Boost(const LDouble_t *beta)
{
   TLorentzBoost boostOp(beta);
   return Boost(boostOp);
}

TFourVectorComplex &TFourVectorComplex::Boost(const TThreeVectorReal &beta)
{
   TLorentzBoost boostOp(beta);
   return Boost(boostOp);
}

TFourVectorComplex &TFourVectorComplex::Boost
                   (const TUnitVector &bhat, const LDouble_t beta)
{
   TLorentzBoost boostOp(bhat,beta);
   return Boost(boostOp);
}

TFourVectorComplex &TFourVectorComplex::BoostToRest(const TFourVector &p)
{
   TLorentzBoost boostOp(p/p[0]);
   return Boost(boostOp);
}

TFourVectorComplex &TFourVectorComplex::BoostFromRest(const TFourVector &p)
{
   TLorentzBoost boostOp(-p/p[0]);
   return Boost(boostOp);
}

void TFourVectorComplex::Streamer(TBuffer &buf)
{
   // Put/get a complex four-vector to/from stream buffer buf.
   // This method assumes that complex is stored in memory as LDouble_t[2].

   Double_t vector[8];
   if (buf.IsReading()) {
      buf.ReadStaticArray(vector);
      fVector[0] = Complex_t(vector[0], vector[1]);
      fVector[1] = Complex_t(vector[2], vector[3]);
      fVector[2] = Complex_t(vector[4], vector[5]);
      fVector[3] = Complex_t(vector[6], vector[7]);
   } else {
      vector[0] = fVector[0].real();
      vector[1] = fVector[0].imag();
      vector[2] = fVector[1].real();
      vector[3] = fVector[1].imag();
      vector[4] = fVector[2].real();
      vector[5] = fVector[2].imag();
      vector[6] = fVector[3].real();
      vector[7] = fVector[3].imag();
      buf.WriteArray(vector, 8);
   }
}

void TFourVectorComplex::Print(Option_t *option)
{
   // Output an ascii representation for complex four-vector.

   cout << "TFourVectorComplex(" << fVector[0] << "," << fVector[1] << ","
        << fVector[2] << "," << fVector[3] << ")" << endl;
}
 
//______________________________________________________________________________


#ifdef R__HPUX

//______________________________________________________________________________
//  These functions should be inline
//______________________________________________________________________________

#endif
