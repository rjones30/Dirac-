//
// TThreeVectorReal.cxx
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
//
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
 
#include <iostream>

#include "TThreeVectorReal.h"
#include "TThreeRotation.h"

ClassImp(TThreeVectorReal)


Double_t TThreeVectorReal::fResolution = 1e-12;

TThreeVectorReal &TThreeVectorReal::Rotate(const TThreeRotation &rotOp)
{
   // Operate on vector with rotation operator rotOp.

   TThreeVectorReal temp = rotOp*(*this);
   return (*this = temp);
}

TThreeVectorReal &TThreeVectorReal::Rotate
                 (const Double_t phi, const Double_t theta, const Double_t psi)
{
   // Rotate  vector using Euler angles z:phi,y':theta,z'':psi.

   TThreeRotation rotOp(phi,theta,psi);
   return Rotate(rotOp);
}

TThreeVectorReal &TThreeVectorReal::Rotate
                 (const TUnitVector &ahat, const Double_t angle)
{
   // Rotate vector by angle degrees using ahat as the axis.

   TThreeRotation rotOp(ahat,angle);
   return Rotate(rotOp);
}

void TThreeVectorReal::Streamer(TBuffer &buf)
{
   // Put/get three Real_t values to/from stream buffer buf.

   if (buf.IsReading()) {
      buf.ReadStaticArray((Real_t *)&fVector[1]);
   } else {
      buf.WriteArray((Real_t *)&fVector[1], 3);
   }
}

void TThreeVectorReal::Print(Option_t *option)
{
   // Output contents of vector in ascii.

   std::cout << "TThreeVectorReal(" << fVector[1] << "," << fVector[2]
        << "," << fVector[3] << ")" << std::endl;
}
 
//______________________________________________________________________________


#ifdef R__HPUX

//______________________________________________________________________________
//  These functions should be inline

//______________________________________________________________________________

#endif
