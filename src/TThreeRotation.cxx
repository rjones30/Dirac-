//
// TThreeRotation.cxx
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
 
#include <math.h>
#include "TFourVectorReal.h"
#include "TFourVectorComplex.h"
#include "TThreeRotation.h"

#include <iostream>
using namespace std;

ClassImp(TThreeRotation)

 
TThreeVectorReal TThreeRotation::Axis() const
{
   LDouble_t angle(0);
   TThreeVectorReal axis;
   GetAxis(axis,angle);
   return (axis *= angle);
}

void TThreeRotation::GetAxis(TUnitVector &ahat, LDouble_t &angle) const
{
//
// A rotation is specified either by a angle about an axis in SetAxis()) or
// via Euler angles in SetEuler().  GetAxis() works backwards from a rotation
// matrix to find the rotation axis and angle that would generate it, if
// passed to SetAxis().  It is based upon the vector identity for a rotation
// of vector V by w radians about a unit vector axis Ahat.
//
//        V' = V  -  sin(w) [Ahat x V]  +  (1-cos(w)) [(Ahat.V)Ahat - V]
//
// In matrix form, the first and third terms are symmetric, whereas the second
// is antisymmetric.  This fact can be used to extract a value for sin(w) and
// the direction of Ahat.  Because of numerical errors in the Asin around 90
// degrees, rather than setting w = Asin( sin(w) ) it is better to use the 
// observation that Trace(Matrix) = 1 + 2cos(w) to extract cos(w), and use
// atan2 to extract w with good precision.
//
   LDouble_t traceM = fMatrix[1][1] + fMatrix[2][2] + fMatrix[3][3];
   ahat.fVector[1] = fMatrix[2][3] - fMatrix[3][2];
   ahat.fVector[2] = fMatrix[3][1] - fMatrix[1][3];
   ahat.fVector[3] = fMatrix[1][2] - fMatrix[2][1];
   angle = atan2(ahat.Length(),traceM-1);
   ahat.Normalize(1);
}

void TThreeRotation::GetEuler
            (LDouble_t &phi, LDouble_t &theta, LDouble_t &psi) const
{
//
// A rotation is specified either by a angle about an axis in SetAxis()) or
// via Euler angles in SetEuler().  GetEuler() works backwards from a rotation
// matrix to find the Euler angles that would generate it, if passed to
// SetEuler().  The sequence of rotations is first by phi about the z axis,
// then by theta about the y' axis, and finally by psi about the z'' axis.
//
   TThreeVectorReal v(fMatrix[3][1], fMatrix[3][2], fMatrix[3][3]);
   phi = atan2(v.fVector[2],v.fVector[1]);
   theta = atan2(v.Rho(),v.fVector[3]);
   psi = atan2(fMatrix[2][3],-fMatrix[1][3]);
   LDouble_t cosPhi = cos(phi);
   LDouble_t cosTheta = cos(theta);
   LDouble_t cosPsi = cos(psi);
   LDouble_t sinPhi = sin(phi);
   LDouble_t sinPsi = sin(psi);
   LDouble_t m11 = cosPhi*cosTheta*cosPsi - sinPhi*sinPsi;
   LDouble_t dif = abs(m11 - fMatrix[1][1]);
   if (dif > Resolution()) {
      LDouble_t sumPhiPsi = atan2(-fMatrix[2][1],fMatrix[1][1]);
      psi = sumPhiPsi - phi;
      psi += (psi > -M_PI) ? 0 : 2*M_PI;
      psi -= (psi <= M_PI) ? 0 : 2*M_PI;
   }
}

TThreeRotation &TThreeRotation::operator=(const TThreeRotation &source)
{
   *(TLorentzTransform *)this = (TLorentzTransform)source;
   return *this;
}

TThreeRotation &TThreeRotation::operator*=(const TThreeRotation &source)
{
   TThreeRotation temp(*this);
   for (Int_t i=1; i<4; i++) {
      for (Int_t j=1; j<4; j++) {
         LDouble_t sum=0;
         for (Int_t k=1; k<4; k++) {
            sum += temp.fMatrix[i][k]*source.fMatrix[k][j];
         }
         fMatrix[i][j] = sum;
      }
   }
   fMatrix[0][0] = 1;
   fMatrix[0][1] = fMatrix[1][0] = 0;
   fMatrix[0][2] = fMatrix[2][0] = 0;
   fMatrix[0][3] = fMatrix[3][0] = 0;
   return *this;
}
 
TThreeRotation &TThreeRotation::Null()
{
   TLorentzTransform::Null();
   return *this;
}
 
TThreeRotation &TThreeRotation::Transpose()
{
   TLorentzTransform::Transpose();
   return *this;
}
 
TThreeRotation &TThreeRotation::Invert()
{
   TLorentzTransform::Invert();
   return *this;
}

TThreeRotation &TThreeRotation::SetAxis(const TThreeVectorReal &axis)
{
   return SetAxis(axis,axis.Length());
}

TThreeRotation &TThreeRotation::SetAxis(const TUnitVector &ahat,
                                        const LDouble_t angle)
{
//
// This function is based upon the vector identity for a rotation of vector
// V by w radians about a unit vector axis Ahat:
//
//        V' = V  -  sin(w) [Ahat x V]  +  (1-cos(w)) [(Ahat.V)Ahat - V]
//
   TUnitVector axis(ahat);
   axis.Normalize(1);
   LDouble_t c0 = cos(angle);
   LDouble_t c1 = -sin(angle);
   LDouble_t c2 = 1-c0;
   fMatrix[0][0] = 1;
   fMatrix[0][1] = fMatrix[1][0] = 0;
   fMatrix[0][2] = fMatrix[2][0] = 0;
   fMatrix[0][3] = fMatrix[3][0] = 0;
   fMatrix[1][1] = c0 + c2*axis.fVector[1]*axis.fVector[1];
   fMatrix[2][2] = c0 + c2*axis.fVector[2]*axis.fVector[2];
   fMatrix[3][3] = c0 + c2*axis.fVector[3]*axis.fVector[3];
   fMatrix[1][2] = -c1*axis.fVector[3]
                            + c2*axis.fVector[1]*axis.fVector[2];
   fMatrix[1][3] = c1*axis.fVector[2]
                            + c2*axis.fVector[1]*axis.fVector[3];
   fMatrix[2][1] = c1*axis.fVector[3]
                            + c2*axis.fVector[1]*axis.fVector[2];
   fMatrix[2][3] = -c1*axis.fVector[1]
                            + c2*axis.fVector[2]*axis.fVector[3];
   fMatrix[3][1] = -c1*axis.fVector[2]
                            + c2*axis.fVector[1]*axis.fVector[3];
   fMatrix[3][2] = c1*axis.fVector[1]
                            + c2*axis.fVector[2]*axis.fVector[3];
   return *this;
}

TThreeRotation &TThreeRotation::SetEuler(const LDouble_t &phi,
                               const LDouble_t &theta,
                               const LDouble_t &psi)
{
//
// This function implements the y-convention for Euler angles, in common
// use in quantum mechanics.  Note that this differs from the engineering
// convention and the x-convention commonly used in classical mechanics.
// The sequence of rotations is first by phi about the z axis, then by theta
// about the y' axis, and finally by psi about the z'' axis.  For discussion
// of the relative advantages of different conventions, see Goldstein,
// Classical Mechanics, p.147 and references therein.
//
   LDouble_t cosPhi = cos(phi);
   LDouble_t cosTheta = cos(theta);
   LDouble_t cosPsi = cos(psi);
   LDouble_t sinPhi = sin(phi);
   LDouble_t sinTheta = sin(theta);
   LDouble_t sinPsi = sin(psi);
   fMatrix[0][0] = 1;
   fMatrix[0][1] = fMatrix[1][0] = 0;
   fMatrix[0][2] = fMatrix[2][0] = 0;
   fMatrix[0][3] = fMatrix[3][0] = 0;
   fMatrix[1][1] = cosPhi*cosTheta*cosPsi - sinPhi*sinPsi;
   fMatrix[1][2] = sinPhi*cosTheta*cosPsi + cosPhi*sinPsi;
   fMatrix[1][3] = -sinTheta*cosPsi;
   fMatrix[2][1] = -cosPhi*cosTheta*sinPsi - sinPhi*cosPsi;
   fMatrix[2][2] = -sinPhi*cosTheta*sinPsi + cosPhi*cosPsi;
   fMatrix[2][3] = sinTheta*sinPsi;
   fMatrix[3][1] = cosPhi*sinTheta;
   fMatrix[3][2] = sinPhi*sinTheta;
   fMatrix[3][3] = cosTheta;
   return *this;
}

TThreeVectorReal TThreeRotation::operator*(const TThreeVectorReal &vec) const
{
   TThreeVectorReal result;
   for (Int_t i=1; i<4; i++) {
      LDouble_t sum=0;
      for (Int_t j=1; j<4; j++) {
         sum += fMatrix[i][j]*vec.fVector[j];
      }
      result.fVector[i] = sum;
   }
   return result;
}
 
TThreeVectorComplex TThreeRotation::operator*
                           (const TThreeVectorComplex &vec) const
{
   TThreeVectorComplex result;
   for (Int_t i=1; i<4; i++) {
      Complex_t sum=0;
      for (Int_t j=1; j<4; j++) {
         sum += fMatrix[i][j]*vec.fVector[j];
      }
      result.fVector[i] = sum;
   }
   return result;
}

TThreeRotation TThreeRotation::operator*(const TThreeRotation &rotOp) const
{
   TThreeRotation result;
   for (Int_t i=1; i<4; i++) {
      for (Int_t j=1; j<4; j++) {
         LDouble_t sum=0;
         for (Int_t k=1; k<4; k++) {
            sum += fMatrix[i][k]*rotOp.fMatrix[k][j];
         }
         result.fMatrix[i][j] = sum;
      }
   }
   result.fMatrix[0][0] = 1;
   result.fMatrix[0][1] = result.fMatrix[1][0] = 0;
   result.fMatrix[0][2] = result.fMatrix[2][0] = 0;
   result.fMatrix[0][3] = result.fMatrix[3][0] = 0;
   return result;
}

void TThreeRotation::Print(Option_t *option)
{
   cout << "TThreeRotation matrix" << endl;
   cout << "(" << fMatrix[1][1] << ","
               << fMatrix[1][2] << "," << fMatrix[1][3] << ")" << endl;
   cout << "(" << fMatrix[2][1] << ","
               << fMatrix[2][2] << "," << fMatrix[2][3] << ")" << endl;
   cout << "(" << fMatrix[3][1] << ","
               << fMatrix[3][2] << "," << fMatrix[3][3] << ")" << endl;
}
 
//______________________________________________________________________________


#ifdef R__HPUX

//______________________________________________________________________________
//  These functions should be inline
//______________________________________________________________________________

#endif
