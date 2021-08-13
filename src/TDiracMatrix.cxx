//
// TDiracMatrix.cxx
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
// Dirac Spinor Algebra Package
//
// The present package implements all the basic algorithms dealing
// with Dirac spinors, which form a fundamental representation of the
// group SL(2,2).  The basic classes are DiracSpinor and DiracMatrix,
// which are 4x1 and 4x4 complex matrices, respectively. In this
// context any complex 4x4 matrix that operates on Dirac spinors is
// called a Dirac matrix, and not simply the four or five Dirac gamma
// matrices.  The standard representation of Dirac is used for the
// gamma matrices 0-5.  The generators of the Lorentz group are the
// Sigma (rotation generators) and Kappa (boost generators) matrices.
//
// The standard matrices are identified by a discrete index of enum
// type EDiracIndex.  A EDiracIndex can take on a value from the list
//    kDiracOne,    kDiracGamma1,   kDiracGamma2,   kDiracGamma3,
//    kDiracGamma4, kDiracGamma5,   kDiracSigma1,   kDiracSigma2,
//    kDiracSigma3, kDiracKappa1,   kDiracKappa2,   kDiracKappa3.
// The constructor invoked with two EDiracIndex values i,j returns
// i_/2 [TDiracMatrix(i),TDiracMatrix(j)] where [a,b] denotes the com-
// utator of matrices a and b, and i_ is the positive square root of
// -1.  In general Dirac matrices describe operators and Dirac spinors
// describe relativistic fermion states.  Dirac matrices are also used
// to describe mixed states, ensembles that contain mixtures of
// particles described by more than one Dirac spinor.
//
// Spinors and matrices can be transformed under rotations and boosts
// according to the commutation rules for the group.  The most general
// transformation combining rotations and boosts is described by the
// LorentzTransform group defined in TFourVector.h.  All angles are
// assumed to be in radians.
//
// This package was developed at the University of Connecticut by
// Richard T. Jones
//
//////////////////////////////////////////////////////////////////////////
 
#include <iostream>
using namespace std;

#include "TThreeVectorComplex.h"
#include "TThreeRotation.h"
#include "TLorentzBoost.h"
#include "TDiracSpinor.h"
#include "TDiracMatrix.h"

ClassImp(TDiracMatrix)


LDouble_t TDiracMatrix::fResolution = 1e-12;

TDiracMatrix::TDiracMatrix(const EDiracIndex i)
{
// Initializes a standard Dirac matrix from the following list:
//    i=kDiracOne    : 1
//    i=kDiracGamma0 : gamma_0
//    i=kDiracGamma1 : gamma_1
//    i=kDiracGamma2 : gamma_2
//    i=kDiracGamma3 : gamma_3
//    i=kDiracGamma4 : gamma_4 (=gamma_0)
//    i=kDiracGamma5 : gamma_5 (=i_*gamma_0*gamma_1*gamma_2*gamma_3)
//    i=kDiracSigma1 : Sigma_1
//    i=kDiracSigma2 : Sigma_2
//    i=kDiracSigma3 : Sigma_3
//    i=kDiracKappa1 : Kappa_1
//    i=kDiracKappa2 : Kappa_2
//    i=kDiracKappa3 : Kappa_3

   const Complex_t i_(0,1);
   Zero();
   switch (i) {
   case kDiracOne:
    fMatrix[0][0] = 1;
            fMatrix[1][1] = 1;
                    fMatrix[2][2] = 1;
                            fMatrix[3][3] = 1;
      break;
   case kDiracGamma0:
    fMatrix[0][0] = 1;
            fMatrix[1][1] = 1;
                    fMatrix[2][2] = -1;
                            fMatrix[3][3] = -1;
      break;
   case kDiracGamma1:
                            fMatrix[0][3] = 1;
                    fMatrix[1][2] = 1;
            fMatrix[2][1] = -1;
    fMatrix[3][0] = -1;
      break;
   case kDiracGamma2:
                            fMatrix[0][3] = -i_;
                    fMatrix[1][2] = i_;
            fMatrix[2][1] = i_;
    fMatrix[3][0] = -i_;
      break;
   case kDiracGamma3:
                    fMatrix[0][2] = 1;
                            fMatrix[1][3] = -1;
    fMatrix[2][0] = -1;
            fMatrix[3][1] = 1;
      break;
   case kDiracGamma5:
                    fMatrix[0][2] = 1;
                            fMatrix[1][3] = 1;
    fMatrix[2][0] = 1;
            fMatrix[3][1] = 1;
      break;
   case kDiracSigma1:
            fMatrix[0][1] = 1;
    fMatrix[1][0] = 1;
                            fMatrix[2][3] = 1;
                    fMatrix[3][2] = 1;
      break;
   case kDiracSigma2:
            fMatrix[0][1] = -i_;
    fMatrix[1][0] = i_;
                            fMatrix[2][3] = -i_;
                    fMatrix[3][2] = i_;
      break;
   case kDiracSigma3:
    fMatrix[0][0] = 1;
            fMatrix[1][1] = -1;
                    fMatrix[2][2] = 1;
                            fMatrix[3][3] = -1;
      break;
   case kDiracKappa1:
                            fMatrix[0][3] = i_;
                    fMatrix[1][2] = i_;
            fMatrix[2][1] = i_;
    fMatrix[3][0] = i_;
      break;
   case kDiracKappa2:
                            fMatrix[0][3] = 1;
                    fMatrix[1][2] = -1;
            fMatrix[2][1] = 1;
    fMatrix[3][0] = -1;
      break;
   case kDiracKappa3:
                    fMatrix[0][2] = i_;
                            fMatrix[1][3] = -i_;
    fMatrix[2][0] = i_;
            fMatrix[3][1] = -i_;
      break;
   default:
      break;
   }
}

TDiracMatrix::TDiracMatrix(const EDiracIndex i, const EDiracIndex j)
{
// Initializes the commutator of standard Dirac matrices i and j, as
//             i_/2 [gamma_i,gamma_j]
// where i,j come from the list of values for EDiracIndex.

   const Complex_t i_(0,1);
   const TDiracMatrix a(i),b(j);
   *this = a*b - b*a;
   *this *= i_/2.L;
}

TDiracMatrix &TDiracMatrix::Transpose()
{
   for (Int_t i=0; i<3; i++)
      for (Int_t j=i+1; j<4; j++) {
         Complex_t temp = fMatrix[i][j];
         fMatrix[i][j] = fMatrix[j][i];
         fMatrix[j][i] = temp;
      }
   return *this;
}

Complex_t TDiracMatrix::Determ() const
{
   TDeterminor deter(4);
   return deter.Minor((const Complex_t *)&fMatrix,3);
}

TDiracMatrix &TDiracMatrix::Invert()
{
   if (Determ() == 0.L) {
      Error("TDiracMatrix::Invert()","matrix is singular!");
      return *this;
   }
   TInvertor obj(4);
   obj.Invert((Complex_t *)&fMatrix);
   return *this;
}

TDiracMatrix &TDiracMatrix::SetUUbar
             (const TFourVectorReal &p, const TDiracSpinor &u)
{
   const TDiracMatrix gamma0(kDiracGamma0);
   for (Int_t i=0; i<4; i++)
      for (Int_t j=0; j<4; j++)
         fMatrix[i][j] = u[i]*conj(u[j]);
   return (*this *= gamma0);
}

TDiracMatrix &TDiracMatrix::SetVVbar
             (const TFourVectorReal &p, const TDiracSpinor &v)
{
   return SetUUbar(p,v);
}

TDiracMatrix &TDiracMatrix::SetUUbar
             (const TFourVectorReal &p, const TThreeVectorReal &polar)
{
   LDouble_t degree=polar.Length();
   if (degree < polar.Resolution())
      SetUUbar(p);
   else {
      TDiracSpinor u1, u2;
      u1.SetStateU(p,polar);
      u2.SetStateU(p,-polar);
      TDiracMatrix m1,m2;
      *this = (1+degree)*m1.SetUUbar(p,u1) + (1-degree)*m2.SetUUbar(p,u2);
   }
   return (*this /= Complex_t(2));
}

TDiracMatrix &TDiracMatrix::SetVVbar
             (const TFourVectorReal &p, const TThreeVectorReal &polar)
{
   LDouble_t degree=polar.Length();
   if (degree < polar.Resolution())
      SetVVbar(p);
   else {
      TDiracSpinor v1, v2;
      v1.SetStateV(p,polar);
      v2.SetStateV(p,-polar);
      TDiracMatrix m1,m2;
      *this = (1+degree)*m1.SetVVbar(p,v1) + (1-degree)*m2.SetVVbar(p,v2);
   }
   return (*this /= Complex_t(2));
}

TDiracMatrix &TDiracMatrix::SetRotation
             (const TUnitVector &axis, const LDouble_t angle)
{
   const Complex_t i_(0,1);
   TThreeVectorReal rotator(axis);
   LDouble_t cosHalfAngle = cos(angle/2);
   LDouble_t sinHalfAngle = sin(angle/2);
   rotator.Normalize(sinHalfAngle);
   *this = TDiracMatrix(cosHalfAngle);
   *this += i_*rotator[1]*TDiracMatrix(kDiracSigma1);
   *this += i_*rotator[2]*TDiracMatrix(kDiracSigma2);
   *this += i_*rotator[3]*TDiracMatrix(kDiracSigma3);
   return *this;
}

TDiracMatrix &TDiracMatrix::SetRotation(const TThreeRotation &rotOp)
{
   TThreeVectorReal axis=rotOp.Axis();
   return SetRotation(axis);
}

TDiracMatrix &TDiracMatrix::SetRotation(const TThreeVectorReal &axis)
{
   LDouble_t angle = axis.Length();
   return SetRotation(axis,angle);
}

TDiracMatrix &TDiracMatrix::SetRotation(const LDouble_t phi,
                                        const LDouble_t theta,
                                        const LDouble_t psi)
{
   TThreeRotation rotOp(phi,theta,psi);
   return SetRotation(rotOp);
}

TDiracMatrix &TDiracMatrix::SetBoost(const LDouble_t betaX,
                                     const LDouble_t betaY,
                                     const LDouble_t betaZ)
{
   TThreeVectorReal beta(betaX,betaY,betaZ);
   return SetBoost(beta);
}

TDiracMatrix &TDiracMatrix::SetBoost(const TThreeVectorReal &beta)
{
   return SetBoost(beta,beta.Length());
}

TDiracMatrix &TDiracMatrix::SetBoost(const TLorentzBoost &boostOp)
{
   TThreeVectorReal beta=boostOp.Beta();
   return SetBoost(beta);
}

TDiracMatrix &TDiracMatrix::SetBoost
             (const TUnitVector &bhat, const LDouble_t &beta)
{
   const Complex_t i_(0,1);
   if (abs(beta) < fResolution) {
      *this = TDiracMatrix(1);
   }
   else
   {
      LDouble_t eta = asinh(beta/sqrt(1-beta*beta));
      TThreeVectorReal booster(bhat);
      LDouble_t coshHalfEta = cosh(eta/2);
      LDouble_t sinhHalfEta = sinh(eta/2);
      booster.Normalize(sinhHalfEta);
      *this = TDiracMatrix(coshHalfEta);
      *this += i_*booster[1]*TDiracMatrix(kDiracKappa1);
      *this += i_*booster[2]*TDiracMatrix(kDiracKappa2);
      *this += i_*booster[3]*TDiracMatrix(kDiracKappa3);
   }
   return *this;
}

TDiracMatrix &TDiracMatrix::SetTransform(const TLorentzTransform &xform)
{
   TLorentzBoost boostOp;
   TThreeRotation rotOp;
   TThreeVectorReal beta,axis;
   xform.Factorize(boostOp,rotOp);
   beta = boostOp.Beta();
   axis = rotOp.Axis();
   TDiracMatrix rotMat;
   SetBoost(beta);
   return (*this *= rotMat.SetRotation(axis));
}

TDiracMatrix &TDiracMatrix::SimTransform(const TDiracMatrix &m)
{
// Similarity transform is defined as A' = M A Minverse

   TDiracMatrix oldOp(*this);
   *this = m;
   *this *= oldOp;
   return (*this /= m);
}

TDiracMatrix &TDiracMatrix::UniTransform(const TDiracMatrix &m)
{
// Unitary transform is defined as A' = M A Mdagger

   TDiracMatrix mtemp(m);
   TDiracMatrix oldOp(*this);
   *this = m;
   *this *= oldOp;
   return (*this *= mtemp.Adjoint());
}

void TDiracMatrix::Streamer(TBuffer &buf)
{
   // Put/get a Dirac matrix to/from stream buffer buf.

   Double_t vector[32];
   if (buf.IsReading()) {
      buf.ReadStaticArray(vector);
      int n=0;
      for (int mu=0; mu < 4; ++mu) {
         for (int nu=0; nu < 4; ++nu, ++n) {
            fMatrix[mu][nu] = Complex_t(vector[n], vector[n+1]);
         }
      }
   } else {
      int n=0;
      for (int mu=0; mu < 4; ++mu) {
         for (int nu=0; nu < 4; ++nu, ++n) {
            vector[n] = fMatrix[mu][nu].real();
            vector[n+1] = fMatrix[mu][nu].imag();
         }
      }
      buf.WriteArray(vector, 32);
   }
}

void TDiracMatrix::Print(Option_t *option)
{
   // Output the Dirac matrix in ascii form.

   cout << "TDiracMatrix is" << endl;
   cout << "(" << fMatrix[0][0] << "   " << fMatrix[0][1] << "   "
               << fMatrix[0][2] << "   " << fMatrix[0][3] << ")" << endl;
   cout << "(" << fMatrix[1][0] << "   " << fMatrix[1][1] << "   "
               << fMatrix[1][2] << "   " << fMatrix[1][3] << ")" << endl;
   cout << "(" << fMatrix[2][0] << "   " << fMatrix[2][1] << "   "
               << fMatrix[2][2] << "   " << fMatrix[2][3] << ")" << endl;
   cout << "(" << fMatrix[3][0] << "   " << fMatrix[3][1] << "   "
               << fMatrix[3][2] << "   " << fMatrix[3][3] << ")" << endl;
}
 
//______________________________________________________________________________


#ifdef R__HPUX

//______________________________________________________________________________
//  These functions should be inline
//______________________________________________________________________________

#endif
