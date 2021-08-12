//
// TDiracSpinor.cxx
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
// type TDiracIndex.  A TDiracIndex can take on a value from the list
//    kDiracOne,    kDiracGamma1, kDiracGamma2, kDiracGamma3,
//    kDiracGamma4, kDiracGamma5, kDiracSigma1, kDiracSigma2,
//    kDiracSigma3, kDiracKappa1, kDiracKappa2, kDiracKappa3.
// The constructor invoked with two TDiracIndex values i,j returns
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

#include "TFourVectorComplex.h"
#include "TLorentzBoost.h"
#include "TThreeRotation.h"
#include "TPauliSpinor.h"
#include "TDiracSpinor.h"
#include "TDiracMatrix.h"

ClassImp(TDiracSpinor)


LDouble_t TDiracSpinor::fResolution = 1e-12;

Complex_t TDiracSpinor::ScalarProd(const TDiracSpinor &other)
{
   const TDiracMatrix gamma0(kDiracGamma0);
   return InnerProd(gamma0 * other);
}
 
TPauliSpinor TDiracSpinor::Upper() const
{
   TPauliSpinor upper;
   upper[0] = fSpinor[0];
   upper[1] = fSpinor[1];
   return upper;
}

TPauliSpinor TDiracSpinor::Lower() const
{
   TPauliSpinor lower;
   lower[0] = fSpinor[2];
   lower[1] = fSpinor[3];
   return lower;
}
 
TDiracSpinor TDiracSpinor::SetUpper(TPauliSpinor &phi)
{
   fSpinor[0] = phi[0];
   fSpinor[1] = phi[1];
   return *this;
}

TDiracSpinor TDiracSpinor::SetLower(TPauliSpinor &chi)
{
   fSpinor[2] = chi[0];
   fSpinor[3] = chi[1];
   return *this;
}

TDiracSpinor &TDiracSpinor::SetStateU(const TFourVectorReal &p,
                                      const Float_t helicity)
{
   // U is the positive-energy solution to the Dirac equation with momentum
   // p and given helicity.  It represents a fermion of momentum p and good
   // helicity.  The axis of spin quantization is the momentum direction
   // unless the particle is at rest, in which case it is the +z axis.

   Int_t h = (helicity < 0 ? -1 : +1);
   TThreeVectorReal quantax(p);
   if (quantax.Length() < quantax.Resolution()) {
      quantax[1]=0;
      quantax[2]=0;
      quantax[3]=1;
   }
   if (p[0] > 0) {
      TPauliSpinor phi(quantax*h);
      TPauliSpinor chi(phi);
      chi *= h*p.Length()/(p[0]+p.Invariant());
      SetUpper(phi);
      SetLower(chi);
      return Normalize(p);
   }
   else if (p[0] < 0) {
      TPauliSpinor chi(quantax*h);
      TPauliSpinor phi(chi);
      phi *= h*p.Length()/(p[0]-p.Invariant());
      SetUpper(phi);
      SetLower(chi);
      return Normalize(p);
   }
   else
      return Zero();
} 

TDiracSpinor &TDiracSpinor::SetStateV(const TFourVectorReal &p,
                                      const Float_t helicity)
{
   // V is the negative-energy solution to the Dirac equation with momentum
   // -p and given helicity.  It represents an antifermion of momentum p and
   // good helicity.  The axis of spin quantization is the momentum direction
   // unless the particle is at rest, in which case it is the +z axis.

   return SetStateU(-p,helicity);
} 

TDiracSpinor &TDiracSpinor::SetStateU(const TFourVectorReal &p,
                                      const TUnitVector &polar)
{
   // Helicity-frame polarization is provided in argument polar, otherwise 
   // the same as SetStateU(TFourVectorReal &p, Float_t helicity).

   TPauliSpinor phi(polar);
   TDiracSpinor u1(p,+1);
   TDiracSpinor u2(p,-1);
   *this  = u1 *= phi[0];
   *this += u2 *= phi[1];
   return *this;
}

TDiracSpinor &TDiracSpinor::SetStateV(const TFourVectorReal &p,
                                      const TUnitVector &polar)
{
   // Helicity-frame polarization is provided in argument polar, otherwise 
   // the same as SetStateV(TFourVectorReal &p, Float_t helicity).

   return SetStateU(-p,polar);
}
 
TDiracSpinor &TDiracSpinor::Operate(const TDiracMatrix &dmOp)
{
   TDiracSpinor temp(*this);
   fSpinor[0] = dmOp.fMatrix[0][0]*temp[0] + dmOp.fMatrix[0][1]*temp[1] +
                dmOp.fMatrix[0][2]*temp[2] + dmOp.fMatrix[0][3]*temp[3] ;
   fSpinor[1] = dmOp.fMatrix[1][0]*temp[0] + dmOp.fMatrix[1][1]*temp[1] +
                dmOp.fMatrix[1][2]*temp[2] + dmOp.fMatrix[1][3]*temp[3] ;
   fSpinor[2] = dmOp.fMatrix[2][0]*temp[0] + dmOp.fMatrix[2][1]*temp[1] +
                dmOp.fMatrix[2][2]*temp[2] + dmOp.fMatrix[2][3]*temp[3] ;
   fSpinor[3] = dmOp.fMatrix[3][0]*temp[0] + dmOp.fMatrix[3][1]*temp[1] +
                dmOp.fMatrix[3][2]*temp[2] + dmOp.fMatrix[3][3]*temp[3] ;
   return *this;
}

TDiracSpinor &TDiracSpinor::Rotate(const TThreeRotation &rotOp)
{
   TDiracMatrix dmR;
   dmR.SetRotation(rotOp);
   return Operate(dmR);
}

TDiracSpinor &TDiracSpinor::Rotate(const LDouble_t &phi,
                                   const LDouble_t &theta,
                                   const LDouble_t &psi)
{
   TDiracMatrix dmR;
   dmR.SetRotation(phi,theta,psi);
   return Operate(dmR);
}

TDiracSpinor &TDiracSpinor::Rotate(const TThreeVectorReal &axis)
{
   TDiracMatrix dmR;
   dmR.SetRotation(axis);
   return Operate(dmR);
}

TDiracSpinor &TDiracSpinor::Rotate
             (const TUnitVector &axis, const LDouble_t angle)
{
   TDiracMatrix dmR;
   dmR.SetRotation(axis,angle);
   return Operate(dmR);
}

TDiracSpinor &TDiracSpinor::Boost(const TLorentzBoost &boostOp)
{
   TDiracMatrix dmB;
   dmB.SetBoost(boostOp);
   return Operate(dmB);
}

TDiracSpinor &TDiracSpinor::Boost(const LDouble_t betaX,
                                  const LDouble_t betaY,
                                  const LDouble_t betaZ)
{
   TDiracMatrix dmB;
   dmB.SetBoost(betaX,betaY,betaZ);
   return Operate(dmB);
}

TDiracSpinor &TDiracSpinor::Boost(const LDouble_t *beta)
{
   TDiracMatrix dmB;
   dmB.SetBoost((TThreeVectorReal)beta);
   return Operate(dmB);
}

TDiracSpinor &TDiracSpinor::Boost(const TThreeVectorReal &beta)
{
   TDiracMatrix dmB;
   dmB.SetBoost(beta);
   return Operate(dmB);
}

TDiracSpinor &TDiracSpinor::Boost
             (const TUnitVector &bhat, const LDouble_t beta)
{
   TDiracMatrix dmB;
   dmB.SetBoost(bhat,beta);
   return Operate(dmB);
}

TDiracSpinor &TDiracSpinor::BoostToRest(const TFourVector &p)
{
   TDiracMatrix boost;
   boost.SetBoost(p/p[0]);
   return Operate(boost);
}

TDiracSpinor &TDiracSpinor::BoostFromRest(const TFourVector &p)
{
   TDiracMatrix boost;
   boost.SetBoost(-p/p[0]);
   return Operate(boost);
}

TDiracSpinor operator*(const TDiracMatrix &dmOp, const TDiracSpinor &vec)
{
   TDiracSpinor result(vec);
   return result.Operate(dmOp);
}

void TDiracSpinor::Streamer(TBuffer &buf)
{
   // Put/get a Dirac spinor to/from stream buffer buf.

   Double_t vector[8];
   if (buf.IsReading()) {
      buf.ReadStaticArray(vector);
      fSpinor[0] = Complex_t(vector[0], vector[1]);
      fSpinor[1] = Complex_t(vector[2], vector[3]);
      fSpinor[2] = Complex_t(vector[4], vector[5]);
      fSpinor[3] = Complex_t(vector[6], vector[7]);
   } else {
      vector[0] = fSpinor[0].real();
      vector[1] = fSpinor[0].imag();
      vector[2] = fSpinor[1].real();
      vector[3] = fSpinor[1].imag();
      vector[4] = fSpinor[2].real();
      vector[5] = fSpinor[2].imag();
      vector[6] = fSpinor[3].real();
      vector[7] = fSpinor[3].imag();
      buf.WriteArray(vector, 8);
   }
}

void TDiracSpinor::Print(Option_t *option)
{
   // Output an ascii representation for the Dirac spinor.

   cout << "TDiracSpinor is {" << endl << fSpinor[0] << endl
        << fSpinor[1] << endl << fSpinor[2] << endl << fSpinor[3]
        << endl << "}" << endl;
}

 
//______________________________________________________________________________


#ifdef R__HPUX

//______________________________________________________________________________
//  These functions should be inline
//______________________________________________________________________________

#endif
