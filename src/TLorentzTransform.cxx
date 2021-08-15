// 
// TLorentzTransform.cxx
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

#include "TLorentzTransform.h"
#include "TLorentzBoost.h"
#include "TThreeRotation.h"
#include "TFourVectorReal.h"
#include "TFourVectorComplex.h"

#include <iostream>
using namespace std;

ClassImp(TLorentzTransform)


LDouble_t TLorentzTransform::fResolution = 1e-12;

Bool_t TLorentzTransform::IsNull()
{
   LDouble_t *p = fMatrix[0];
   for (Int_t i=0; i<4; i++) {
      for (Int_t j=0; j<4; j++) {
         LDouble_t dif = *(p++) - (i==j ? 1 : 0);
         if (abs(dif) > Resolution()) return 0;
      }
   }
   return 1;
}

Bool_t TLorentzTransform::IsRotation()
{
   TThreeRotation rot;
   TLorentzBoost boost;
   Factorize(boost,rot);
   return boost.IsNull();
}

Bool_t TLorentzTransform::IsLorentzBoost()
{
   TThreeRotation rot;
   TLorentzBoost boost;
   Factorize(boost,rot);
   return rot.IsNull();
}

Bool_t TLorentzTransform::IsOrthogonal()
{
   TLorentzTransform copy(*this);
   TLorentzTransform copyT(*this);
   copy.SpaceInv();
   copyT.Transpose();
   copyT.SpaceInv();
   copy *= copyT;
   return copy.IsNull();
}

Bool_t TLorentzTransform::IsIsochronous()
{
   return (fMatrix[0][0] > 0);
}

Bool_t TLorentzTransform::IsProper()
{
   LDouble_t dif = Determ() - 1;
   return (abs(dif) < Resolution());
}

void TLorentzTransform::Factorize(TLorentzBoost &boost, TThreeRotation &rot) const
{
   TThreeVectorReal beta((LDouble_t *)&fMatrix[0][1]);
   beta /= -fMatrix[0][0];
   boost.SetBeta(beta);
   TLorentzBoost xOp(-beta);
   xOp *= *this;
   (TLorentzTransform &)rot = xOp;
}

TFourVectorReal TLorentzTransform::operator*
                       (const TFourVectorReal &vec) const
{
   TFourVectorReal result;
   for (Int_t i=0; i<4; i++) {
      LDouble_t sum=0;
      for (Int_t j=0; j<4; j++) {
         sum += fMatrix[i][j]*vec.fVector[j];
      }
      result.fVector[i] = sum;
   }
   return result;
}

TFourVectorComplex TLorentzTransform::operator*
                          (const TFourVectorComplex &vec) const
{
   TFourVectorComplex result;
   for (Int_t i=0; i<4; i++) {
      Complex_t sum=0;
      for (Int_t j=0; j<4; j++) {
         sum += fMatrix[i][j]*vec.fVector[j];
      }
      result.fVector[i] = sum;
   }
   return result;
}

TLorentzTransform &TLorentzTransform::operator*=
                         (const TLorentzTransform &source)
{
   TLorentzTransform temp(*this);
   for (Int_t i=0; i<4; i++) {
      for (Int_t j=0; j<4; j++) {
         LDouble_t sum=0;
         for (Int_t k=0; k<4; k++) {
            sum += temp.fMatrix[i][k]*source.fMatrix[k][j];
         }
         fMatrix[i][j] = sum;
      }
   }
   return *this;
}
 
TLorentzTransform TLorentzTransform::operator*
                  (const TLorentzTransform &xform) const
{
   TLorentzTransform result;
   for (Int_t i=0; i<4; i++) {
      for (Int_t j=0; j<4; j++) {
         LDouble_t sum=0;
         for (Int_t k=0; k<4; k++) {
            sum += fMatrix[i][k]*xform.fMatrix[k][j];
         }
         result.fMatrix[i][j] = sum;
      }
   }
   return result;
}

Bool_t TLorentzTransform::operator==
              (const TLorentzTransform &other) const
{
   LDouble_t *pThis = (LDouble_t *)&fMatrix[0][0];
   LDouble_t *pOther = (LDouble_t *)&other.fMatrix[0][0];
   for (Int_t i=0; i< 16; i++, pThis++, pOther++) {
      LDouble_t dif = *pThis - *pOther;
      if (abs(dif) > Resolution()) return 0;
   }
   return 1;
}

void TLorentzTransform::Streamer(TBuffer &buf)
{
   // Put/get a Lorentz transform matrix to/from stream buffer buf.

   Double_t matrix[4][4];
   if (buf.IsReading()) {
      buf.ReadStaticArray(&matrix[0][0]);
      for (int mu=0; mu < 4; ++mu)
         for (int nu=0; nu < 4; ++nu)
            fMatrix[mu][nu] = matrix[mu][nu];
   } else {
      for (int mu=0; mu < 4; ++mu)
         for (int nu=0; nu < 4; ++nu)
            matrix[mu][nu] = fMatrix[mu][nu];
      buf.WriteArray(&matrix[0][0], 16);
   }
}

void TLorentzTransform::Print(Option_t *option)
{
   // Output an ascii representation for Lorentz transform matrix.

   cout << "TLorentzTransform matrix" << endl;
   cout << "(" << fMatrix[0][0] << "," << fMatrix[0][1] << ","
               << fMatrix[0][2] << "," << fMatrix[0][3] << ")" << endl;
   cout << "(" << fMatrix[1][0] << "," << fMatrix[1][1] << ","
               << fMatrix[1][2] << "," << fMatrix[1][3] << ")" << endl;
   cout << "(" << fMatrix[2][0] << "," << fMatrix[2][1] << ","
               << fMatrix[2][2] << "," << fMatrix[2][3] << ")" << endl;
   cout << "(" << fMatrix[3][0] << "," << fMatrix[3][1] << ","
               << fMatrix[3][2] << "," << fMatrix[3][3] << ")" << endl;
}
 
Float_t TDeterminor::Minor(const Float_t *matrix, const Int_t ncol)
{
   // Use a recursive method to evaluate the determinant of matrix
   // with ncol columns.  The matrix can be stored row-wise or column-wise.

   if (ncol == 0) return matrix[fRow[0]];
   Float_t result=0;
   for (Int_t col=0; col<ncol; col++) {
      Swap(fRow[col],fRow[ncol]);
      result -= matrix[fRow[ncol]+ncol*fDim]*Minor(matrix,ncol-1);
      Swap(fRow[col],fRow[ncol]);
   }
   result += matrix[fRow[ncol]+ncol*fDim]*Minor(matrix,ncol-1);
   return result;
}

LDouble_t TDeterminor::Minor(const LDouble_t *matrix, const Int_t ncol)
{
   // Use a recursive method to evaluate the determinant of matrix
   // with ncol columns.  The matrix can be stored row-wise or column-wise.

   if (ncol == 0) return matrix[fRow[0]];
   LDouble_t result=0;
   for (Int_t col=0; col<ncol; col++) {
      Swap(fRow[col],fRow[ncol]);
      result -= matrix[fRow[ncol]+ncol*fDim]*Minor(matrix,ncol-1);
      Swap(fRow[col],fRow[ncol]);
   }
   result += matrix[fRow[ncol]+ncol*fDim]*Minor(matrix,ncol-1);
   return result;
}

Complex_t TDeterminor::Minor(const Complex_t *matrix, const Int_t ncol)
{
   // Use a recursive method to evaluate the determinant of matrix
   // with ncol columns.  The matrix can be stored row-wise or column-wise.

   if (ncol == 0) return matrix[fRow[0]];
   Complex_t result=0;
   for (Int_t col=0; col<ncol; col++) {
      Swap(fRow[col],fRow[ncol]);
      result -= matrix[fRow[ncol]+ncol*fDim]*Minor(matrix,ncol-1);
      Swap(fRow[col],fRow[ncol]);
   }
   result += matrix[fRow[ncol]+ncol*fDim]*Minor(matrix,ncol-1);
   return result;
}

TInvertor::TInvertor(const Int_t dim)
{
   if (dim <= 0) {
      Error("TInvertor::TInvertor","dimension nonpositive");
   }
   fDim = dim;
   fPivot = new Int_t[dim];
   for (Int_t i=0; i<dim; i++)
      fPivot[i] = i;
}

TInvertor::~TInvertor()
{
   delete [] fPivot;
}

Float_t *TInvertor::Invert(Float_t *matrix)
{
   return Invert(matrix,matrix);
}

Float_t *TInvertor::Invert(const Float_t *matrix, Float_t *inverse)
{
   // Invert matrix into inverse using a pivoting method. 

   const Int_t nelem = fDim*fDim;
   Float_t *work = new Float_t[nelem];
   Float_t *winv = new Float_t[nelem];
   const Float_t *m = matrix;
   Float_t *w = work;
   Float_t *v = winv;
   for (Int_t i=0; i<nelem; i++) {             // copy initial matrix
      *(w++) = *(m++);                         // to work matrix
      *(v++) = 0;                              // and zero work inverse
   }
   v = winv;                                   // initialize inverse work
   for (Int_t i=0; i<fDim; i++, v+=fDim+1)     // matrix to unity
      *v = 1;

   for (Int_t row=0; row<fDim; row++) {        // pivot down to one non-zero
      SetPivot(row,work);                      // element per row or column
      for (Int_t r=0; r<fDim; r++)
         if (r != row) PivotRow(row,r,work,winv);
   }

   for (Int_t row=0; row<fDim; row++) {        // copy to destination array
      Int_t r = fPivot[row];                   // swapping rows as needed
      LDouble_t norm = work[row*fDim+r];
      for (Int_t i=0; i<fDim; i++)
         inverse[r*fDim+i] = winv[row*fDim+i]/norm;
   }
   delete [] work;
   delete [] winv;
   return inverse;
}

LDouble_t *TInvertor::Invert(LDouble_t *matrix)
{
   return Invert(matrix,matrix);
}

LDouble_t *TInvertor::Invert(const LDouble_t *matrix, LDouble_t *inverse)
{
   // Invert matrix into inverse using a pivoting method. 

   const Int_t nelem = fDim*fDim;
   LDouble_t *work = new LDouble_t[nelem];
   LDouble_t *winv= new LDouble_t[nelem];
   const LDouble_t *m = matrix;
   LDouble_t *w = work;
   LDouble_t *v = winv;
   for (Int_t i=0; i<nelem; i++) {              // copy initial matrix
      *(w++) = *(m++);                          // to work matrix
      *(v++) = 0;                               // and zero work inverse
   }
   v = winv;                                    // initialize inverse work
   for (Int_t i=0; i<fDim; i++, v+=fDim+1)      // matrix to unity
      *v = 1;

   for (Int_t row=0; row<fDim; row++) {         // pivot down to one non-zero
      SetPivot(row,work);                       // element per row or column
      for (Int_t r=0; r<fDim; r++)
         if (r != row) PivotRow(row,r,work,winv);
   }

   for (Int_t row=0; row<fDim; row++) {         // copy to destination array
      Int_t r = fPivot[row];                    // swapping rows as needed
      LDouble_t norm = work[row*fDim+r];
      for (Int_t i=0; i<fDim; i++)
         inverse[r*fDim+i] = winv[row*fDim+i]/norm;
   }
   delete [] work;
   delete [] winv;
   return inverse;
}

Complex_t *TInvertor::Invert(Complex_t *matrix)
{
   return Invert(matrix,matrix);
}

Complex_t *TInvertor::Invert
                 (const Complex_t *matrix, Complex_t *inverse)
{
   // Invert matrix into inverse using a pivoting method. 

   const Int_t nelem = fDim*fDim;
   Complex_t *work = new Complex_t[nelem];
   Complex_t *winv = new Complex_t[nelem];
   const Complex_t *m = matrix;
   Complex_t *w = work;
   Complex_t *v = winv;
   for (Int_t i=0; i<nelem; i++) {                // copy initial matrix
      *(w++) = *(m++);                            // to work matrix
      *(v++) = 0;                                 // and zero work inverse
   }
   v = winv;                                      // initialize inverse work
   for (Int_t i=0; i<fDim; i++, v+=fDim+1)        // matrix to unity
      *v = 1;

   for (Int_t row=0; row<fDim; row++) {           // pivot down to one non-zero
      SetPivot(row,work);                         // element per row or column
      for (Int_t r=0; r<fDim; r++)
         if (r != row) PivotRow(row,r,work,winv);
   }

   for (Int_t row=0; row<fDim; row++) {           // copy to destination array
      Int_t r = fPivot[row];                      // swapping rows as needed
      Complex_t norm = work[row*fDim+r];
      for (Int_t i=0; i<fDim; i++)
         inverse[r*fDim+i] = winv[row*fDim+i]/norm;
   }
   delete [] work;
   delete [] winv;
   return inverse;
}

Int_t TInvertor::SetPivot(const Int_t row, const Float_t *matrix)
{
   Float_t max=0;
   Int_t rmax=-1;
   matrix += row*fDim;
   for (Int_t r=row; r<fDim; r++) {
      Int_t col = fPivot[r];
      if (abs(matrix[col]) > max) {
         max = abs(matrix[col]);
         rmax = r;
      }
   }
   if (rmax < 0) {
      Error("TInvertor::SetPivot","all row elements are zero");
      return fPivot[row];
   }
   Swap(fPivot[row],fPivot[rmax]);
   return fPivot[row];
}

Int_t TInvertor::SetPivot(const Int_t row, const LDouble_t *matrix)
{
   LDouble_t max=0;
   Int_t rmax=-1;
   matrix += row*fDim;
   for (Int_t r=row; r<fDim; r++) {
      Int_t col = fPivot[r];
      if (abs(matrix[col]) > max) {
         max = abs(matrix[col]);
         rmax = r;
      }
   }
   if (rmax < 0) {
      Error("TInvertor::SetPivot","all row elements are zero");
      return fPivot[row];
   }
   Swap(fPivot[row],fPivot[rmax]);
   return fPivot[row];
}

Int_t TInvertor::SetPivot(const Int_t row, const Complex_t *matrix)
{
   LDouble_t max=0;
   Int_t rmax=-1;
   matrix += row*fDim;
   for (Int_t r=row; r<fDim; r++) {
      Int_t col = fPivot[r];
      if (abs(matrix[col]) > max) {
         max = abs(matrix[col]);
         rmax = r;
      }
   }
   if (rmax < 0) {
      Error("TInvertor::SetPivot","all row elements are zero");
      return fPivot[row];
   }
   Swap(fPivot[row],fPivot[rmax]);
   return fPivot[row];
}

void TInvertor::Swap(Int_t &a, Int_t &b)
{
   Int_t temp=a; a=b; b=temp;
}

void TInvertor::PivotRow
            (Int_t row1, Int_t row2, Float_t *matrix, Float_t *inverse)
{
   Float_t *p1 = matrix + fDim*row1;
   Float_t *p2 = matrix + fDim*row2;
   Int_t jpivot = fPivot[row1];
   LDouble_t pfactor = p2[jpivot]/p1[jpivot];
   p2[jpivot] = 0;
   for (Int_t j=row1+1; j<fDim; j++) {
      Int_t jtarget = fPivot[j];
      p2[jtarget] -= p1[jtarget]*pfactor;
   }
   p1 = inverse + fDim*row1;
   p2 = inverse + fDim*row2;
   for (Int_t j=0; j<fDim; j++) {
      p2[j] -= p1[j]*pfactor;
   }
}

void TInvertor::PivotRow
            (Int_t row1, Int_t row2, LDouble_t *matrix, LDouble_t *inverse)
{
   LDouble_t *p1 = matrix + fDim*row1;
   LDouble_t *p2 = matrix + fDim*row2;
   Int_t jpivot = fPivot[row1];
   LDouble_t pfactor = p2[jpivot]/p1[jpivot];
   p2[jpivot] = 0;
   for (Int_t j=row1+1; j<fDim; j++) {
      Int_t jtarget = fPivot[j];
      p2[jtarget] -= p1[jtarget]*pfactor;
   }
   p1 = inverse + fDim*row1;
   p2 = inverse + fDim*row2;
   for (Int_t j=0; j<fDim; j++) {
      p2[j] -= p1[j]*pfactor;
   }
}

void TInvertor::PivotRow
            (Int_t row1, Int_t row2, Complex_t *matrix, Complex_t *inverse)
{
   Complex_t *p1 = matrix + fDim*row1;
   Complex_t *p2 = matrix + fDim*row2;
   Int_t jpivot = fPivot[row1];
   Complex_t pfactor = p2[jpivot]/p1[jpivot];
   p2[jpivot] = 0;
   for (Int_t j=row1+1; j<fDim; j++) {
      Int_t jtarget = fPivot[j];
      p2[jtarget] -= p1[jtarget]*pfactor;
   }
   p1 = inverse + fDim*row1;
   p2 = inverse + fDim*row2;
   for (Int_t j=0; j<fDim; j++) {
      p2[j] -= p1[j]*pfactor;
   }
}

//______________________________________________________________________________


#ifdef R__HPUX

//______________________________________________________________________________
//  These functions should be inline
//______________________________________________________________________________

#endif
