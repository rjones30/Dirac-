//
// TLorentzTransform.h 
//
// This file is distributed as part of the Lorentz++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See LorentzPackage.h for details.

#ifndef ROOT_TLorentzTransform
#define ROOT_TLorentzTransform
 
#include "Double.h"
#include "TBuffer.h"
#include "Complex.h"
#include "TFourVectorComplex.h"
#include "TError.h"
 
#include <math.h>

class TLorentzBoost;
class TThreeRotation;

 
class TLorentzTransform {
 
protected:

   LDouble_t        fMatrix[4][4];    // storage for rotation matrix;
   static LDouble_t fResolution;      // matrix resolving "distance"

   LDouble_t Determ() const;

public:
   TLorentzTransform() { }
   TLorentzTransform(const TLorentzTransform &another);
 
   virtual ~TLorentzTransform() { }
 
   static void SetResolution(const LDouble_t resolution);
   LDouble_t Resolution() const;

   Bool_t IsNull();
   Bool_t IsRotation();
   Bool_t IsLorentzBoost();
   Bool_t IsOrthogonal();
   Bool_t IsIsochronous();
   Bool_t IsProper();

   void Factorize(TLorentzBoost &boost, TThreeRotation &rot) const;

   TLorentzTransform &Null();
   TLorentzTransform &TimeRev();
   TLorentzTransform &SpaceInv();
   TLorentzTransform &Transpose();
   TLorentzTransform &Invert();
 
   TLorentzTransform &operator=(const TLorentzTransform &source);
   TLorentzTransform &operator*=(const TLorentzTransform &source);
 
   Bool_t operator==(const TLorentzTransform &other) const;
   Bool_t operator!=(const TLorentzTransform &other) const;
 
   TFourVectorReal operator*(const TFourVectorReal &vec) const;
   TFourVectorComplex operator*(const TFourVectorComplex &vec) const;
   TLorentzTransform operator*(const TLorentzTransform &xform) const;

   friend TBuffer &operator>>(TBuffer &buf, TLorentzTransform *&xOp);
   friend TBuffer &operator<<(TBuffer &buf, const TLorentzTransform *xOp);
   void Print(Option_t *option="");
 
   ClassDef(TLorentzTransform,1)  // Lorentz transform operator class
};

class TDeterminor {

// The TDeterminor class is a helper for calculating the determinant of
// a square matrix.  It uses recursive procedure to cover any size of
// square matrix, using a permutation iteration method.  Matrices of
// Float_t, LDouble_t and Complex_t elements are currently supported.

private:
   Int_t fDim;        // dimension of square matrix
   Int_t *fRow;       // permutator array of colIndex->rowIndex

public:
   TDeterminor() : fDim(0), fRow(0) { }
   TDeterminor(const Int_t dim);
   ~TDeterminor();
   Float_t Minor(const Float_t *matrix, const Int_t ncol);
   LDouble_t Minor(const LDouble_t *matrix, const Int_t ncol);
   Complex_t Minor(const Complex_t *matrix, const Int_t ncol);
   void Swap(Int_t &a, Int_t &b);
};


class TInvertor {

// The TInvertor class is a helper for calculating the inverse of a
// square matrix.  It uses a pivoting (i.e. Gaussian elimination) method.
// Matrices of Float_t, LDouble_t and Complex_t elements are supported.

private:
   Int_t fDim;            // dimension of square matrix
   Int_t *fPivot;        // array of pivot element indices

public:
   TInvertor() : fDim(0), fPivot(0) { }
   TInvertor(const Int_t dim);
   ~TInvertor();
   Float_t *Invert(Float_t *matrix);
   Float_t *Invert(const Float_t *matrix, Float_t *inverse);
   LDouble_t *Invert(LDouble_t *matrix);
   LDouble_t *Invert(const LDouble_t *matrix, LDouble_t *inverse);
   Complex_t *Invert(Complex_t *matrix);
   Complex_t *Invert(const Complex_t *matrix, Complex_t *inverse);
   Int_t SetPivot(const Int_t row, const Float_t *matrix);
   Int_t SetPivot(const Int_t row, const LDouble_t *matrix);
   Int_t SetPivot(const Int_t row, const Complex_t *matrix);
   void Swap(Int_t &a, Int_t &b);
   void PivotRow
        (Int_t row1, Int_t row2, Float_t *matrix, Float_t *inverse);
   void PivotRow
        (Int_t row1, Int_t row2, LDouble_t *matrix, LDouble_t *inverse);
   void PivotRow
        (Int_t row1, Int_t row2, Complex_t *matrix, Complex_t *inverse);
};

//----- inlines ----------------------------------------------------------------

inline TLorentzTransform::TLorentzTransform(const TLorentzTransform &another)
{
   *this = another;
}
 
inline void TLorentzTransform::SetResolution(const LDouble_t resolution)
{
   fResolution = resolution;
}

inline LDouble_t TLorentzTransform::Resolution() const
{
   LDouble_t scale = ((fMatrix[0][0] > 0) ? fMatrix[0][0] : -fMatrix[0][0]);
   if (scale > 1)
      return fResolution*scale;
   else
      return fResolution;
}
 
inline TLorentzTransform &TLorentzTransform::Null()
{
   LDouble_t *p = fMatrix[0];
   for (Int_t i=0; i<4; i++)
      for (Int_t j=0; j<4; j++)
         *(p++) = (i==j ? 1 : 0);
   return *this;
}

inline TLorentzTransform &TLorentzTransform::TimeRev()
{
   LDouble_t *pRow = fMatrix[0];
   for (Int_t i=0; i<4; i++) *(pRow++) *= -1;
   return *this;
}

inline TLorentzTransform &TLorentzTransform::SpaceInv()
{
   LDouble_t *pRow = fMatrix[1];
   for (Int_t i=0; i<12; i++) *(pRow++) *= -1;
   return *this;
}

inline TLorentzTransform &TLorentzTransform::Transpose()
{
   TLorentzTransform copy(*this);
   for (Int_t i=0; i<4; i++)
      for (Int_t j=0; j<4; j++)
         fMatrix[i][j] = copy.fMatrix[j][i];
   return *this; 
}

inline TLorentzTransform &TLorentzTransform::Invert()
{
   SpaceInv();
   Transpose();
   SpaceInv();
   return *this;
}

inline TLorentzTransform &TLorentzTransform::operator=
                         (const TLorentzTransform &source)
{
   memcpy(fMatrix, source.fMatrix, 16*sizeof(LDouble_t));
   return *this;
}

inline Bool_t TLorentzTransform::operator!=
              (const TLorentzTransform &other) const
{
   if (*this == other) return 0;
   return 1;
}

inline TBuffer &operator>>(TBuffer &buf, TLorentzTransform *&xOp)
{
   Double_t matrix[4][4];
   buf.ReadStaticArray(&matrix[0][0]);
   for (int mu=0; mu < 4; ++mu)
      for (int nu=0; nu < 4; ++nu)
         xOp->fMatrix[mu][nu] = matrix[mu][nu];
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TLorentzTransform *xOp)
{
   Double_t matrix[4][4];
   for (int mu=0; mu < 4; ++mu)
      for (int nu=0; nu < 4; ++nu)
         matrix[mu][nu] = xOp->fMatrix[mu][nu];
   buf.WriteArray(&matrix[0][0], 16);
   return buf;
}

inline LDouble_t TLorentzTransform::Determ() const
{
   TDeterminor deter(4);
   return deter.Minor((const LDouble_t *)&fMatrix,3);
}

inline TDeterminor::TDeterminor(const Int_t dim) : fDim(dim)
{
   fRow = new Int_t[dim];
   for (Int_t i=0; i<dim; i++)
      fRow[i] = i;
}

inline TDeterminor::~TDeterminor() { delete [] fRow; }

inline void TDeterminor::Swap(Int_t &a, Int_t &b)
{
   Int_t temp=a; a=b; b=temp;
}

#endif
