//
// TPauliMatrix.h 
//
// This file is distributed as part of the Pauli++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See PauliPackage.h for details.

#ifndef ROOT_TPauliMatrix
#define ROOT_TPauliMatrix
 
#include "Double.h"
#include "TThreeVectorComplex.h"
#include "TError.h"
 
#include <math.h>

// special constants for indicating standard matrix forms

enum EPauliIndex
{
   kPauliOne=0,
   kPauliSigma1=1,
   kPauliSigma2=2,
   kPauliSigma3=3
};

class TPauliMatrix {

friend class TPauliSpinor;

protected:
   Complex_t       fMatrix[2][2];   // complex matrix allocated on stack
   static LDouble_t fResolution;	    // matrix resolving "distance"

public:
   TPauliMatrix() { }
   TPauliMatrix(const EPauliIndex i);
   TPauliMatrix(const Int_t a);
   TPauliMatrix(const LDouble_t a);
   TPauliMatrix(const Complex_t &a);
   TPauliMatrix(const Complex_t &a, const TThreeVectorComplex &b);
   TPauliMatrix(const TPauliMatrix &another);
 
   virtual ~TPauliMatrix() { }

   Complex_t *operator[](Int_t row) const;

   static void SetResolution(const LDouble_t resolution);
   LDouble_t Resolution() const;
 
   Bool_t IsIdentity() const;
   Bool_t IsUnitary() const;
   Bool_t IsDiagonal() const;
   Bool_t IsHermetian() const;
   Bool_t IsIdempotent() const;

   Complex_t Trace() const;
   Complex_t Determ() const;
   void Decompose(LDouble_t &a, TThreeVectorReal &b) const;
   void Decompose(Complex_t &a, TThreeVectorComplex &b) const;
   void GetDiagonal(Complex_t &a11, Complex_t &a22) const;

   TPauliMatrix &operator=(const TPauliMatrix &source);
   TPauliMatrix &operator+=(const TPauliMatrix &source);
   TPauliMatrix &operator+=(const Complex_t &factor);
   TPauliMatrix &operator+=(const LDouble_t &factor);
   TPauliMatrix &operator-=(const TPauliMatrix &source);
   TPauliMatrix &operator-=(const Complex_t &factor);
   TPauliMatrix &operator-=(const LDouble_t &factor);
   TPauliMatrix &operator*=(const TPauliMatrix &source);
   TPauliMatrix &operator*=(const Complex_t &factor);
   TPauliMatrix &operator*=(const LDouble_t &factor);
   TPauliMatrix &operator/=(const TPauliMatrix &source);
   TPauliMatrix &operator/=(const Complex_t &factor);
   TPauliMatrix &operator/=(const LDouble_t &factor);
 
   Bool_t operator==(const TPauliMatrix &other) const;
   Bool_t operator!=(const TPauliMatrix &other) const;
 
   TPauliMatrix &Zero();
   TPauliMatrix &Conj();
   TPauliMatrix &Invert();
   TPauliMatrix &Adjoint();
   TPauliMatrix &Transpose();
   TPauliMatrix &Compose
                (const LDouble_t a, const TThreeVectorReal &polar);
   TPauliMatrix &Compose
                (const Complex_t &a, const TThreeVectorComplex &polar);
   TPauliMatrix &SetDiagonal(const Complex_t &a);
   TPauliMatrix &SetDiagonal(const Complex_t &a11, const Complex_t &a22);
   TPauliMatrix &SetDensity(const TThreeVectorReal &polar);
   TPauliMatrix &SetRotation(const TThreeRotation &rotOp);
   TPauliMatrix &SetRotation(const TThreeVectorReal &axis);
   TPauliMatrix &SetRotation(const TUnitVector &axis, const LDouble_t angle);
   TPauliMatrix &SetRotation(const LDouble_t phi,
                             const LDouble_t theta,
                             const LDouble_t psi);
   TPauliMatrix &SimTransform(const TPauliMatrix &m);  // A' = M A Minverse
   TPauliMatrix &UniTransform(const TPauliMatrix &m);  // A' = M A Mdagger
 
   TPauliMatrix operator-() const;
   friend TPauliMatrix operator+(const TPauliMatrix &v1,
                                 const TPauliMatrix &v2);
   friend TPauliMatrix operator+(const Complex_t &a1,
                                 const TPauliMatrix &v2);
   friend TPauliMatrix operator+(const TPauliMatrix &v1,
                                 const Complex_t &a2);
   friend TPauliMatrix operator-(const TPauliMatrix &v1,
                                 const TPauliMatrix &v2);
   friend TPauliMatrix operator-(const Complex_t &a1,
                                 const TPauliMatrix &v2);
   friend TPauliMatrix operator-(const TPauliMatrix &v1,
                                 const Complex_t &a2);
   friend TPauliMatrix operator*(const TPauliMatrix &v1,
                                 const TPauliMatrix &v2);
   friend TPauliMatrix operator*(const Complex_t &a1,
                                 const TPauliMatrix &v2);
   friend TPauliMatrix operator*(const TPauliMatrix &v1,
                                 const Complex_t &a2);
   friend TPauliMatrix operator/(const TPauliMatrix &v1,
                                 const TPauliMatrix &v2);
   friend TPauliMatrix operator/(const Complex_t &a1,
                                 const TPauliMatrix &v2);
   friend TPauliMatrix operator/(const TPauliMatrix &v1,
                                 const Complex_t &a2);

   friend TBuffer &operator>>(TBuffer &buf, TPauliMatrix *&obj);
   friend TBuffer &operator<<(TBuffer &buf, const TPauliMatrix *obj);
   void Print(Option_t *option="");
 
   ClassDef(TPauliMatrix,1)  // Complex Pauli Matrix class
};


//----- inlines ----------------------------------------------------------------

inline TPauliMatrix::TPauliMatrix(const Int_t a)
{
   fMatrix[0][0] = a;	fMatrix[0][1] = 0;
   fMatrix[1][0] = 0;	fMatrix[1][1] = a;
}

inline TPauliMatrix::TPauliMatrix(const LDouble_t a)
{
   fMatrix[0][0] = a;	fMatrix[0][1] = 0;
   fMatrix[1][0] = 0;	fMatrix[1][1] = a;
}

inline TPauliMatrix::TPauliMatrix(const Complex_t &a)
{
   fMatrix[0][0] = a;	fMatrix[0][1] = 0;
   fMatrix[1][0] = 0;	fMatrix[1][1] = a;
}

inline Complex_t *TPauliMatrix::operator[](Int_t row) const
{
   return ((Complex_t *)fMatrix + row*2);
}

inline TPauliMatrix::TPauliMatrix
        (const Complex_t &a, const TThreeVectorComplex &b)
{
   Compose(a,b);
}

inline TPauliMatrix::TPauliMatrix(const TPauliMatrix &another)
{
   *this = another;
}

inline void TPauliMatrix::SetResolution(const LDouble_t resolution)
{
   fResolution = resolution;
}

inline LDouble_t TPauliMatrix::Resolution() const
{
   LDouble_t scale = abs(fMatrix[0][0]) + abs(fMatrix[0][1]) +
                    abs(fMatrix[1][0]) + abs(fMatrix[1][1]);
   if (scale > 0)
      return fResolution*scale;
   else
      return fResolution;
}

inline Bool_t TPauliMatrix::IsIdentity() const
{
   const TPauliMatrix one(1.0L);
   return (*this == one);
}

inline Bool_t TPauliMatrix::IsUnitary() const
{
   TPauliMatrix copy1(*this), copy2(*this);
   copy2 *= copy1.Adjoint();
   return copy2.IsIdentity();
}

inline Bool_t TPauliMatrix::IsHermetian() const
{
   TPauliMatrix adj(*this);
   return (*this == adj.Adjoint());
}

inline Bool_t TPauliMatrix::IsIdempotent() const
{
   TPauliMatrix thisSqr(*this);
   thisSqr *= *this;
   return (*this == thisSqr);
}

inline Complex_t TPauliMatrix::Trace() const
{
   return (fMatrix[0][0] + fMatrix[1][1]);
}

inline Complex_t TPauliMatrix::Determ() const
{
   return (fMatrix[0][0]*fMatrix[1][1] - fMatrix[1][0]*fMatrix[0][1]);
}

inline void TPauliMatrix::GetDiagonal(Complex_t &a11, Complex_t &a22) const
{
   a11 = fMatrix[0][0];
   a22 = fMatrix[1][1];
}

inline TPauliMatrix &TPauliMatrix::operator=(const TPauliMatrix &source)
{
   fMatrix[0][0] = source.fMatrix[0][0];
   fMatrix[0][1] = source.fMatrix[0][1];
   fMatrix[1][0] = source.fMatrix[1][0];
   fMatrix[1][1] = source.fMatrix[1][1];
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator+=(const TPauliMatrix &source)
{
   fMatrix[0][0] += source.fMatrix[0][0];
   fMatrix[0][1] += source.fMatrix[0][1];
   fMatrix[1][0] += source.fMatrix[1][0];
   fMatrix[1][1] += source.fMatrix[1][1];
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator+=(const Complex_t &factor)
{
   fMatrix[0][0] += factor;
   fMatrix[1][1] += factor;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator+=(const LDouble_t &factor)
{
   fMatrix[0][0] += factor;
   fMatrix[1][1] += factor;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator-=(const TPauliMatrix &source)
{
   fMatrix[0][0] -= source.fMatrix[0][0];
   fMatrix[0][1] -= source.fMatrix[0][1];
   fMatrix[1][0] -= source.fMatrix[1][0];
   fMatrix[1][1] -= source.fMatrix[1][1];
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator-=(const Complex_t &factor)
{
   fMatrix[0][0] -= factor;
   fMatrix[1][1] -= factor;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator-=(const LDouble_t &factor)
{
   fMatrix[0][0] -= factor;
   fMatrix[1][1] -= factor;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator*=(const Complex_t &factor)
{
   fMatrix[0][0] *= factor;
   fMatrix[0][1] *= factor;
   fMatrix[1][0] *= factor;
   fMatrix[1][1] *= factor;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator*=(const LDouble_t &factor)
{
   fMatrix[0][0] *= factor;
   fMatrix[0][1] *= factor;
   fMatrix[1][0] *= factor;
   fMatrix[1][1] *= factor;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator/=(const TPauliMatrix &source)
{
   TPauliMatrix sourceCopy(source);
   return (*this *= sourceCopy.Invert());
}

inline TPauliMatrix &TPauliMatrix::operator/=(const Complex_t &factor)
{
   fMatrix[0][0] /= factor;
   fMatrix[0][1] /= factor;
   fMatrix[1][0] /= factor;
   fMatrix[1][1] /= factor;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::operator/=(const LDouble_t &factor)
{
   fMatrix[0][0] /= factor;
   fMatrix[0][1] /= factor;
   fMatrix[1][0] /= factor;
   fMatrix[1][1] /= factor;
   return *this;
}
 
inline Bool_t TPauliMatrix::operator!=(const TPauliMatrix &other) const
{
   if (*this == other)
      return 0;
   else
      return 1;
}
 
inline TPauliMatrix &TPauliMatrix::Zero()
{
   fMatrix[0][0] = fMatrix[0][1] = fMatrix[1][0] = fMatrix[1][1] = 0;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::Conj()
{
   fMatrix[0][0] = conj(fMatrix[0][0]);
   fMatrix[0][1] = conj(fMatrix[0][1]);
   fMatrix[1][0] = conj(fMatrix[1][0]);
   fMatrix[1][1] = conj(fMatrix[1][1]);
   return *this;
}

inline TPauliMatrix &TPauliMatrix::Adjoint()
{
   Transpose(); Conj();
   return *this;
}

inline TPauliMatrix &TPauliMatrix::Transpose()
{
   Complex_t temp = fMatrix[0][1];
   fMatrix[0][1] = fMatrix[1][0];
   fMatrix[1][0] = temp;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::SetDiagonal(const Complex_t &a)
{
   fMatrix[0][0] = fMatrix[1][1] = a;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::SetDiagonal
                    (const Complex_t &a11, const Complex_t &a22)
{
   fMatrix[0][0] = a11; fMatrix[1][1] = a22;
   return *this;
}

inline TPauliMatrix &TPauliMatrix::SetDensity(const TThreeVectorReal &polar)
{
   Compose(1,polar);
   return (*this /= Complex_t(2));
}

inline TPauliMatrix &TPauliMatrix::SetRotation(const TThreeVectorReal &axis)
{
   LDouble_t angle = axis.Length();
   return SetRotation(axis,angle);
}

inline TPauliMatrix TPauliMatrix::operator-() const
{
   TPauliMatrix result;
   result.fMatrix[0][0] = -fMatrix[0][0];
   result.fMatrix[0][1] = -fMatrix[0][1];
   result.fMatrix[1][0] = -fMatrix[1][0];
   result.fMatrix[1][1] = -fMatrix[1][1];
   return result;
}

inline TPauliMatrix operator+(const TPauliMatrix &v1, const TPauliMatrix &v2)
{
   TPauliMatrix result(v1);
   return (result += v2);
}

inline TPauliMatrix operator+(const Complex_t &a1, const TPauliMatrix &v2)
{
   TPauliMatrix result(v2);
   return (result += a1);
}

inline TPauliMatrix operator+(const TPauliMatrix &v1, const Complex_t &a2)
{
   TPauliMatrix result(v1);
   return (result += a2);
}

inline TPauliMatrix operator-(const TPauliMatrix &v1, const TPauliMatrix &v2)
{
   TPauliMatrix result(v1);
   return (result -= v2);
}

inline TPauliMatrix operator-(const Complex_t &a1, const TPauliMatrix &v2)
{
   TPauliMatrix result(a1);
   return (result -= v2);
}

inline TPauliMatrix operator-(const TPauliMatrix &v1, const Complex_t &a2)
{
   TPauliMatrix result(v1);
   return (result -= a2);
}

inline TPauliMatrix operator*(const TPauliMatrix &v1, const TPauliMatrix &v2)
{
   TPauliMatrix result(v1);
   return (result *= v2);
}

inline TPauliMatrix operator*(const Complex_t &a1, const TPauliMatrix &v2)
{
   TPauliMatrix result(v2);
   return (result *= a1);
}

inline TPauliMatrix operator*(const TPauliMatrix &v1, const Complex_t &a2)
{
   TPauliMatrix result(v1);
   return (result *= a2);
}

inline TPauliMatrix operator/(const TPauliMatrix &v1, const TPauliMatrix &v2)
{
   TPauliMatrix result(v1);
   return (result /= v2);
}

inline TPauliMatrix operator/(const Complex_t &a1, const TPauliMatrix &v2)
{
   TPauliMatrix result(v2);
   result.Invert();
   return (result *= a1);
}

inline TPauliMatrix operator/(const TPauliMatrix &v1, const Complex_t &a2)
{
   TPauliMatrix result(v1);
   return (result /= a2);
}

inline TBuffer &operator>>(TBuffer &buf, TPauliMatrix *&obj)
{
   for (Int_t i=0; i<2; i++) {
      for (Int_t j=0; j<2; j++) {
         Double_t real,imag;
         buf >> real >> imag;
         obj->fMatrix[i][j] = Complex_t(real,imag);
      }
   }
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TPauliMatrix *obj)
{
   for (Int_t i=0; i<2; i++) {
      for (Int_t j=0; j<2; j++) {
         Double_t real = obj->fMatrix[i][j].real();
         Double_t imag = obj->fMatrix[i][j].imag();
         buf << real << imag;
      }
   }
   return buf;
}

#endif
