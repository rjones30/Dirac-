//
// TDiracMatrix.h 
//
// This file is distributed as part of the Dirac++ package,
// a general toolkit for computing the amplitudes for Feynman
// graphs. See DiracPackage.h for details.

#ifndef ROOT_TDiracMatrix
#define ROOT_TDiracMatrix
 
#include "Double.h"
#include "TFourVectorComplex.h"
#include "TPauliMatrix.h"
#include "TDiracSpinor.h"
#include "TError.h"
 
#include <math.h>

// special constants to indicate standard Dirac matrix forms

enum EDiracIndex
{
   kDiracOne=0,
   kDiracGamma0=4,
   kDiracGamma1=1,
   kDiracGamma2=2,
   kDiracGamma3=3,
   kDiracGamma4=4,
   kDiracGamma5=5,
   kDiracSigma1=23,
   kDiracSigma2=31,
   kDiracSigma3=12,
   kDiracKappa1=10,
   kDiracKappa2=20,
   kDiracKappa3=30
};

class TDeterminor;


class TDiracMatrix {

friend class TDiracSpinor;
 
protected:
   Complex_t        fMatrix[4][4];   // complex matrix allocated on stack
   static LDouble_t fResolution;     // matrix resolving "distance"
 
public:
   TDiracMatrix() { }
   explicit TDiracMatrix(const EDiracIndex i);
   explicit TDiracMatrix(const EDiracIndex i, const EDiracIndex j);
   explicit TDiracMatrix(const Int_t a);
   explicit TDiracMatrix(const LDouble_t a);
   explicit TDiracMatrix(const Complex_t &a);
   explicit TDiracMatrix(const TLorentzTransform &xOp);
   TDiracMatrix(const TDiracMatrix &another);
 
   virtual ~TDiracMatrix() { }

   Complex_t *operator[](Int_t row) const;

   static void SetResolution(const LDouble_t resolution);
   LDouble_t Resolution() const;
 
   Bool_t IsIdentity() const;
   Bool_t IsUnitary() const;
   Bool_t IsDiagonal() const;
   Bool_t IsAntiUnitary() const;
   Bool_t IsHermetian() const;
   Bool_t IsIdempotent() const;
   Complex_t Trace() const;
   Complex_t Determ() const;
   Complex_t Component(const EDiracIndex i) const;
   Complex_t Component(const EDiracIndex i, const EDiracIndex j) const;
   void GetDiagonal(Complex_t &a11, Complex_t &a22,
                    Complex_t &a33, Complex_t &a44);

   TDiracMatrix &operator=(const TDiracMatrix &source);
   TDiracMatrix &operator+=(const TDiracMatrix &source);
   TDiracMatrix &operator+=(const Complex_t &factor);
   TDiracMatrix &operator+=(const LDouble_t &factor);
   TDiracMatrix &operator-=(const TDiracMatrix &source);
   TDiracMatrix &operator-=(const Complex_t &factor);
   TDiracMatrix &operator-=(const LDouble_t &factor);
   TDiracMatrix &operator*=(const TDiracMatrix &source);
   TDiracMatrix &operator*=(const Complex_t &factor);
   TDiracMatrix &operator*=(const LDouble_t &factor);
   TDiracMatrix &operator/=(const TDiracMatrix &source);
   TDiracMatrix &operator/=(const Complex_t &factor);
   TDiracMatrix &operator/=(const LDouble_t &factor);
 
   Bool_t operator==(const TDiracMatrix &other) const;
   Bool_t operator!=(const TDiracMatrix &other) const;
 
   TDiracMatrix &Zero();
   TDiracMatrix &Conj();
   TDiracMatrix &Invert();
   TDiracMatrix &Adjoint();
   TDiracMatrix &Transpose();
   TDiracMatrix &SetDiagonal(const Complex_t &a);
   TDiracMatrix &SetDiagonal(const Complex_t &a11, const Complex_t &a22,
                             const Complex_t &a33, const Complex_t &a44);

   TDiracMatrix &SetUUbar(const TFourVectorReal &p);
   TDiracMatrix &SetVVbar(const TFourVectorReal &p);
   TDiracMatrix &SetUUbar(const TFourVectorReal &p, const Float_t helicity);
   TDiracMatrix &SetVVbar(const TFourVectorReal &p, const Float_t helicity);
   TDiracMatrix &SetUUbar(const TFourVectorReal &p, const TDiracSpinor &u);
   TDiracMatrix &SetVVbar(const TFourVectorReal &p, const TDiracSpinor &u);
   TDiracMatrix &SetUUbar(const TFourVectorReal &p,
                          const TThreeVectorReal &polar);
   TDiracMatrix &SetVVbar(const TFourVectorReal &p,
                          const TThreeVectorReal &polar);
   TDiracMatrix &SetUUbar(const TFourVectorReal &p,
                          const TPauliMatrix &density);
   TDiracMatrix &SetVVbar(const TFourVectorReal &p,
                          const TPauliMatrix &density);

   TDiracMatrix &SetRotation(const TThreeRotation &rotOp);
   TDiracMatrix &SetRotation(const TThreeVectorReal &axis);
   TDiracMatrix &SetRotation(const TUnitVector &axis, const LDouble_t angle);
   TDiracMatrix &SetRotation(const LDouble_t phi,
                             const LDouble_t theta,
                             const LDouble_t psi);
   TDiracMatrix &SetBoost(const TLorentzBoost &boostOp);
   TDiracMatrix &SetBoost(const LDouble_t betaX,
                          const LDouble_t betaY,
                          const LDouble_t betaZ);
   TDiracMatrix &SetBoost(const TThreeVectorReal &beta);
   TDiracMatrix &SetBoost(const TUnitVector &bhat, const LDouble_t &beta);
   TDiracMatrix &SetTransform(const TLorentzTransform &xform);
   TDiracMatrix &SimTransform(const TDiracMatrix &M); // A' = M A Minverse
   TDiracMatrix &UniTransform(const TDiracMatrix &U); // A' = U A Udagger
   TDiracMatrix &Slash(const TFourVectorReal &p);
   TDiracMatrix &Slash(const TFourVectorComplex &A);
 
   TDiracMatrix &operator-();
   friend TDiracMatrix operator+(const TDiracMatrix &v1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator+(const Complex_t &a1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator+(const TDiracMatrix &v1,
                                 const Complex_t &a2);
   friend TDiracMatrix operator+(const LDouble_t &a1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator+(const TDiracMatrix &v1,
                                 const LDouble_t &a2);
   friend TDiracMatrix operator-(const TDiracMatrix &v1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator-(const Complex_t &a1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator-(const TDiracMatrix &v1,
                                 const Complex_t &a2);
   friend TDiracMatrix operator-(const LDouble_t &a1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator-(const TDiracMatrix &v1,
                                 const LDouble_t &a2);
   friend TDiracMatrix operator*(const TDiracMatrix &v1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator*(const Complex_t &a1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator*(const TDiracMatrix &v1,
                                 const Complex_t &a2);
   friend TDiracMatrix operator*(const LDouble_t &a1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator*(const TDiracMatrix &v1,
                                 const LDouble_t &a2);
   friend TDiracMatrix operator/(const TDiracMatrix &v1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator/(const Complex_t &a1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator/(const TDiracMatrix &v1,
                                 const Complex_t &a2);
   friend TDiracMatrix operator/(const LDouble_t &a1,
                                 const TDiracMatrix &v2);
   friend TDiracMatrix operator/(const TDiracMatrix &v1,
                                 const LDouble_t &a2);

   friend TBuffer &operator>>(TBuffer &buf, TDiracMatrix *&obj);
   friend TBuffer &operator<<(TBuffer &buf, const TDiracMatrix *obj);
   void Print(Option_t *option="");
 
   ClassDef(TDiracMatrix,1)  // Complex Dirac Matrix class
};

//----- inlines ----------------------------------------------------------------

inline TDiracMatrix::TDiracMatrix(const Int_t a)
{
   Zero();
   SetDiagonal(a);
}

inline TDiracMatrix::TDiracMatrix(const LDouble_t a)
{
   Zero();
   SetDiagonal(a);
}

inline TDiracMatrix::TDiracMatrix(const Complex_t &a)
{
   Zero();
   SetDiagonal(a);
}

inline TDiracMatrix::TDiracMatrix(const TDiracMatrix &another)
{
   *this = another;
}

inline TDiracMatrix::TDiracMatrix(const TLorentzTransform &xOp)
{
   SetTransform(xOp);
}

inline Complex_t *TDiracMatrix::operator[](Int_t row) const
{
   return ((Complex_t *)fMatrix + row*4);
}
 
inline void TDiracMatrix::SetResolution(const LDouble_t resolution)
{
   fResolution = resolution;
}

inline LDouble_t TDiracMatrix::Resolution() const
{
   LDouble_t scale = abs(fMatrix[0][0]) + abs(fMatrix[0][1]) +
                    abs(fMatrix[0][2]) + abs(fMatrix[0][3]);
   if (scale > 0)
      return fResolution*scale;
   else
      return fResolution;
}

inline Bool_t TDiracMatrix::IsIdentity() const
{
   const TDiracMatrix one(1.0L);
   return (*this == one);
}

inline Bool_t TDiracMatrix::IsUnitary() const
{
   TDiracMatrix copy1(*this), copy2(*this);
   copy2 *= copy1.Adjoint();
   return copy2.IsIdentity();
}

inline Bool_t TDiracMatrix::IsDiagonal() const
{
   LDouble_t limit = Resolution();
   if ( abs(fMatrix[0][1]) >= limit ||
        abs(fMatrix[0][2]) >= limit ||
        abs(fMatrix[0][3]) >= limit ||
        abs(fMatrix[1][2]) >= limit ||
        abs(fMatrix[1][3]) >= limit ||
        abs(fMatrix[2][3]) >= limit ||
        abs(fMatrix[1][0]) >= limit ||
        abs(fMatrix[2][0]) >= limit ||
        abs(fMatrix[3][0]) >= limit ||
        abs(fMatrix[2][1]) >= limit ||
        abs(fMatrix[3][1]) >= limit || 
        abs(fMatrix[3][2]) >= limit )
      return 0;
   else
      return 1;
}

inline Bool_t TDiracMatrix::IsAntiUnitary() const
{
   TDiracMatrix copy1(*this), copy2(*this);
   copy2 *= -copy1.Adjoint();
   return copy2.IsIdentity();
}

inline Bool_t TDiracMatrix::IsHermetian() const
{
   TDiracMatrix adj(*this);
   return (*this == adj.Adjoint());
}

inline Bool_t TDiracMatrix::IsIdempotent() const
{
   TDiracMatrix thisSqr(*this);
   thisSqr *= *this;
   return (*this == thisSqr);
}

inline Complex_t TDiracMatrix::Trace() const
{
   return (fMatrix[0][0] + fMatrix[1][1] + fMatrix[2][2] + fMatrix[3][3]);
}

inline Complex_t TDiracMatrix::Component(const EDiracIndex i) const
{
// Returns coefficient b_i from the general expansion for Dirac matrix M
//  M = a + sum{i=1,5}(b_i * gamma_i) + sum{i=1,5;j=i+1,5}(c_ij * gamma_ij)
// but where index 0 and 4 are interchangeable.  For the definition of
// the gamma matrices gamma_i and commutators gamma_ij, see constructor.

   TDiracMatrix g(i);
   g *= *this;
   return (g.Trace() / ((i>0 && i<4) ? -4.L : 4.L));
}

inline Complex_t TDiracMatrix::Component(const EDiracIndex i,
                                         const EDiracIndex j) const
{
// Returns coefficient c_ij from the general expansion for Dirac matrix M
//  M = a + sum{i=1,5}(b_i * gamma_i) + sum{i=1,5;j=i+1,5}(c_ij * gamma_ij)
// but where index 0 and 4 are interchangeable.  For the definition of
// the gamma matrices gamma_i and commutators gamma_ij, see constructor.

   TDiracMatrix g(i,j);
   g *= *this;
   return (g.Trace() / ((i>0 && i<4) ? -2.L : 2.L) / ((j>0 && j<4) ? -2.L : 2.L));
}

inline void TDiracMatrix::GetDiagonal(Complex_t &a11, Complex_t &a22,
                                      Complex_t &a33, Complex_t &a44)
{
   a11 = fMatrix[0][0];
   a22 = fMatrix[1][1];
   a33 = fMatrix[2][2];
   a44 = fMatrix[3][3];
}

inline TDiracMatrix &TDiracMatrix::operator=(const TDiracMatrix &source)
{
   Complex_t *s = (Complex_t *)&source.fMatrix[0][0];
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++, s++)
      *d = *s;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator+=(const TDiracMatrix &source)
{
   Complex_t *s = (Complex_t *)&source.fMatrix[0][0];
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++, s++)
      *d += *s;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator+=(const Complex_t &factor)
{
   fMatrix[0][0] += factor;
   fMatrix[1][1] += factor;
   fMatrix[2][2] += factor;
   fMatrix[3][3] += factor;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator+=(const LDouble_t &factor)
{
   fMatrix[0][0] += factor;
   fMatrix[1][1] += factor;
   fMatrix[2][2] += factor;
   fMatrix[3][3] += factor;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator-=(const TDiracMatrix &source)
{
   Complex_t *s = (Complex_t *)&source.fMatrix[0][0];
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++, s++)
      *d -= *s;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator-=(const Complex_t &factor)
{
   fMatrix[0][0] -= factor;
   fMatrix[1][1] -= factor;
   fMatrix[2][2] -= factor;
   fMatrix[3][3] -= factor;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator-=(const LDouble_t &factor)
{
   fMatrix[0][0] -= factor;
   fMatrix[1][1] -= factor;
   fMatrix[2][2] -= factor;
   fMatrix[3][3] -= factor;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator*=(const TDiracMatrix &source)
{
   TDiracMatrix copy(*this);
   for (Int_t i=0; i<4; i++)
      for (Int_t j=0; j<4; j++) {
         Complex_t sum=0;
         for (Int_t k=0; k<4; k++)
            sum += copy.fMatrix[i][k]*source.fMatrix[k][j];
         fMatrix[i][j] = sum;
      }
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator*=(const Complex_t &factor)
{
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++)
      *d *= factor;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator*=(const LDouble_t &factor)
{
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++)
      *d *= factor;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::operator/=(const TDiracMatrix &source)
{
   TDiracMatrix sourceCopy(source);
   return (*this *= sourceCopy.Invert());
}

inline TDiracMatrix &TDiracMatrix::operator/=(const Complex_t &factor)
{
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++)
      *d /= factor;
   return *this;
}
 
inline TDiracMatrix &TDiracMatrix::operator/=(const LDouble_t &factor)
{
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++)
      *d /= factor;
   return *this;
}
 
inline Bool_t TDiracMatrix::operator==(const TDiracMatrix &other) const
{
   Complex_t *me = (Complex_t *)&fMatrix[0][0];
   Complex_t *her = (Complex_t *)&other.fMatrix[0][0];
   for (Int_t i=0; i<16; i++, me++, her++)
      if (abs(*me - *her) >= Resolution())
         return 0;
   return 1;
}

inline Bool_t TDiracMatrix::operator!=(const TDiracMatrix &other) const
{
   if (*this == other)
      return 0;
   else
      return 1;
}
 
inline TDiracMatrix &TDiracMatrix::Zero()
{
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++)
      *d = 0;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::Conj()
{
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++)
      *d = conj(*d);
   return *this;
}

inline TDiracMatrix &TDiracMatrix::Adjoint()
{
   Transpose(); Conj();
   return *this;
}

inline TDiracMatrix &TDiracMatrix::SetDiagonal(const Complex_t &a)
{
   fMatrix[0][0] = fMatrix[1][1] = fMatrix[2][2] = fMatrix[3][3] = a;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::SetDiagonal
                    (const Complex_t &a11, const Complex_t &a22,
                     const Complex_t &a33, const Complex_t &a44)
{
   fMatrix[0][0] = a11;
   fMatrix[1][1] = a22;
   fMatrix[2][2] = a33;
   fMatrix[3][3] = a44;
   return *this;
}

inline TDiracMatrix &TDiracMatrix::SetUUbar(const TFourVectorReal &p)
{
   Slash(p);
   return (*this += p.Invariant());
}

inline TDiracMatrix &TDiracMatrix::SetVVbar(const TFourVectorReal &p)
{
   Slash(p);
   return (*this -= p.Invariant());
}

inline TDiracMatrix &TDiracMatrix::SetUUbar
                    (const TFourVectorReal &p, const Float_t helicity)
{
   TDiracSpinor u(p,helicity);
   return SetUUbar(p,u);
}

inline TDiracMatrix &TDiracMatrix::SetVVbar
                    (const TFourVectorReal &p, const Float_t helicity)
{
   TDiracSpinor v(-p,helicity);
   return SetVVbar(p,v);
}

inline TDiracMatrix &TDiracMatrix::SetUUbar
                    (const TFourVectorReal &p, const TPauliMatrix &density)
{
   LDouble_t a;
   TThreeVectorReal polar;
   density.Decompose(a,polar);
   SetUUbar(p,polar/=a);
   return *this *= 2*a;
}
  
inline TDiracMatrix &TDiracMatrix::SetVVbar
                    (const TFourVectorReal &p, const TPauliMatrix &density)
{
   LDouble_t a;
   TThreeVectorReal polar;
   density.Decompose(a,polar);
   SetVVbar(p,polar/=a);
   return *this *= 2*a;
}

inline TDiracMatrix &TDiracMatrix::Slash(const TFourVectorReal &p)
{
   return (*this = p[0]*TDiracMatrix(kDiracGamma0) -
                   p[1]*TDiracMatrix(kDiracGamma1) - 
                   p[2]*TDiracMatrix(kDiracGamma2) -
                   p[3]*TDiracMatrix(kDiracGamma3) );
}

inline TDiracMatrix &TDiracMatrix::Slash(const TFourVectorComplex &a)
{
   return (*this = a[0]*TDiracMatrix(kDiracGamma0) - 
                   a[1]*TDiracMatrix(kDiracGamma1) -
                   a[2]*TDiracMatrix(kDiracGamma2) -
                   a[3]*TDiracMatrix(kDiracGamma3) );
}

inline TDiracMatrix &TDiracMatrix::operator-()
{
   Complex_t *d = (Complex_t *)&fMatrix[0][0];
   for (Int_t i=0; i<16; i++, d++)
      *d = -(*d);
   return *this;
}

inline TDiracMatrix operator+(const TDiracMatrix &v1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(v1);
   return (result += v2);
}

inline TDiracMatrix operator+(const Complex_t &a1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(a1);
   return (result += v2);
}

inline TDiracMatrix operator+(const TDiracMatrix &v1,
                              const Complex_t &a2)
{
   TDiracMatrix result(v1);
   return result += a2;
}

inline TDiracMatrix operator+(const LDouble_t &a1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(a1);
   return (result += v2);
}

inline TDiracMatrix operator+(const TDiracMatrix &v1,
                              const LDouble_t &a2)
{
   TDiracMatrix result(v1);
   return result += a2;
}

inline TDiracMatrix operator-(const TDiracMatrix &v1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(v1);
   return result -= v2;
}

inline TDiracMatrix operator-(const Complex_t &a1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(a1);
   return result -= v2;
}

inline TDiracMatrix operator-(const TDiracMatrix &v1,
                              const Complex_t &a2)
{
   TDiracMatrix result(v1);
   return result -= a2;
}

inline TDiracMatrix operator-(const LDouble_t &a1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(a1);
   return result -= v2;
}

inline TDiracMatrix operator-(const TDiracMatrix &v1,
                              const LDouble_t &a2)
{
   TDiracMatrix result(v1);
   return result -= a2;
}

inline TDiracMatrix operator*(const TDiracMatrix &v1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(v1);
   return result *= v2;
}

inline TDiracMatrix operator*(const Complex_t &a1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(a1);
   return result *= v2;
}

inline TDiracMatrix operator*(const TDiracMatrix &v1,
                              const Complex_t &a2)
{
   TDiracMatrix result(v1);
   return result *= a2;
}

inline TDiracMatrix operator*(const LDouble_t &a1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(a1);
   return result *= v2;
}

inline TDiracMatrix operator*(const TDiracMatrix &v1,
                              const LDouble_t &a2)
{
   TDiracMatrix result(v1);
   return result *= a2;
}

inline TDiracMatrix operator/(const TDiracMatrix &v1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(v1);
   return result /= v2;
}

inline TDiracMatrix operator/(const Complex_t &a1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(a1);
   return result /= v2;
}

inline TDiracMatrix operator/(const TDiracMatrix &v1,
                              const Complex_t &a2)
{
   TDiracMatrix result(v1);
   return result /= a2;
}

inline TDiracMatrix operator/(const LDouble_t &a1,
                              const TDiracMatrix &v2)
{
   TDiracMatrix result(a1);
   return result /= v2;
}

inline TDiracMatrix operator/(const TDiracMatrix &v1,
                              const LDouble_t &a2)
{
   TDiracMatrix result(v1);
   return result /= a2;
}

inline TBuffer &operator>>(TBuffer &buf, TDiracMatrix *&obj)
{
   for (Int_t i=0; i<4; i++) {
      for (Int_t j=0; j<4; j++) {
         Double_t real,imag;
         buf >> real >> imag;
         obj->fMatrix[i][j] = Complex_t(real,imag);
      }
   }
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TDiracMatrix *obj)
{
   for (Int_t i=0; i<4; i++) {
      for (Int_t j=0; j<4; j++) {
         Double_t real = obj->fMatrix[i][j].real();
         Double_t imag = obj->fMatrix[i][j].imag();
         buf << real << imag;
      }
   }
   return buf;
}

#endif
