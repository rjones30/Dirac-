//
// TFourVectorComplex.h 
//
// This file is distributed as part of the Lorentz++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See LorentzPackage.h for details.

#ifndef ROOT_TFourVectorComplex
#define ROOT_TFourVectorComplex
 
#include "Double.h"
#include "TBuffer.h"
#include "TThreeVectorComplex.h" 
#include "TFourVectorReal.h" 
#include "TError.h"
 
#include <math.h>

class TLorentzTransform;
class TLorentzBoost;
class TThreeRotation;
 

class TFourVectorComplex : public TThreeVectorComplex {
 
friend class TLorentzTransform;
friend class TLorentzBoost;
friend class TThreeRotation;

public:
   TFourVectorComplex() : TThreeVectorComplex() { }
   explicit TFourVectorComplex(const Complex_t &t,
                               const Complex_t &x,
                               const Complex_t &y,
                               const Complex_t &z);
   explicit TFourVectorComplex(const Float_t *array);
   explicit TFourVectorComplex(const LDouble_t *array);
   explicit TFourVectorComplex(const Complex_t *array);
   explicit TFourVectorComplex(const Complex_t &t, const TThreeVectorComplex &another);
   TFourVectorComplex(const TFourVectorReal &another);
   TFourVectorComplex(const TFourVectorComplex &another);
 
   virtual ~TFourVectorComplex() { }

   LDouble_t Resolution() const;
 
   Complex_t &operator[](const Int_t index) const;
 
   LDouble_t Invariant() const;       // use +--- metric and
   LDouble_t InvariantSqr() const;    // sum over absolute-squares
   TFourVectorReal RealPart() const;
   TFourVectorReal ImagPart() const;
   void GetCoord(Complex_t &t, Complex_t &x, Complex_t &y, Complex_t &z) const;
   void GetCoord(Complex_t *array) const;
   LDouble_t DistanceTo(const Complex_t t, const Complex_t x,
                       const Complex_t y, const Complex_t z) const;
   LDouble_t DistanceTo(const Float_t *array) const;
   LDouble_t DistanceTo(const LDouble_t *array) const;
   LDouble_t DistanceTo(const Complex_t *array) const;
   LDouble_t DistanceTo(const TFourVectorComplex &vec2) const;
 
   TFourVectorComplex &operator=(const TFourVectorReal &source);
   TFourVectorComplex &operator=(const TFourVectorComplex &source);
   TFourVectorComplex &operator=(const Float_t *array);
   TFourVectorComplex &operator=(const LDouble_t *array);
   TFourVectorComplex &operator=(const Complex_t *array);
   TFourVectorComplex &operator+=(const TFourVectorReal &source);
   TFourVectorComplex &operator+=(const TFourVectorComplex &source);
   TFourVectorComplex &operator+=(const Float_t *array);
   TFourVectorComplex &operator+=(const LDouble_t *array);
   TFourVectorComplex &operator+=(const Complex_t *array);
   TFourVectorComplex &operator-=(const TFourVectorReal &source);
   TFourVectorComplex &operator-=(const TFourVectorComplex &source);
   TFourVectorComplex &operator-=(const Float_t *array);
   TFourVectorComplex &operator-=(const LDouble_t *array);
   TFourVectorComplex &operator-=(const Complex_t *array);
   TFourVectorComplex &operator*=(const Complex_t &factor);
   TFourVectorComplex &operator/=(const Complex_t &factor);
 
   Bool_t operator==(const TFourVectorComplex &other) const;
   Bool_t operator!=(const TFourVectorComplex &other) const;
 
   TFourVectorComplex &Zero();
   TFourVectorComplex &Conj();
   TFourVectorComplex &Transform(const TLorentzTransform &xformOp);
   TFourVectorComplex &Boost(const TLorentzBoost &boostOp);
   TFourVectorComplex &Boost(const LDouble_t betaX,
                             const LDouble_t betaY,
                             const LDouble_t betaZ);
   TFourVectorComplex &Boost(const LDouble_t *beta);
   TFourVectorComplex &Boost(const TThreeVectorReal &beta);
   TFourVectorComplex &Boost(const TUnitVector &bhat, const LDouble_t beta);
   TFourVectorComplex &BoostToRest(const TFourVector &p);
   TFourVectorComplex &BoostFromRest(const TFourVector &p);
   Complex_t ScalarProd(const TFourVectorComplex &other);
   Complex_t ScalarProd(const TFourVectorComplex &v1,
                        const TFourVectorComplex &v2);
 
   TFourVectorComplex operator-() const;
   friend TFourVectorComplex operator+(const TFourVectorComplex &v1,
                                       const TFourVectorComplex &v2);
   friend TFourVectorComplex operator-(const TFourVectorComplex &v1,
                                       const TFourVectorComplex &v2);
   friend TFourVectorComplex operator*(const TFourVectorComplex &vec,
                                       const Complex_t &factor);
   friend TFourVectorComplex operator*(const Complex_t &factor,
                                       const TFourVectorComplex &vec);
   friend TFourVectorComplex operator/(const TFourVectorComplex &vec,
                                       const Complex_t &factor);

   friend TBuffer &operator>>(TBuffer &buf, TFourVectorComplex *&vec);
   friend TBuffer &operator<<(TBuffer &buf, const TFourVectorComplex *vec);
   void Print(Option_t *option="");
 
   ClassDef(TFourVectorComplex,1)  // Complex four-vector class
};

//----- inlines ----------------------------------------------------------------
 
inline TFourVectorComplex::TFourVectorComplex
       (const Complex_t &t, const Complex_t &x,
        const Complex_t &y, const Complex_t &z)
 : TThreeVectorComplex()
{
   fVector[0] = t;
   fVector[1] = x;
   fVector[2] = y;
   fVector[3] = z;
}

inline TFourVectorComplex::TFourVectorComplex(const Float_t *array)
 : TThreeVectorComplex()
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TFourVectorComplex::TFourVectorComplex(const LDouble_t *array)
 : TThreeVectorComplex()
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TFourVectorComplex::TFourVectorComplex(const Complex_t *array)
 : TThreeVectorComplex()
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TFourVectorComplex::TFourVectorComplex(const TFourVectorComplex &another)
 : TThreeVectorComplex()
{
   *this = another;
}

inline TFourVectorComplex::TFourVectorComplex(const TFourVectorReal &another)
 : TThreeVectorComplex()
{
   *this = another;
}
 
 
inline TFourVectorComplex::TFourVectorComplex
       (const Complex_t &t, const TThreeVectorComplex &another)
 : TThreeVectorComplex(another)
{
   fVector[0] = t;
}
 
inline LDouble_t TFourVectorComplex::Resolution() const
{
   LDouble_t scale = sqrt(norm(fVector[0]) + LengthSqr());
   if (scale > 0)
      return fResolution*scale;
   else
      return fResolution;
}
 
inline Complex_t &TFourVectorComplex::operator[](const Int_t index) const
{
   if (index < 0 || index > 3) {
      Error("TFourVectorComplex::operator[]","index out of range");
      return (Complex_t &)fVector[0];
   }
   return (Complex_t &)fVector[index];
}

inline LDouble_t TFourVectorComplex::Invariant() const
{
   LDouble_t inv2 = InvariantSqr();
   if (inv2 > 0)
      return sqrt(inv2);
   else if (inv2 > -Resolution())
      return 0;
   else
      Error("TFourVectorComplex::Invariant()",
            "invoked on a vector with negative norm");
   return -1;
}

inline LDouble_t TFourVectorComplex::InvariantSqr() const
{
   return (norm(fVector[0]) - LengthSqr());
}

inline TFourVectorReal TFourVectorComplex::RealPart() const
{
   TFourVectorReal R;
   R.fVector[0] = real(fVector[0]);
   R.fVector[1] = real(fVector[1]);
   R.fVector[2] = real(fVector[2]);
   R.fVector[3] = real(fVector[3]);
   return R;
}

inline TFourVectorReal TFourVectorComplex::ImagPart() const
{
   TFourVectorReal I;
   I.fVector[0] = imag(fVector[0]);
   I.fVector[1] = imag(fVector[1]);
   I.fVector[2] = imag(fVector[2]);
   I.fVector[3] = imag(fVector[3]);
   return I;
}

inline void TFourVectorComplex::GetCoord
            (Complex_t &t, Complex_t &x, Complex_t &y, Complex_t &z) const
{
   t = fVector[0];
   x = fVector[1];
   y = fVector[2];
   z = fVector[3];
}

inline void TFourVectorComplex::GetCoord(Complex_t *array) const
{
   *(array++) = fVector[0];
   *(array++) = fVector[1];
   *(array++) = fVector[2];
   *array     = fVector[3];
}
 
inline LDouble_t TFourVectorComplex::DistanceTo
                (const Complex_t t, const Complex_t x,
                 const Complex_t y, const Complex_t z) const
{
   Complex_t tloc(t), xloc(x), yloc(y), zloc(z);
   tloc -= fVector[0];
   xloc -= fVector[1];
   yloc -= fVector[2];
   zloc -= fVector[3];
   return sqrt(norm(tloc) + norm(xloc) + norm(yloc) + norm(zloc));
}

inline LDouble_t TFourVectorComplex::DistanceTo(const Float_t *array) const
{
   Complex_t t = fVector[0] - (LDouble_t)*(array++);
   Complex_t x = fVector[1] - (LDouble_t)*(array++);
   Complex_t y = fVector[2] - (LDouble_t)*(array++);
   Complex_t z = fVector[3] - (LDouble_t)*array;
   return sqrt(norm(t) + norm(x) + norm(y) + norm(z));
}

inline LDouble_t TFourVectorComplex::DistanceTo(const LDouble_t *array) const
{
   Complex_t t = fVector[0] - *(array++);
   Complex_t x = fVector[1] - *(array++);
   Complex_t y = fVector[2] - *(array++);
   Complex_t z = fVector[3] - *array;
   return sqrt(norm(t) + norm(x) + norm(y) + norm(z));
}

inline LDouble_t TFourVectorComplex::DistanceTo(const Complex_t *array) const
{
   Complex_t t = fVector[0] - *(array++);
   Complex_t x = fVector[1] - *(array++);
   Complex_t y = fVector[2] - *(array++);
   Complex_t z = fVector[3] - *array;
   return sqrt(norm(t) + norm(x) + norm(y) + norm(z));
}

inline LDouble_t TFourVectorComplex::DistanceTo
                (const TFourVectorComplex &vec2) const
{
   Complex_t t = fVector[0] - vec2.fVector[0];
   Complex_t x = fVector[1] - vec2.fVector[1];
   Complex_t y = fVector[2] - vec2.fVector[2];
   Complex_t z = fVector[3] - vec2.fVector[3];
   return sqrt(norm(t) + norm(x) + norm(y) + norm(z));
}

inline TFourVectorComplex &TFourVectorComplex::operator=
                       (const TFourVectorReal &source)
{
   fVector[0] = source.fVector[0];
   fVector[1] = source.fVector[1];
   fVector[2] = source.fVector[2];
   fVector[3] = source.fVector[3];
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator=
                       (const TFourVectorComplex &source)
{
   fVector[0] = source.fVector[0];
   fVector[1] = source.fVector[1];
   fVector[2] = source.fVector[2];
   fVector[3] = source.fVector[3];
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator=(const Float_t *array)
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator=(const LDouble_t *array)
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator=(const Complex_t *array)
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator+=
                        (const TFourVectorComplex &source)
{
   fVector[0] += source.fVector[0];
   fVector[1] += source.fVector[1];
   fVector[2] += source.fVector[2];
   fVector[3] += source.fVector[3];
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator+=
                        (const TFourVectorReal &source)
{
   fVector[0] += source.fVector[0];
   fVector[1] += source.fVector[1];
   fVector[2] += source.fVector[2];
   fVector[3] += source.fVector[3];
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator+=(const Float_t *array)
{
   fVector[0] += *(array++);
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator+=(const LDouble_t *array)
{
   fVector[0] += *(array++);
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator+=(const Complex_t *array)
{
   fVector[0] += *(array++);
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator-=
                        (const TFourVectorComplex &source)
{
   fVector[0] -= source.fVector[0];
   fVector[1] -= source.fVector[1];
   fVector[2] -= source.fVector[2];
   fVector[3] -= source.fVector[3];
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator-=
                        (const TFourVectorReal &source)
{
   fVector[0] -= source.fVector[0];
   fVector[1] -= source.fVector[1];
   fVector[2] -= source.fVector[2];
   fVector[3] -= source.fVector[3];
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator-=
                          (const Float_t *array)
{
   fVector[0] -= *(array++);
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator-=
                          (const LDouble_t *array)
{
   fVector[0] -= *(array++);
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator-=
                          (const Complex_t *array)
{
   fVector[0] -= *(array++);
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator*=
                          (const Complex_t &factor)
{
   fVector[0] *= factor;
   fVector[1] *= factor;
   fVector[2] *= factor;
   fVector[3] *= factor;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::operator/=
                          (const Complex_t &factor)
{
   fVector[0] /= factor;
   fVector[1] /= factor;
   fVector[2] /= factor;
   fVector[3] /= factor;
   return *this;
}
 
inline Bool_t TFourVectorComplex::operator==
              (const TFourVectorComplex &other) const
{
   return (DistanceTo(other) < Resolution());
}

inline Bool_t TFourVectorComplex::operator!=
              (const TFourVectorComplex &other) const
{
   return !(*this == other);
}

inline TFourVectorComplex &TFourVectorComplex::Zero()
{
   fVector[0] = 0;
   fVector[1] = 0;
   fVector[2] = 0;
   fVector[3] = 0;
   return *this;
}

inline TFourVectorComplex &TFourVectorComplex::Conj()
{
   fVector[0] = conj(fVector[0]);
   fVector[1] = conj(fVector[1]);
   fVector[2] = conj(fVector[2]);
   fVector[3] = conj(fVector[3]);
   return *this;
}

inline Complex_t TFourVectorComplex::ScalarProd
                 (const TFourVectorComplex &other)
{
   return Complex_t(fVector[0] * other.fVector[0] -
                    fVector[1] * other.fVector[1] -
                    fVector[2] * other.fVector[2] -
                    fVector[3] * other.fVector[3]);
}

inline Complex_t TFourVectorComplex::ScalarProd
                 (const TFourVectorComplex &v1, const TFourVectorComplex &v2)
{
   return Complex_t(v1.fVector[0] * v2.fVector[0] -
                    v1.fVector[1] * v2.fVector[1] -
                    v1.fVector[2] * v2.fVector[2] -
                    v1.fVector[3] * v2.fVector[3]);
}

inline TFourVectorComplex TFourVectorComplex::operator-() const
{
   TFourVectorComplex minusThis;
   minusThis.fVector[0] = -fVector[0];
   minusThis.fVector[1] = -fVector[1];
   minusThis.fVector[2] = -fVector[2];
   minusThis.fVector[3] = -fVector[3];
   return minusThis;
}

inline TFourVectorComplex operator+
       (const TFourVectorComplex &v1, const TFourVectorComplex &v2)
{
   TFourVectorComplex result(v1);
   return (result += v2);
}

inline TFourVectorComplex operator-
       (const TFourVectorComplex &v1, const TFourVectorComplex &v2)
{
   TFourVectorComplex result(v1);
   return (result -= v2);
}

inline TFourVectorComplex operator*
       (const TFourVectorComplex &vec, const Complex_t &factor)
{
   TFourVectorComplex result(vec);
   return (result *= factor);
}

inline TFourVectorComplex operator*
       (const Complex_t &factor, const TFourVectorComplex &vec)
{
   TFourVectorComplex result(vec);
   return (result *= factor);
}

inline TFourVectorComplex operator/
       (const TFourVectorComplex &vec, const Complex_t &factor)
{
   TFourVectorComplex result(vec);
   return (result /= factor);
}

inline TBuffer &operator>>(TBuffer &buf, TFourVectorComplex *&obj)
{
   for (Int_t i=0; i<4; i++) {
      Double_t real,imag;
      buf >> real >> imag;
      obj->fVector[i] = Complex_t(real,imag);
   }
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TFourVectorComplex *obj)
{
// This method assumes that complex is stored in memory as LDouble_t[2]
   for (Int_t i=0; i<4; i++) {
      Double_t real = obj->fVector[i].real();
      Double_t imag = obj->fVector[i].imag();
      buf << real << imag;
   }
   return buf;
}

#endif
