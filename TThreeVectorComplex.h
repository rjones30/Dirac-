//
// TThreeVectorComplex.h 
//
// This file is distributed as part of the Lorentz++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See LorentzPackage.h for details.

#ifndef ROOT_TThreeVectorComplex
#define ROOT_TThreeVectorComplex
 
#include "Double.h"
#include "TBuffer.h"
#include "Complex.h"
#include "TThreeVectorReal.h"
 
#include <math.h>

class TThreeRotation;
 
 
class TThreeVectorComplex {
 
friend class TLorentzTransform;
friend class TLorentzBoost;
friend class TThreeRotation;
 
protected:
   Complex_t        fVector[4];        // Complex vector allocated on stack
   static LDouble_t fResolution;       // vector resolving "distance"
 
public:
   TThreeVectorComplex() { }
   explicit TThreeVectorComplex(const Complex_t &x,
                               const Complex_t &y,
                               const Complex_t &z);
   explicit TThreeVectorComplex(const Float_t *array);
   explicit TThreeVectorComplex(const LDouble_t *array);
   explicit TThreeVectorComplex(const Complex_t *array);
   TThreeVectorComplex(const TThreeVectorReal &another);
   TThreeVectorComplex(const TThreeVectorComplex &another);
 
   virtual ~TThreeVectorComplex() { }
 
   Complex_t &operator[](const Int_t index) const;
 
   static void SetResolution(const LDouble_t resolution);
   LDouble_t Resolution() const;

   LDouble_t Length() const;        // sqrt of sum of absolute-squares
   LDouble_t LengthSqr() const;        // sum of absolute-squares
   TThreeVectorReal RealPart() const;
   TThreeVectorReal ImagPart() const;
   void GetCartesian(Complex_t &x, Complex_t &y, Complex_t &z) const;
   void GetCartesian(Complex_t *array) const;
   LDouble_t DistanceTo
            (const Complex_t x, const Complex_t y, const Complex_t z) const;
   LDouble_t DistanceTo(const Complex_t *array) const;
   LDouble_t DistanceTo(const TThreeVectorComplex &vec2) const;
 
   TThreeVectorComplex &operator=(const TThreeVectorReal &source);
   TThreeVectorComplex &operator=(const TThreeVectorComplex &source);
   TThreeVectorComplex &operator=(const Float_t *array);
   TThreeVectorComplex &operator=(const LDouble_t *array);
   TThreeVectorComplex &operator=(const Complex_t *array);
   TThreeVectorComplex &operator+=(const TThreeVectorComplex &source);
   TThreeVectorComplex &operator+=(const Float_t *array);
   TThreeVectorComplex &operator+=(const LDouble_t *array);
   TThreeVectorComplex &operator+=(const Complex_t *array);
   TThreeVectorComplex &operator-=(const TThreeVectorComplex &source);
   TThreeVectorComplex &operator-=(const Float_t *array);
   TThreeVectorComplex &operator-=(const LDouble_t *array);
   TThreeVectorComplex &operator-=(const Complex_t *array);
   TThreeVectorComplex &operator*=(const Complex_t &factor);
   TThreeVectorComplex &operator/=(const Complex_t &factor);
 
   Bool_t operator==(const TThreeVectorComplex &other) const;
   Bool_t operator!=(const TThreeVectorComplex &other) const;
 
   TThreeVectorComplex &Zero();
   TThreeVectorComplex &Conj();
   TThreeVectorComplex &SpaceInv();
   TThreeVectorComplex &Normalize(const LDouble_t length);
   TThreeVectorComplex &Rotate(const TThreeRotation &rotOp);
   TThreeVectorComplex &Rotate(const LDouble_t phi,
                               const LDouble_t theta,
                               const LDouble_t psi);
   TThreeVectorComplex &Rotate(const TUnitVector &ahat, const LDouble_t angle);
   TThreeVectorComplex &Cross(const TThreeVectorComplex &other);
   TThreeVectorComplex &Cross(const TThreeVectorComplex &va,
                              const TThreeVectorComplex &vb);
   Complex_t Dot(const TThreeVectorComplex &other);
 
   TThreeVectorComplex operator-() const;
   friend TThreeVectorComplex operator+(const TThreeVectorComplex &v1,
                                        const TThreeVectorComplex &v2);
   friend TThreeVectorComplex operator-(const TThreeVectorComplex &v1,
                                        const TThreeVectorComplex &v2);
   friend TThreeVectorComplex operator*(const TThreeVectorComplex &vec,
                                        const Complex_t &factor);
   friend TThreeVectorComplex operator*(const Complex_t &factor,
                                        const TThreeVectorComplex &vec);
   friend TThreeVectorComplex operator/(const TThreeVectorComplex &vec,
                                        const Complex_t &factor);

   friend TBuffer &operator>>(TBuffer &buf, TThreeVectorComplex *&vec);
   friend TBuffer &operator<<(TBuffer &buf, const TThreeVectorComplex *vec);
   void Print(Option_t *option="");
 
   ClassDef(TThreeVectorComplex,1)  // Complex three-vector class
};

//----- inlines ----------------------------------------------------------------

inline TThreeVectorComplex::TThreeVectorComplex
       (const Complex_t &x, const Complex_t &y, const Complex_t &z)
{
   fVector[1] = x;
   fVector[2] = y;
   fVector[3] = z;
}

inline TThreeVectorComplex::TThreeVectorComplex(const Float_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TThreeVectorComplex::TThreeVectorComplex(const LDouble_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TThreeVectorComplex::TThreeVectorComplex(const Complex_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TThreeVectorComplex::TThreeVectorComplex
       (const TThreeVectorComplex &another)
{
   *this = another;
}

inline TThreeVectorComplex::TThreeVectorComplex
       (const TThreeVectorReal &another)
{
   *this = another;
}
 
inline Complex_t &TThreeVectorComplex::operator[](const Int_t index) const
{
   if (index < 1 || index > 3) {
      Error("TThreeVectorComplex::operator[]","index out of range");
      return (Complex_t &)fVector[1];
   }
   return (Complex_t &)fVector[index];
}

inline void TThreeVectorComplex::SetResolution(const LDouble_t resolution)
{
   fResolution = resolution;
}

inline LDouble_t TThreeVectorComplex::Resolution() const
{
   LDouble_t scale = Length();
   if (scale > 0)
      return fResolution*scale;
   else
      return fResolution;
}
 
inline LDouble_t TThreeVectorComplex::Length() const
{
   return sqrt(LengthSqr());
}

inline LDouble_t TThreeVectorComplex::LengthSqr() const
{
   return LDouble_t(norm(fVector[1]) + norm(fVector[2]) + norm(fVector[3]));
}

inline TThreeVectorReal TThreeVectorComplex::RealPart() const
{
   TThreeVectorReal R;
   R.fVector[1] = real(fVector[1]);
   R.fVector[2] = real(fVector[2]);
   R.fVector[3] = real(fVector[3]);
   return R;
}

inline TThreeVectorReal TThreeVectorComplex::ImagPart() const
{
   TThreeVectorReal I;
   I.fVector[1] = imag(fVector[1]);
   I.fVector[2] = imag(fVector[2]);
   I.fVector[3] = imag(fVector[3]);
   return I;
}

inline void TThreeVectorComplex::GetCartesian
            (Complex_t &x, Complex_t &y, Complex_t &z) const
{
   x = fVector[1];
   y = fVector[2];
   z = fVector[3];
}

inline void TThreeVectorComplex::GetCartesian(Complex_t *array) const
{
   *(array++) = fVector[1];
   *(array++) = fVector[2];
   *array     = fVector[3];
}
 
inline LDouble_t TThreeVectorComplex::DistanceTo
                (const Complex_t x, const Complex_t y, const Complex_t z)
       const
{
   Complex_t xloc(x), yloc(y), zloc(z);
   xloc -= fVector[1];
   yloc -= fVector[2];
   zloc -= fVector[3];
   return sqrt(norm(xloc) + norm(yloc) + norm(zloc));
}

inline LDouble_t TThreeVectorComplex::DistanceTo(const Complex_t *array) const
{
   Complex_t x = fVector[1] - *(array++);
   Complex_t y = fVector[2] - *(array++);
   Complex_t z = fVector[3] - *array;
   return sqrt(norm(x) + norm(y) + norm(z));
}

inline LDouble_t TThreeVectorComplex::DistanceTo
       (const TThreeVectorComplex &vec2) const
{
   Complex_t x = fVector[1] - vec2.fVector[1];
   Complex_t y = fVector[2] - vec2.fVector[2];
   Complex_t z = fVector[3] - vec2.fVector[3];
   return sqrt(norm(x) + norm(y) + norm(z));
}

inline TThreeVectorComplex &TThreeVectorComplex::operator=
                           (const TThreeVectorReal &source)
{
   fVector[1] = source.fVector[1];
   fVector[2] = source.fVector[2];
   fVector[3] = source.fVector[3];
   return *this;
}

inline TThreeVectorComplex
      &TThreeVectorComplex::operator=(const TThreeVectorComplex &source)
{
   fVector[1] = source.fVector[1];
   fVector[2] = source.fVector[2];
   fVector[3] = source.fVector[3];
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator=
                           (const Float_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator=
                           (const LDouble_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator=
                           (const Complex_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator+=
                           (const TThreeVectorComplex &source)
{
   fVector[1] += source.fVector[1];
   fVector[2] += source.fVector[2];
   fVector[3] += source.fVector[3];
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator+=
                           (const Float_t *array)
{
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator+=
                           (const LDouble_t *array)
{
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator+=
                           (const Complex_t *array)
{
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator-=
                           (const TThreeVectorComplex &source)
{
   fVector[1] -= source.fVector[1];
   fVector[2] -= source.fVector[2];
   fVector[3] -= source.fVector[3];
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator-=
                           (const Float_t *array)
{
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator-=
                           (const LDouble_t *array)
{
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator-=
                           (const Complex_t *array)
{
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator*=
                           (const Complex_t &factor)
{
   fVector[1] *= factor;
   fVector[2] *= factor;
   fVector[3] *= factor;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::operator/=
                           (const Complex_t &factor)
{
   fVector[1] /= factor;
   fVector[2] /= factor;
   fVector[3] /= factor;
   return *this;
}
 
inline Bool_t TThreeVectorComplex::operator==
              (const TThreeVectorComplex &other) const
{
   return (DistanceTo(other) < Resolution());
}

inline Bool_t TThreeVectorComplex::operator!=
              (const TThreeVectorComplex &other) const
{
   return !(*this == other);
}

inline TThreeVectorComplex &TThreeVectorComplex::Zero()
{
   fVector[1] = 0;
   fVector[2] = 0;
   fVector[3] = 0;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::Conj()
{
   fVector[1] = conj(fVector[1]);
   fVector[2] = conj(fVector[2]);
   fVector[3] = conj(fVector[3]);
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::SpaceInv()
{
   fVector[1] = -fVector[1];
   fVector[2] = -fVector[2];
   fVector[3] = -fVector[3];
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::Normalize
                           (const LDouble_t length)
{
   LDouble_t r = Length();
   *this *= length/r;
   return *this;
}

inline TThreeVectorComplex &TThreeVectorComplex::Cross
                           (const TThreeVectorComplex &other)
{
   TThreeVectorComplex temp(*this);
   return Cross(temp,other);
}

inline TThreeVectorComplex &TThreeVectorComplex::Cross
       (const TThreeVectorComplex &va, const TThreeVectorComplex &vb)
{
   fVector[1] = va.fVector[2]*vb.fVector[3] - va.fVector[3]*vb.fVector[2];
   fVector[2] = va.fVector[3]*vb.fVector[1] - va.fVector[1]*vb.fVector[3];
   fVector[3] = va.fVector[1]*vb.fVector[2] - va.fVector[2]*vb.fVector[1];
   return *this;
}

inline Complex_t TThreeVectorComplex::Dot(const TThreeVectorComplex &other)
{
   return (fVector[1]*other.fVector[1] +
           fVector[2]*other.fVector[2] +
           fVector[3]*other.fVector[3]);
}
 
inline TThreeVectorComplex TThreeVectorComplex::operator-() const
{
   TThreeVectorComplex minusThis;
   minusThis.fVector[1] = -fVector[1];
   minusThis.fVector[2] = -fVector[2];
   minusThis.fVector[3] = -fVector[3];
   return minusThis;
}

inline TThreeVectorComplex operator+
       (const TThreeVectorComplex &v1, const TThreeVectorComplex &v2)
{
   TThreeVectorComplex result(v1);
   return (result += v2);
}

inline TThreeVectorComplex operator-
       (const TThreeVectorComplex &v1, const TThreeVectorComplex &v2)
{
   TThreeVectorComplex result(v1);
   return (result -= v2);
}

inline TThreeVectorComplex operator*
       (const TThreeVectorComplex &vec, const Complex_t &factor)
{
   TThreeVectorComplex result(vec);
   return (result *= factor);
}

inline TThreeVectorComplex operator*
       (const Complex_t &factor, const TThreeVectorComplex &vec)
{
   TThreeVectorComplex result(vec);
   return (result *= factor);
}

inline TThreeVectorComplex operator/
       (const TThreeVectorComplex &vec, const Complex_t &factor)
{
   TThreeVectorComplex result(vec);
   return (result /= factor);
}

inline TBuffer &operator>>(TBuffer &buf, TThreeVectorComplex *&obj)
{
   for (Int_t i=1; i<4; i++) {
      Double_t real,imag;
      buf >> real >> imag;
      obj->fVector[i] = Complex_t(real,imag);
   }
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TThreeVectorComplex *obj)
{
   for (Int_t i=1; i<4; i++) {
      Double_t real = obj->fVector[i].real();
      Double_t imag = obj->fVector[i].imag();
      buf << real << imag;
   }
   return buf;
}

#endif
