//
// TFourVectorReal.h 
//
// This file is distributed as part of the Lorentz++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See LorentzPackage.h for details.

#ifndef ROOT_TFourVectorReal
#define ROOT_TFourVectorReal
 
#include "TBuffer.h"
#include "TThreeVectorReal.h"
#include "TError.h"
 
#include <math.h>

class TFourVectorReal;
typedef TFourVectorReal TFourVector;

class TFourVectorComplex;
class TLorentzTransform;
class TLorentzBoost;
class TThreeRotation;


class TFourVectorReal : public TThreeVectorReal {

friend class TFourVectorComplex;
friend class TLorentzTransform;
friend class TLorentzBoost;
friend class TThreeRotation;
 
public:
   TFourVectorReal() : TThreeVectorReal() { }
   TFourVectorReal(const Double_t t, const Double_t x,
                   const Double_t y, const Double_t z);
   TFourVectorReal(const Float_t *array);
   TFourVectorReal(const Double_t *array);
   TFourVectorReal(const TFourVectorReal &another);
   TFourVectorReal(const Double_t t, const TThreeVectorReal &r);
 
   virtual ~TFourVectorReal() { }

   Double_t Resolution() const;
 
   Double_t &operator[](const Int_t index) const;
 
   Double_t Invariant() const;		// uses the +--- metric
   Double_t InvariantSqr() const;
   void GetCoord(Double_t &t, Double_t &x, Double_t &y, Double_t &z) const;
   void GetCoord(Double_t *array) const;
   Double_t DistanceTo(const Double_t t, const Double_t x,
                       const Double_t y, const Double_t z) const;
   Double_t DistanceTo(const Double_t *array) const;
   Double_t DistanceTo(const TFourVectorReal &vec2) const;
 
   TFourVectorReal &operator=(const TFourVectorReal &source);
   TFourVectorReal &operator=(const Float_t *array);
   TFourVectorReal &operator=(const Double_t *array);
   TFourVectorReal &operator+=(const TFourVectorReal &source);
   TFourVectorReal &operator+=(const Float_t *array);
   TFourVectorReal &operator+=(const Double_t *array);
   TFourVectorReal &operator-=(const TFourVectorReal &source);
   TFourVectorReal &operator-=(const Float_t *array);
   TFourVectorReal &operator-=(const Double_t *array);
   TFourVectorReal &operator*=(const Double_t factor);
   TFourVectorReal &operator/=(const Double_t factor);
 
   Bool_t operator==(const TFourVectorReal &other) const;
   Bool_t operator!=(const TFourVectorReal &other) const;
 
   TFourVectorReal &Zero();
   TFourVectorReal &Transform(const TLorentzTransform &xformOp);
   TFourVectorReal &Boost(const TLorentzBoost &boostOp);
   TFourVectorReal &Boost(const Double_t betaX,
                          const Double_t betaY,
                          const Double_t betaZ);
   TFourVectorReal &Boost(const TUnitVector &bhat, const Double_t beta);
   Double_t ScalarProd(const TFourVectorReal &other);
 
   TFourVectorReal operator-() const;
   friend TFourVectorReal operator+(const TFourVectorReal &v1,
                                    const TFourVectorReal &v2);
   friend TFourVectorReal operator-(const TFourVectorReal &v1,
                                    const TFourVectorReal &v2);
   friend TFourVectorReal operator*(const TFourVectorReal &vec,
                                    const Double_t factor);
   friend TFourVectorReal operator*(const Double_t factor,
                                    const TFourVectorReal &vec);
   friend TFourVectorReal operator/(const TFourVectorReal &vec,
                                    const Double_t factor);

   friend TBuffer &operator>>(TBuffer &buf, TFourVectorReal *&vec);
   friend TBuffer &operator<<(TBuffer &buf, const TFourVectorReal *vec);
   void Print(Option_t *option="");
 
   ClassDef(TFourVectorReal,1)  // Real four-vector class
};

//----- inlines ----------------------------------------------------------------
 
inline TFourVectorReal::TFourVectorReal
       (const Double_t t, const Double_t x, const Double_t y, const Double_t z)
 : TThreeVectorReal()
{
   fVector[0] = t;
   fVector[1] = x;
   fVector[2] = y;
   fVector[3] = z;
}

inline TFourVectorReal::TFourVectorReal(const Float_t *array)
 : TThreeVectorReal()
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TFourVectorReal::TFourVectorReal(const Double_t *array)
 : TThreeVectorReal()
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TFourVectorReal::TFourVectorReal(const TFourVectorReal &another)
 : TThreeVectorReal()
{
   *this = another;
}
 
inline TFourVectorReal::TFourVectorReal
                        (const Double_t t, const TThreeVectorReal &another)
 : TThreeVectorReal(another)
{
   fVector[0] = t;
}
 
inline Double_t TFourVectorReal::Resolution() const
{
   Double_t scale = sqrt(fVector[0]*fVector[0] + LengthSqr());
   if (scale > 0)
      return fResolution*scale;
   else
      return fResolution;
}
 
inline Double_t &TFourVectorReal::operator[](const Int_t index) const
{
   if (index < 0 || index > 3) {
      Error("TFourVectorReal::operator[]","index out of range");
      return (Double_t &)fVector[0];
   }
   return (Double_t &)fVector[index];
}
 
inline Double_t TFourVectorReal::Invariant() const
{
   Double_t inv2 = InvariantSqr();
   if (inv2 > 0)
      return sqrt(inv2);
   else if (inv2 > -Resolution())
      return 0;
   else
      Error("TFourVectorReal::Invariant()",
            "invoked on a vector with negative norm");
   return -1;
}

inline Double_t TFourVectorReal::InvariantSqr() const
{
   return (fVector[0]*fVector[0] - LengthSqr());
}

inline void TFourVectorReal::GetCoord
            (Double_t &t, Double_t &x, Double_t &y, Double_t &z) const
{
   t = fVector[0];
   x = fVector[1];
   y = fVector[2];
   z = fVector[3];
}

inline void TFourVectorReal::GetCoord(Double_t *array) const
{
   *(array++) = fVector[0];
   *(array++) = fVector[1];
   *(array++) = fVector[2];
   *array     = fVector[3];
}
 
inline Double_t TFourVectorReal::DistanceTo
                (const Double_t t, const Double_t x,
                 const Double_t y, const Double_t z) const
{
   Double_t tloc(t), xloc(x), yloc(y), zloc(z);
   tloc -= fVector[0];
   xloc -= fVector[1];
   yloc -= fVector[2];
   zloc -= fVector[3];
   return sqrt(tloc*tloc + xloc*xloc + yloc*yloc + zloc*zloc);
}

inline Double_t TFourVectorReal::DistanceTo(const Double_t *array) const
{
   Double_t t = fVector[0] - *(array++);
   Double_t x = fVector[1] - *(array++);
   Double_t y = fVector[2] - *(array++);
   Double_t z = fVector[3] - *array;
   return sqrt(t*t + x*x + y*y + z*z);
}

inline Double_t TFourVectorReal::DistanceTo
                (const TFourVectorReal &vec2) const
{
   Double_t t = fVector[0] - vec2.fVector[0];
   Double_t x = fVector[1] - vec2.fVector[1];
   Double_t y = fVector[2] - vec2.fVector[2];
   Double_t z = fVector[3] - vec2.fVector[3];
   return sqrt(t*t + x*x + y*y + z*z);
}

inline TFourVectorReal &TFourVectorReal::operator=
                       (const TFourVectorReal &source)
{
   fVector[0] = source.fVector[0];
   fVector[1] = source.fVector[1];
   fVector[2] = source.fVector[2];
   fVector[3] = source.fVector[3];
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator=(const Float_t *array)
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator=(const Double_t *array)
{
   fVector[0] = *(array++);
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator+=
                        (const TFourVectorReal &source)
{
   fVector[0] += source.fVector[0];
   fVector[1] += source.fVector[1];
   fVector[2] += source.fVector[2];
   fVector[3] += source.fVector[3];
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator+=(const Float_t *array)
{
   fVector[0] += *(array++);
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator+=(const Double_t *array)
{
   fVector[0] += *(array++);
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator-=
                        (const TFourVectorReal &source)
{
   fVector[0] -= source.fVector[0];
   fVector[1] -= source.fVector[1];
   fVector[2] -= source.fVector[2];
   fVector[3] -= source.fVector[3];
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator-=(const Float_t *array)
{
   fVector[0] -= *(array++);
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator-=(const Double_t *array)
{
   fVector[0] -= *(array++);
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator*=(const Double_t factor)
{
   fVector[0] *= factor;
   fVector[1] *= factor;
   fVector[2] *= factor;
   fVector[3] *= factor;
   return *this;
}

inline TFourVectorReal &TFourVectorReal::operator/=(const Double_t factor)
{
   fVector[0] /= factor;
   fVector[1] /= factor;
   fVector[2] /= factor;
   fVector[3] /= factor;
   return *this;
}
 
inline Bool_t TFourVectorReal::operator==(const TFourVectorReal &other) const
{
   return (DistanceTo(other) < Resolution());
}

inline Bool_t TFourVectorReal::operator!=(const TFourVectorReal &other) const
{
   return !(*this == other);
}

inline TFourVectorReal &TFourVectorReal::Zero()
{
   fVector[0] = 0;
   fVector[1] = 0;
   fVector[2] = 0;
   fVector[3] = 0;
   return *this;
}

inline Double_t TFourVectorReal::ScalarProd(const TFourVectorReal &other)
{
   return  Double_t(fVector[0] * other.fVector[0] -
                    fVector[1] * other.fVector[1] -
                    fVector[2] * other.fVector[2] -
                    fVector[3] * other.fVector[3]);
}

inline TFourVectorReal TFourVectorReal::operator-() const
{
   TFourVectorReal minusThis;
   minusThis.fVector[0] = -fVector[0];
   minusThis.fVector[1] = -fVector[1];
   minusThis.fVector[2] = -fVector[2];
   minusThis.fVector[3] = -fVector[3];
   return minusThis;
}

inline TFourVectorReal operator+
       (const TFourVectorReal &v1, const TFourVectorReal &v2)
{
   TFourVectorReal result(v1);
   return (result += v2);
}

inline TFourVectorReal operator-
       (const TFourVectorReal &v1, const TFourVectorReal &v2)
{
   TFourVectorReal result(v1);
   return (result -= v2);
}

inline TFourVectorReal operator*
       (const TFourVectorReal &vec, const Double_t factor)
{
   TFourVectorReal result(vec);
   return (result *= factor);
}

inline TFourVectorReal operator*
       (const Double_t factor, const TFourVectorReal &vec)
{
   TFourVectorReal result(vec);
   return (result *= factor);
}

inline TFourVectorReal operator/
       (const TFourVectorReal &vec, const Double_t factor)
{
   TFourVectorReal result(vec);
   return (result /= factor);
}

inline TBuffer &operator>>(TBuffer &buf, TFourVectorReal *&obj)
{
   buf.ReadStaticArray(obj->fVector);
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TFourVectorReal *obj)
{
   buf.WriteArray(obj->fVector, 4);
   return buf;
}

#endif
