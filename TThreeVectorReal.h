//
// TThreeVectorReal.h 
//
// This file is distributed as part of the Lorentz++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See LorentzPackage.h for details.

#ifndef ROOT_TThreeVectorReal
#define ROOT_TThreeVectorReal
 
#include "Double.h"
#include "TBuffer.h"
#include "TError.h"
 
#include <math.h>

class TThreeVectorReal;
typedef TThreeVectorReal TUnitVector;
typedef TThreeVectorReal TThreeVector;
class TThreeRotation;
 
class TThreeVectorReal {
 
friend class TThreeVectorComplex;
friend class TLorentzTransform;
friend class TLorentzBoost;
friend class TThreeRotation;
 
protected:
   LDouble_t        fVector[4];       // real vector allocated on stack
   static LDouble_t fResolution;      // vector resolving "distance"
 
public:
   TThreeVectorReal() { }
   explicit TThreeVectorReal(const LDouble_t x, const LDouble_t y, const LDouble_t z);
   explicit TThreeVectorReal(const Float_t *array);
   explicit TThreeVectorReal(const LDouble_t *array);
   TThreeVectorReal(const TThreeVectorReal &another);
 
   virtual ~TThreeVectorReal() { }
 
   LDouble_t &operator[](const Int_t index) const;
 
   static void SetResolution(const LDouble_t resolution);
   LDouble_t Resolution() const;

   LDouble_t Length() const;
   LDouble_t LengthSqr() const;
   LDouble_t Rho() const;
   LDouble_t RhoSqr() const;
   LDouble_t Theta() const;
   LDouble_t CosTheta() const;
   LDouble_t Phi() const;
   void GetPolar(LDouble_t &r, LDouble_t &theta, LDouble_t &phi) const;
   void GetCartesian(LDouble_t &x, LDouble_t &y, LDouble_t &z) const;
   void GetCartesian(LDouble_t *array) const;
   LDouble_t DistanceTo(const LDouble_t x,
                       const LDouble_t y,
                       const LDouble_t z) const;
   LDouble_t DistanceTo(const LDouble_t *array) const;
   LDouble_t DistanceTo(const TThreeVectorReal &vec2) const;
 
   TThreeVectorReal &operator=(const TThreeVectorReal &source);
   TThreeVectorReal &operator=(const Float_t *array);
   TThreeVectorReal &operator=(const LDouble_t *array);
   TThreeVectorReal &operator+=(const TThreeVectorReal &source);
   TThreeVectorReal &operator+=(const Float_t *array);
   TThreeVectorReal &operator+=(const LDouble_t *array);
   TThreeVectorReal &operator-=(const TThreeVectorReal &source);
   TThreeVectorReal &operator-=(const Float_t *array);
   TThreeVectorReal &operator-=(const LDouble_t *array);
   TThreeVectorReal &operator*=(const LDouble_t factor);
   TThreeVectorReal &operator/=(const LDouble_t factor);
 
   Bool_t operator==(const TThreeVectorReal &other) const;
   Bool_t operator!=(const TThreeVectorReal &other) const;
 
   TThreeVectorReal &Zero();
   TThreeVectorReal &SpaceInv();
   TThreeVectorReal &Normalize(const LDouble_t length);
   TThreeVectorReal &SetPolar(const LDouble_t r,
                              const LDouble_t theta,
                              const LDouble_t phi);
   TThreeVectorReal &Rotate(const TThreeRotation &rotOp);
   TThreeVectorReal &Rotate(const LDouble_t phi,
                            const LDouble_t theta,
                            const LDouble_t psi);
   TThreeVectorReal &Rotate(const TUnitVector &ahat, const LDouble_t angle);
   TThreeVectorReal &Cross(const TThreeVectorReal &other);
   TThreeVectorReal &Cross(const TThreeVectorReal &va,
                           const TThreeVectorReal &vb);
   LDouble_t Dot(const TThreeVectorReal &other);
 
   TThreeVectorReal operator-() const;
   friend TThreeVectorReal operator+(const TThreeVectorReal &v1,
                                     const TThreeVectorReal &v2);
   friend TThreeVectorReal operator-(const TThreeVectorReal &v1,
                                     const TThreeVectorReal &v2);
   friend TThreeVectorReal operator*(const TThreeVectorReal &vec,
                                     const LDouble_t factor);
   friend TThreeVectorReal operator*(const LDouble_t factor,
                                     const TThreeVectorReal &vec);
   friend TThreeVectorReal operator/(const TThreeVectorReal &vec,
                                     const LDouble_t factor);

   friend TBuffer &operator>>(TBuffer &buf, TThreeVectorReal *&vec);
   friend TBuffer &operator<<(TBuffer &buf, const TThreeVectorReal *vec);
   void Print(Option_t *option="");
 
   ClassDef(TThreeVectorReal,1)  // Real three-vector class
};


//----- inlines ----------------------------------------------------------------
 
inline TThreeVectorReal::TThreeVectorReal
       (const LDouble_t x, const LDouble_t y, const LDouble_t z)
{
   fVector[1] = x;
   fVector[2] = y;
   fVector[3] = z;
}

inline TThreeVectorReal::TThreeVectorReal(const Float_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TThreeVectorReal::TThreeVectorReal(const LDouble_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
}

inline TThreeVectorReal::TThreeVectorReal(const TThreeVectorReal &another)
{
   *this = another;
}
 
inline LDouble_t &TThreeVectorReal::operator[](const Int_t index) const
{
   if (index <= 0 || index > 3) {
      Error("TThreeVectorReal::operator[]","index out of range");
      return (LDouble_t &)fVector[1];
   }
   return (LDouble_t &)fVector[index];
}

inline void TThreeVectorReal::SetResolution(const LDouble_t resolution)
{
   TThreeVectorReal::fResolution = resolution;
}

inline LDouble_t TThreeVectorReal::Resolution() const
{
   LDouble_t scale = Length();
   if (scale > 0)
      return fResolution*scale;
   else
      return fResolution;
}
 
inline LDouble_t TThreeVectorReal::Length() const
{
   return sqrt(LengthSqr());
}

inline LDouble_t TThreeVectorReal::LengthSqr() const
{
   return fVector[1]*fVector[1] +
          fVector[2]*fVector[2] +
          fVector[3]*fVector[3];
}

inline LDouble_t TThreeVectorReal::Rho() const
{
   return sqrt(RhoSqr());
}

inline LDouble_t TThreeVectorReal::RhoSqr() const
{
   return fVector[1]*fVector[1] + fVector[2]*fVector[2];
}

inline LDouble_t TThreeVectorReal::Theta() const
{
   return atan2(Rho(),fVector[3]);
}

inline LDouble_t TThreeVectorReal::CosTheta() const
{
   return fVector[3]/Length();
}

inline LDouble_t TThreeVectorReal::Phi() const
{
   return atan2(fVector[2],fVector[1]);
}

inline void TThreeVectorReal::GetPolar
            (LDouble_t &r, LDouble_t &theta, LDouble_t &phi) const
{
   r = Length();
   theta = Theta();
   phi = Phi();
}

inline void TThreeVectorReal::GetCartesian
            (LDouble_t &x, LDouble_t &y, LDouble_t &z) const
{
   x = fVector[1];
   y = fVector[2];
   z = fVector[3];
}

inline void TThreeVectorReal::GetCartesian(LDouble_t *array) const
{
   *(array++) = fVector[1];
   *(array++) = fVector[2];
   *array     = fVector[3];
}
 
inline LDouble_t TThreeVectorReal::DistanceTo
                (const LDouble_t x, const LDouble_t y, const LDouble_t z) const
{
   LDouble_t xloc(x), yloc(y), zloc(z);
   xloc -= fVector[1];
   yloc -= fVector[2];
   zloc -= fVector[3];
   return sqrt(xloc*xloc + yloc*yloc + zloc*zloc);
}

inline LDouble_t TThreeVectorReal::DistanceTo(const LDouble_t *array) const
{
   LDouble_t x = fVector[1] - *(array++);
   LDouble_t y = fVector[2] - *(array++);
   LDouble_t z = fVector[3] - *array;
   return sqrt(x*x + y*y + z*z);
}

inline LDouble_t TThreeVectorReal::DistanceTo
                (const TThreeVectorReal &vec2) const
{
   LDouble_t x = fVector[1] - vec2.fVector[1];
   LDouble_t y = fVector[2] - vec2.fVector[2];
   LDouble_t z = fVector[3] - vec2.fVector[3];
   return sqrt(x*x + y*y + z*z);
}

inline TThreeVectorReal &TThreeVectorReal::operator=
                        (const TThreeVectorReal &source)
{
   fVector[1] = source.fVector[1];
   fVector[2] = source.fVector[2];
   fVector[3] = source.fVector[3];
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator=(const Float_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator=(const LDouble_t *array)
{
   fVector[1] = *(array++);
   fVector[2] = *(array++);
   fVector[3] = *array;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator+=
                        (const TThreeVectorReal &source)
{
   fVector[1] += source.fVector[1];
   fVector[2] += source.fVector[2];
   fVector[3] += source.fVector[3];
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator+=(const Float_t *array)
{
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator+=(const LDouble_t *array)
{
   fVector[1] += *(array++);
   fVector[2] += *(array++);
   fVector[3] += *array;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator-=
                        (const TThreeVectorReal &source)
{
   fVector[1] -= source.fVector[1];
   fVector[2] -= source.fVector[2];
   fVector[3] -= source.fVector[3];
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator-=(const Float_t *array)
{
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator-=(const LDouble_t *array)
{
   fVector[1] -= *(array++);
   fVector[2] -= *(array++);
   fVector[3] -= *array;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator*=(const LDouble_t factor)
{
   fVector[1] *= factor;
   fVector[2] *= factor;
   fVector[3] *= factor;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::operator/=(const LDouble_t factor)
{
   fVector[1] /= factor;
   fVector[2] /= factor;
   fVector[3] /= factor;
   return *this;
}
 
inline Bool_t TThreeVectorReal::operator==(const TThreeVectorReal &other) const
{
   return (DistanceTo(other) < Resolution());
}

inline Bool_t TThreeVectorReal::operator!=(const TThreeVectorReal &other) const
{
   return !(*this == other);
}

inline TThreeVectorReal &TThreeVectorReal::Zero()
{
   fVector[1] = 0;
   fVector[2] = 0;
   fVector[3] = 0;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::SpaceInv()
{
   fVector[1] = -fVector[1];
   fVector[2] = -fVector[2];
   fVector[3] = -fVector[3];
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::Normalize(const LDouble_t length)
{
   LDouble_t r = Length();
   *this *= length/r;
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::SetPolar(const LDouble_t r,
                                                    const LDouble_t theta,
                                                    const LDouble_t phi)
{
   fVector[1] = r*sin(theta);
   fVector[2] = fVector[1]*sin(phi);
   fVector[3] = r*cos(theta);
   fVector[1] *= cos(phi);
   return *this;
}

inline TThreeVectorReal &TThreeVectorReal::Cross(const TThreeVectorReal &other)
{
   TThreeVectorReal temp(*this);
   return Cross(temp,other);
}

inline TThreeVectorReal &TThreeVectorReal::Cross
                        (const TThreeVectorReal &va, const TThreeVectorReal &vb)
{
   fVector[1] = va.fVector[2]*vb.fVector[3] - va.fVector[3]*vb.fVector[2];
   fVector[2] = va.fVector[3]*vb.fVector[1] - va.fVector[1]*vb.fVector[3];
   fVector[3] = va.fVector[1]*vb.fVector[2] - va.fVector[2]*vb.fVector[1];
   return *this;
}

inline LDouble_t TThreeVectorReal::Dot(const TThreeVectorReal &other)
{
   return (fVector[1]*other.fVector[1] +
           fVector[2]*other.fVector[2] +
           fVector[3]*other.fVector[3]);
}
 
inline TThreeVectorReal TThreeVectorReal::operator-() const
{
   TThreeVectorReal minusThis;
   minusThis.fVector[1] = -fVector[1];
   minusThis.fVector[2] = -fVector[2];
   minusThis.fVector[3] = -fVector[3];
   return minusThis;
}

inline TThreeVectorReal operator+
                        (const TThreeVectorReal &v1, const TThreeVectorReal &v2)
{
   TThreeVectorReal result(v1);
   return (result += v2);
}

inline TThreeVectorReal operator-
                        (const TThreeVectorReal &v1, const TThreeVectorReal &v2)
{
   TThreeVectorReal result(v1);
   return (result -= v2);
}

inline TThreeVectorReal operator*
                        (const TThreeVectorReal &vec, const LDouble_t factor)
{
   TThreeVectorReal result(vec);
   return (result *= factor);
}

inline TThreeVectorReal operator*
                        (const LDouble_t factor, const TThreeVectorReal &vec)
{
   TThreeVectorReal result(vec);
   return (result *= factor);
}

inline TThreeVectorReal operator/
                        (const TThreeVectorReal &vec, const LDouble_t factor)
{
   TThreeVectorReal result(vec);
   return (result /= factor);
}

inline TBuffer &operator>>(TBuffer &buf, TThreeVectorReal *&obj)
{
   Double_t vector[3];
   buf.ReadStaticArray(vector);
   obj->fVector[1] = vector[0];
   obj->fVector[2] = vector[1];
   obj->fVector[3] = vector[2];
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TThreeVectorReal *obj)
{
   Double_t vector[3];
   vector[0] = obj->fVector[1];
   vector[1] = obj->fVector[2];
   vector[2] = obj->fVector[3];
   buf.WriteArray(vector, 3);
   return buf;
}

#endif
