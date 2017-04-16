//
// TDiracSpinor.h 
//
// This file is distributed as part of the Dirac++ package,
// a general toolkit for computing the amplitudes for Feynman
// graphs. See DiracPackage.h for details.

#ifndef ROOT_TDiracSpinor
#define ROOT_TDiracSpinor
 
#include "Double.h"
#include "TBuffer.h"
#include "TFourVectorComplex.h"
#include "TPauliSpinor.h"
#include "TError.h"
 
#include <math.h>

class TDiracMatrix;
class TLorentzTransform;
class TLorentzBoost;
class TThreeRotation;


class TDiracSpinor {
 
friend class TDiracMatrix;
friend class TPauliSpinor;
 
protected:
   Complex_t   fSpinor[4];        // complex vector allocated on stack
   static LDouble_t fResolution;  // resolution "distance" between objects
 
public:
   TDiracSpinor() { }
   explicit TDiracSpinor(const Complex_t &a1, const Complex_t &a2,
                         const Complex_t &a3, const Complex_t &a4);
   explicit TDiracSpinor(const Complex_t *array);
   explicit TDiracSpinor(const TFourVectorReal &p, const Float_t helicity);
   explicit TDiracSpinor(const TFourVectorReal &p, const TUnitVector &polar);
   TDiracSpinor(const TDiracSpinor &another);
 
   virtual ~TDiracSpinor() { }
 
   Complex_t &operator[](const Int_t index) const;
 
   void static SetResolution(const LDouble_t resolution);
   LDouble_t Resolution() const;

   LDouble_t Norm() const;
   LDouble_t NormSqr() const;
   TPauliSpinor Upper() const;
   TPauliSpinor Lower() const;
   LDouble_t DistanceTo(const TDiracSpinor &another) const;
   LDouble_t DistanceTo(const Complex_t *array) const;
 
   TDiracSpinor &operator=(const TDiracSpinor &source);
   TDiracSpinor &operator=(const Complex_t *array);
   TDiracSpinor &operator+=(const TDiracSpinor &source);
   TDiracSpinor &operator+=(const Complex_t *array);
   TDiracSpinor &operator-=(const TDiracSpinor &source);
   TDiracSpinor &operator-=(const Complex_t *array);
   TDiracSpinor &operator*=(const Complex_t factor);
   TDiracSpinor &operator/=(const Complex_t factor);
 
   Bool_t operator==(const TDiracSpinor &other) const;
   Bool_t operator==(const Complex_t *array) const;
   friend Bool_t operator==(const Complex_t *array, const TDiracSpinor &v1);
   Bool_t operator!=(const TDiracSpinor &other) const;
   Bool_t operator!=(const Complex_t *array) const;
   friend Bool_t operator!=(const Complex_t *array, const TDiracSpinor &v1);
 
   TDiracSpinor &Zero();
   TDiracSpinor &Conj();
   TDiracSpinor &Bar();
   TDiracSpinor SetUpper(TPauliSpinor &phi);
   TDiracSpinor SetLower(TPauliSpinor &chi);
   TDiracSpinor &Normalize(const LDouble_t &norm);
   TDiracSpinor &Normalize(const TFourVectorReal &p);
   TDiracSpinor &SetStateU(const TFourVectorReal &p,
                           const Float_t helicity);
   TDiracSpinor &SetStateV(const TFourVectorReal &p,
                           const Float_t helicity);
   TDiracSpinor &SetStateU(const TFourVectorReal &p,
                           const TUnitVector &polar);
   TDiracSpinor &SetStateV(const TFourVectorReal &p,
                           const TUnitVector &polar);
   TDiracSpinor &Operate(const TDiracMatrix &dmOp);
   TDiracSpinor &Rotate(const TThreeRotation &rotOp);
   TDiracSpinor &Rotate(const LDouble_t &phi,
                        const LDouble_t &theta,
                        const LDouble_t &psi);
   TDiracSpinor &Rotate(const TThreeVectorReal &axis);
   TDiracSpinor &Rotate(const TUnitVector &axis, const LDouble_t angle);
   TDiracSpinor &Boost(const TLorentzBoost &boostOp);
   TDiracSpinor &Boost(const LDouble_t betaX,
                       const LDouble_t betaY,
                       const LDouble_t betaZ);
   TDiracSpinor &Boost(const LDouble_t *beta);
   TDiracSpinor &Boost(const TThreeVectorReal &beta);
   TDiracSpinor &Boost(const TUnitVector &bhat, const LDouble_t beta);
   TDiracSpinor &BoostToRest(const TFourVector &p);
   TDiracSpinor &BoostFromRest(const TFourVector &p);
   Complex_t InnerProd(const TDiracSpinor &other);
   Complex_t ScalarProd(const TDiracSpinor &other);
 
   TDiracSpinor operator-() const;
   friend TDiracSpinor operator+(const TDiracSpinor &v1,
                                 const TDiracSpinor &v2);
   friend TDiracSpinor operator+(const TDiracSpinor &v1,
                                 const Complex_t *a2);
   friend TDiracSpinor operator+(const Complex_t *a1,
                                 const TDiracSpinor &v2);
   friend TDiracSpinor operator-(const TDiracSpinor &v1,
                                 const TDiracSpinor &v2);
   friend TDiracSpinor operator-(const TDiracSpinor &v1,
                                 const Complex_t *a2);
   friend TDiracSpinor operator-(const Complex_t *a1,
                                 const TDiracSpinor &v2);
   friend TDiracSpinor operator*(const TDiracSpinor &vec,
                                 const LDouble_t &factor);
   friend TDiracSpinor operator*(const LDouble_t &factor,
                                 const TDiracSpinor &vec);
   friend TDiracSpinor operator*(const TDiracSpinor &vec,
                                 const Complex_t &factor);
   friend TDiracSpinor operator*(const Complex_t &factor,
                                 const TDiracSpinor &vec);
   friend TDiracSpinor operator*(const TDiracMatrix &dmOp,
                                 const TDiracSpinor &vec);
   friend TDiracSpinor operator/(const TDiracSpinor &vec,
                                 const Complex_t &factor);

   friend TBuffer &operator>>(TBuffer &buf, TDiracSpinor *&obj);
   friend TBuffer &operator<<(TBuffer &buf, const TDiracSpinor *obj);
   void Print(Option_t *option="");
 
   ClassDef(TDiracSpinor,1)  // Complex Dirac spinor class
};


//----- inlines ----------------------------------------------------------------
 
inline TDiracSpinor::TDiracSpinor(const Complex_t &a1, const Complex_t &a2,
                                  const Complex_t &a3, const Complex_t &a4)
{
   fSpinor[0] = a1;
   fSpinor[1] = a2;
   fSpinor[2] = a3;
   fSpinor[3] = a4;
}

inline TDiracSpinor::TDiracSpinor(const Complex_t *array)
{
   fSpinor[0] = *(array++);
   fSpinor[1] = *(array++);
   fSpinor[2] = *(array++);
   fSpinor[3] = *array;
}

inline TDiracSpinor::TDiracSpinor
       (const TFourVectorReal &p, const Float_t helicity)
{
   SetStateU(p,helicity);
}

inline TDiracSpinor::TDiracSpinor
       (const TFourVectorReal &p, const TUnitVector &polar)
{
   SetStateU(p,polar);
}

inline TDiracSpinor::TDiracSpinor(const TDiracSpinor &another)
{
   *this = another;
}
 
inline Complex_t &TDiracSpinor::operator[](const Int_t index) const
{
   if (index < 0 || index > 3) {
      Error("TDiracSpinor::operator[]","index out of range");
      return (Complex_t &)fSpinor[0];
   }
   return (Complex_t &)fSpinor[index];
}
 
inline void TDiracSpinor::SetResolution(const LDouble_t resolution)
{
   fResolution = resolution;
}

inline LDouble_t TDiracSpinor::Resolution() const
{
   LDouble_t scale = Norm();
   if (scale > 0)
      return fResolution*scale;
   else
      return fResolution;
}

inline LDouble_t TDiracSpinor::Norm() const
{
   return sqrt(NormSqr());
}

inline LDouble_t TDiracSpinor::NormSqr() const
{
   return (norm(fSpinor[0]) + norm(fSpinor[1]) +
           norm(fSpinor[2]) + norm(fSpinor[3]) );
}

inline LDouble_t TDiracSpinor::DistanceTo(const TDiracSpinor &another) const
{
   return sqrt(norm(fSpinor[0] - another.fSpinor[0]) +
               norm(fSpinor[1] - another.fSpinor[1]) +
               norm(fSpinor[2] - another.fSpinor[2]) +
               norm(fSpinor[3] - another.fSpinor[3]) );
}

inline LDouble_t TDiracSpinor::DistanceTo(const Complex_t *array) const
{
   return sqrt(norm(fSpinor[0] - array[0]) +
               norm(fSpinor[1] - array[1]) +
               norm(fSpinor[2] - array[2]) +
               norm(fSpinor[3] - array[3]) );
}

inline TDiracSpinor &TDiracSpinor::operator=(const TDiracSpinor &source)
{
   fSpinor[0] = source.fSpinor[0];
   fSpinor[1] = source.fSpinor[1];
   fSpinor[2] = source.fSpinor[2];
   fSpinor[3] = source.fSpinor[3];
   return *this;
}

inline TDiracSpinor &TDiracSpinor::operator=(const Complex_t *array)
{
   fSpinor[0] = *(array++);
   fSpinor[1] = *(array++);
   fSpinor[2] = *(array++);
   fSpinor[3] = *array;
   return *this;
}

inline TDiracSpinor &TDiracSpinor::operator+=(const TDiracSpinor &source)
{
   fSpinor[0] += source.fSpinor[0];
   fSpinor[1] += source.fSpinor[1];
   fSpinor[2] += source.fSpinor[2];
   fSpinor[3] += source.fSpinor[3];
   return *this;
}

inline TDiracSpinor &TDiracSpinor::operator+=(const Complex_t *array)
{
   fSpinor[0] += *(array++);
   fSpinor[1] += *(array++);
   fSpinor[2] += *(array++);
   fSpinor[3] += *array;
   return *this;
}

inline TDiracSpinor &TDiracSpinor::operator-=(const TDiracSpinor &source)
{
   fSpinor[0] -= source.fSpinor[0];
   fSpinor[1] -= source.fSpinor[1];
   fSpinor[2] -= source.fSpinor[2];
   fSpinor[3] -= source.fSpinor[3];
   return *this;
}

inline TDiracSpinor &TDiracSpinor::operator-=(const Complex_t *array)
{
   fSpinor[0] -= *(array++);
   fSpinor[1] -= *(array++);
   fSpinor[2] -= *(array++);
   fSpinor[3] -= *array;
   return *this;
}

inline TDiracSpinor &TDiracSpinor::operator*=(const Complex_t factor)
{
   fSpinor[0] *= factor;
   fSpinor[1] *= factor;
   fSpinor[2] *= factor;
   fSpinor[3] *= factor;
   return *this;
}

inline TDiracSpinor &TDiracSpinor::operator/=(const Complex_t factor)
{
   fSpinor[0] /= factor;
   fSpinor[1] /= factor;
   fSpinor[2] /= factor;
   fSpinor[3] /= factor;
   return *this;
}

inline Bool_t TDiracSpinor::operator==(const TDiracSpinor &other) const
{
   return (DistanceTo(other) < Resolution());
}

inline Bool_t TDiracSpinor::operator==(const Complex_t *array) const
{
   return (DistanceTo(array) < Resolution());
}

inline Bool_t operator==(const Complex_t *array, const TDiracSpinor &v1)
{
   return (v1.DistanceTo(array) < v1.Resolution());
}

inline Bool_t TDiracSpinor::operator!=(const TDiracSpinor &other) const
{
   return !(*this == other);
}

inline Bool_t TDiracSpinor::operator!=(const Complex_t *array) const
{
   return !(*this == array);
}

inline Bool_t operator!=(const Complex_t *array, const TDiracSpinor &v1)
{
   return !(v1 == array);
}
 
inline TDiracSpinor &TDiracSpinor::Zero()
{
   fSpinor[0] = fSpinor[1] = fSpinor[2] = fSpinor[3] = 0;
   return *this;
}

inline TDiracSpinor &TDiracSpinor::Conj()
{
   fSpinor[0] = conj(fSpinor[0]);
   fSpinor[1] = conj(fSpinor[1]);
   fSpinor[2] = conj(fSpinor[2]);
   fSpinor[3] = conj(fSpinor[3]);
   return *this;
}

inline TDiracSpinor &TDiracSpinor::Bar()
{
   fSpinor[0] = conj(fSpinor[0]);
   fSpinor[1] = conj(fSpinor[1]);
   fSpinor[2] = -conj(fSpinor[2]);
   fSpinor[3] = -conj(fSpinor[3]);
   return *this;
}

inline TDiracSpinor &TDiracSpinor::Normalize(const LDouble_t &norm)
{
   return (*this *= norm/Norm());
}

inline TDiracSpinor &TDiracSpinor::Normalize(const TFourVectorReal &p)
{
// Multiply the spinor by a real number that adjusts the norm of the state
// to its standard relativistic value of sqrt(2E).  Note that this standard
// convention differs from the sqrt(E/m) used by Bjorken and Drell.  The B&D
// choice suffers from the fact that it is ill-defined for m = 0.

   return Normalize(sqrt(2*abs(p[0])));
}

inline Complex_t TDiracSpinor::InnerProd(const TDiracSpinor &other)
{
   return ( conj(fSpinor[0])*other.fSpinor[0] +
            conj(fSpinor[1])*other.fSpinor[1] +
            conj(fSpinor[2])*other.fSpinor[2] +
            conj(fSpinor[3])*other.fSpinor[3] );
}

inline TDiracSpinor TDiracSpinor::operator-() const
{
   TDiracSpinor result;
   result.fSpinor[0] = -fSpinor[0];
   result.fSpinor[1] = -fSpinor[1];
   result.fSpinor[2] = -fSpinor[2];
   result.fSpinor[3] = -fSpinor[3];
   return result;
}

inline TDiracSpinor operator+(const TDiracSpinor &v1, const TDiracSpinor &v2)
{
   TDiracSpinor result(v1);
   result.fSpinor[0] += v2.fSpinor[0];
   result.fSpinor[1] += v2.fSpinor[1];
   result.fSpinor[2] += v2.fSpinor[2];
   result.fSpinor[3] += v2.fSpinor[3];
   return result;
}

inline TDiracSpinor operator+(const TDiracSpinor &v1, const Complex_t *a2)
{
   TDiracSpinor result(v1);
   result.fSpinor[0] += *(a2++);
   result.fSpinor[1] += *(a2++);
   result.fSpinor[2] += *(a2++);
   result.fSpinor[3] += *a2;
   return result;
}

inline TDiracSpinor operator+(const Complex_t *a1, const TDiracSpinor &v2)
{
   TDiracSpinor result(v2);
   result.fSpinor[0] += *(a1++);
   result.fSpinor[1] += *(a1++);
   result.fSpinor[2] += *(a1++);
   result.fSpinor[3] += *a1;
   return result;
}

inline TDiracSpinor operator-(const TDiracSpinor &v1, const TDiracSpinor &v2)
{
   TDiracSpinor result(v1);
   result.fSpinor[0] -= v2.fSpinor[0];
   result.fSpinor[1] -= v2.fSpinor[1];
   result.fSpinor[2] -= v2.fSpinor[2];
   result.fSpinor[3] -= v2.fSpinor[3];
   return result;
}

inline TDiracSpinor operator-(const TDiracSpinor &v1, const Complex_t *a2)
{
   TDiracSpinor result(v1);
   result.fSpinor[0] -= *(a2++);
   result.fSpinor[1] -= *(a2++);
   result.fSpinor[2] -= *(a2++);
   result.fSpinor[3] -= *a2;
   return result;
}

inline TDiracSpinor operator-(const Complex_t *a1, const TDiracSpinor &v2)
{
   TDiracSpinor result(v2);
   result.fSpinor[0] -= *(a1++);
   result.fSpinor[1] -= *(a1++);
   result.fSpinor[2] -= *(a1++);
   result.fSpinor[3] -= *a1;
   return result;
}

inline TDiracSpinor operator*(const TDiracSpinor &vec, const LDouble_t &factor)
{
   TDiracSpinor result(vec);
   result.fSpinor[0] *= factor;
   result.fSpinor[1] *= factor;
   result.fSpinor[2] *= factor;
   result.fSpinor[3] *= factor;
   return result;
}

inline TDiracSpinor operator*(const LDouble_t &factor, const TDiracSpinor &vec)
{
   TDiracSpinor result(vec);
   result.fSpinor[0] *= factor;
   result.fSpinor[1] *= factor;
   result.fSpinor[2] *= factor;
   result.fSpinor[3] *= factor;
   return result;
}

inline TDiracSpinor operator*(const TDiracSpinor &vec, const Complex_t &factor)
{
   TDiracSpinor result(vec);
   result.fSpinor[0] *= factor;
   result.fSpinor[1] *= factor;
   result.fSpinor[2] *= factor;
   result.fSpinor[3] *= factor;
   return result;
}

inline TDiracSpinor operator*(const Complex_t &factor, const TDiracSpinor &vec)
{
   TDiracSpinor result(vec);
   result.fSpinor[0] *= factor;
   result.fSpinor[1] *= factor;
   result.fSpinor[2] *= factor;
   result.fSpinor[3] *= factor;
   return result;
}

inline TDiracSpinor operator/(const TDiracSpinor &vec, const Complex_t &factor)
{
   TDiracSpinor result(vec);
   result.fSpinor[0] /= factor;
   result.fSpinor[1] /= factor;
   result.fSpinor[2] /= factor;
   result.fSpinor[3] /= factor;
   return result;
}

inline TBuffer &operator>>(TBuffer &buf, TDiracSpinor *&obj)
{
   for (Int_t i=0; i<4; i++) {
      Double_t real,imag;
      buf >> real >> imag;
      obj->fSpinor[i] = Complex_t(real,imag);
   }
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TDiracSpinor *obj)
{
   for (Int_t i=0; i<4; i++) {
      Double_t real = obj->fSpinor[i].real();
      Double_t imag = obj->fSpinor[i].imag();
      buf << real << imag;
   }
   return buf;
}

#endif
