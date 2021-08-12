//
// TPauliSpinor.h 
//
// This file is distributed as part of the Pauli++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See PauliPackage.h for details.

#ifndef ROOT_TPauliSpinor
#define ROOT_TPauliSpinor
 
#include "Double.h"
#include "TBuffer.h"
#include "TThreeVectorComplex.h"
#include "TError.h"
 
#include <math.h>

class TPauliMatrix;
class TThreeRotation;


class TPauliSpinor {
 
friend class TPauliMatrix;
 
protected:
   Complex_t        fSpinor[2];      // complex state vector allocated on stack
   static LDouble_t fResolution;     // vector resolving "distance"
 
public:
   TPauliSpinor() { }
   explicit TPauliSpinor(const Complex_t &a1, const Complex_t &a2);
   explicit TPauliSpinor(const Complex_t *array);
   explicit TPauliSpinor(const LDouble_t theta, const LDouble_t phi);
   explicit TPauliSpinor(const TUnitVector &pol);
   TPauliSpinor(const TPauliSpinor &another);
 
   virtual ~TPauliSpinor() { }
 
   Complex_t &operator[](Int_t index) const;

   static void SetResolution(const LDouble_t resolution);
   LDouble_t Resolution() const;
 
   LDouble_t Norm() const;
   LDouble_t NormSqr() const;
   void GetPolar(LDouble_t &theta, LDouble_t &phi) const;
   TUnitVector Polar() const;
   LDouble_t DistanceTo(const TPauliSpinor &another) const;
   LDouble_t DistanceTo(const Complex_t *array) const;
 
   TPauliSpinor &operator=(const TPauliSpinor &source);
   TPauliSpinor &operator=(const Complex_t *array);
   TPauliSpinor &operator+=(const TPauliSpinor &source);
   TPauliSpinor &operator+=(const Complex_t *array);
   TPauliSpinor &operator-=(const TPauliSpinor &source);
   TPauliSpinor &operator-=(const Complex_t *array);
   TPauliSpinor &operator*=(const Complex_t &factor);
   TPauliSpinor &operator/=(const Complex_t &factor);
 
   Bool_t operator==(const TPauliSpinor &other) const;
   Bool_t operator==(const Complex_t *array) const;
   friend Bool_t operator==(const Complex_t *array, const TPauliSpinor &v1);
   Bool_t operator!=(const TPauliSpinor &other) const;
   Bool_t operator!=(const Complex_t *array) const;
   friend Bool_t operator!=(const Complex_t *array, const TPauliSpinor &v1);
 
   TPauliSpinor &Zero();
   TPauliSpinor &Conj();
   TPauliSpinor &Normalize();
   TPauliSpinor &Normalize(const LDouble_t &norm);
   TPauliSpinor &SetPolar(const LDouble_t &theta, const LDouble_t &phi);
   TPauliSpinor &SetPolar(const TUnitVector &pol);
   TPauliSpinor &Operate(const TPauliMatrix &xOp);
   TPauliSpinor &Rotate(const TThreeRotation &rotOp);
   TPauliSpinor &Rotate(const LDouble_t &phi,
                        const LDouble_t &theta,
                        const LDouble_t &psi);
   TPauliSpinor &Rotate(const TThreeVectorReal &axis);
   TPauliSpinor &Rotate(const TUnitVector &axis, const LDouble_t angle);
   Complex_t InnerProd(const TPauliSpinor &other);
   Complex_t ScalarProd(const TPauliSpinor &other);
 
   TPauliSpinor operator-() const;
   friend TPauliSpinor operator+(const TPauliSpinor &v1,
                                 const TPauliSpinor &v2);
   friend TPauliSpinor operator+(const TPauliSpinor &v1,
                                 const Complex_t *a2);
   friend TPauliSpinor operator+(const Complex_t *a1,
                                 const TPauliSpinor &v2);
   friend TPauliSpinor operator-(const TPauliSpinor &v1,
                                 const TPauliSpinor &v2);
   friend TPauliSpinor operator-(const TPauliSpinor &v1,
                                 const Complex_t *a2);
   friend TPauliSpinor operator-(const Complex_t *a1,
                                 const TPauliSpinor &v2);
   friend TPauliSpinor operator*(const TPauliSpinor &vec,
                                 const LDouble_t &factor);
   friend TPauliSpinor operator*(const LDouble_t &factor,
                                 const TPauliSpinor &vec);
   friend TPauliSpinor operator*(const TPauliSpinor &vec,
                                 const Complex_t &factor);
   friend TPauliSpinor operator*(const Complex_t &factor,
                                 const TPauliSpinor &vec);
   friend TPauliSpinor operator*(const TPauliMatrix &pmOp,
                                 const TPauliSpinor &vec);
   friend TPauliSpinor operator/(const TPauliSpinor &vec,
                                 const Complex_t &factor);

   friend TBuffer &operator>>(TBuffer &buf, TPauliSpinor *&obj);
   friend TBuffer &operator<<(TBuffer &buf, const TPauliSpinor *obj);
   void Print(Option_t *option="");
 
   ClassDef(TPauliSpinor,1)  // Complex Pauli spinor class
};

//----- inlines ----------------------------------------------------------------
 
inline  TPauliSpinor::TPauliSpinor(const Complex_t &a1, const Complex_t &a2)
{
   fSpinor[0] = a1;
   fSpinor[1] = a2;
}

inline  TPauliSpinor::TPauliSpinor(const Complex_t *array)
{
   fSpinor[0] = *(array++);
   fSpinor[1] = *array;
}
   
inline  TPauliSpinor::TPauliSpinor(const LDouble_t theta, const LDouble_t phi)
{
   SetPolar(theta,phi);
}

inline  TPauliSpinor::TPauliSpinor(const TUnitVector &pol)
{
   SetPolar(pol);
}

inline  TPauliSpinor::TPauliSpinor(const TPauliSpinor &another)
{
   *this = another;
}
 
inline Complex_t &TPauliSpinor::operator[](Int_t index) const
{
   if (index < 0 || index > 1) {
      Error("TPauliSpinor::operator[]","index out of range");
      return (Complex_t &)fSpinor[0];
   }
   return (Complex_t &)fSpinor[index];
}

inline void TPauliSpinor::SetResolution(const LDouble_t resolution)
{
   fResolution = resolution;
}

inline LDouble_t TPauliSpinor::Resolution() const
{
   LDouble_t scale = Norm();
   if (scale > 0)
      return fResolution*scale;
   else
      return fResolution;
}

inline LDouble_t TPauliSpinor::Norm() const
{
   return sqrt(NormSqr());
}

inline LDouble_t TPauliSpinor::NormSqr() const
{
   return (norm(fSpinor[0]) + norm(fSpinor[1]));
}

inline void TPauliSpinor::GetPolar(LDouble_t &theta, LDouble_t &phi) const
{
   LDouble_t aUp = abs(fSpinor[0]);
   LDouble_t aDown = abs(fSpinor[1]);
   theta = 2*atan2(aDown,aUp);
   if (aDown < Resolution())
      phi = -2*arg(fSpinor[0]);
   else if (aUp < Resolution())
      phi = 2*arg(fSpinor[1]);
   else
      phi = arg(fSpinor[1]/fSpinor[0]);
}

inline TUnitVector TPauliSpinor::Polar() const
{
   TUnitVector pol;
   LDouble_t theta=0, phi=0;
   GetPolar(theta,phi);
   return pol.SetPolar(1,theta,phi);
}
 
inline LDouble_t TPauliSpinor::DistanceTo(const TPauliSpinor &another) const
{
   return sqrt(norm(fSpinor[0] - another.fSpinor[0]) +
               norm(fSpinor[1] - another.fSpinor[1]) );
}

inline LDouble_t TPauliSpinor::DistanceTo(const Complex_t *array) const
{
   return sqrt(norm(fSpinor[0] - array[0]) +
               norm(fSpinor[1] - array[1]) );
}

inline TPauliSpinor &TPauliSpinor::operator=(const TPauliSpinor &source)
{
   fSpinor[0] = source.fSpinor[0];
   fSpinor[1] = source.fSpinor[1];
   return *this;
}

inline TPauliSpinor &TPauliSpinor::operator=(const Complex_t *array)
{
   fSpinor[0] = *(array++);
   fSpinor[1] = *array;
   return *this;
}

inline TPauliSpinor &TPauliSpinor::operator+=(const TPauliSpinor &source)
{
   fSpinor[0] += source.fSpinor[0];
   fSpinor[1] += source.fSpinor[1];
   return *this;
}

inline TPauliSpinor &TPauliSpinor::operator+=(const Complex_t *array)
{
   fSpinor[0] += *(array++);
   fSpinor[1] += *array;
   return *this;
}

inline TPauliSpinor &TPauliSpinor::operator-=(const TPauliSpinor &source)
{
   fSpinor[0] -= source.fSpinor[0];
   fSpinor[1] -= source.fSpinor[1];
   return *this;
}

inline TPauliSpinor &TPauliSpinor::operator-=(const Complex_t *array)
{
   fSpinor[0] -= *(array++);
   fSpinor[1] -= *array;
   return *this;
}

inline TPauliSpinor &TPauliSpinor::operator*=(const Complex_t &factor)
{
   fSpinor[0] *= factor;
   fSpinor[1] *= factor;
   return *this;
}

inline TPauliSpinor &TPauliSpinor::operator/=(const Complex_t &factor)
{
   fSpinor[0] /= factor;
   fSpinor[1] /= factor;
   return *this;
}

inline Bool_t TPauliSpinor::operator==(const TPauliSpinor &other) const
{
   return (DistanceTo(other) < Resolution());
}

inline Bool_t TPauliSpinor::operator==(const Complex_t *array) const
{
   return (DistanceTo(array) < Resolution());
}

inline Bool_t operator==(const Complex_t *array, const TPauliSpinor &v1)
{
   return (v1.DistanceTo(array) < v1.Resolution());
}

inline Bool_t TPauliSpinor::operator!=(const TPauliSpinor &other) const
{
   return !(*this == other);
}

inline Bool_t TPauliSpinor::operator!=(const Complex_t *array) const
{
   return !(*this == array);
}

inline Bool_t operator!=(const Complex_t *array, const TPauliSpinor &v1)
{
   return !(v1 == array);
}
 
inline TPauliSpinor &TPauliSpinor::Zero()
{
   fSpinor[0] = fSpinor[1] = 0;
   return *this;
}

inline TPauliSpinor &TPauliSpinor::Conj()
{
   fSpinor[0] = conj(fSpinor[0]);
   fSpinor[1] = conj(fSpinor[1]);
   return *this;
}

inline TPauliSpinor &TPauliSpinor::Normalize()
{
   Normalize(1);
   return *this;
}

inline TPauliSpinor &TPauliSpinor::Normalize(const LDouble_t &norm)
{
   return (*this *= norm/Norm());
}

inline TPauliSpinor &TPauliSpinor::SetPolar(const TUnitVector &pol)
{
   LDouble_t r=0, theta=0, phi=0;
   pol.GetPolar(r,theta,phi);
   SetPolar(theta,phi);
   return *this;
}

inline Complex_t TPauliSpinor::InnerProd(const TPauliSpinor &other)
{
   return ( conj(fSpinor[0])*other.fSpinor[0] + 
            conj(fSpinor[1])*other.fSpinor[1] );
}
 
inline Complex_t TPauliSpinor::ScalarProd(const TPauliSpinor &other)
{
   return InnerProd(other);
}
 
inline TPauliSpinor TPauliSpinor::operator-() const
{
   TPauliSpinor result;
   result.fSpinor[0] = -fSpinor[0];
   result.fSpinor[1] = -fSpinor[1];
   return result;
}

inline TPauliSpinor operator+
                    (const TPauliSpinor &v1, const TPauliSpinor &v2)
{
   TPauliSpinor result(v1);
   result.fSpinor[0] += v2.fSpinor[0];
   result.fSpinor[1] += v2.fSpinor[1];
   return result;
}

inline TPauliSpinor operator+
                    (const TPauliSpinor &v1, const Complex_t *a2)
{
   TPauliSpinor result(v1);
   result.fSpinor[0] += *(a2++);
   result.fSpinor[1] += *a2;
   return result;
}

inline TPauliSpinor operator+(const Complex_t *a1, const TPauliSpinor &v2)
{
   TPauliSpinor result(v2);
   result.fSpinor[0] += *(a1++);
   result.fSpinor[1] += *a1;
   return result;
}

inline TPauliSpinor operator-(const TPauliSpinor &v1, const TPauliSpinor &v2)
{
   TPauliSpinor result(v1);
   result.fSpinor[0] -= v2.fSpinor[0];
   result.fSpinor[1] -= v2.fSpinor[1];
   return result;
}

inline TPauliSpinor operator-(const TPauliSpinor &v1, const Complex_t *a2)
{
   TPauliSpinor result(v1);
   result.fSpinor[0] -= *(a2++);
   result.fSpinor[1] -= *a2;
   return result;
}

inline TPauliSpinor operator-(const Complex_t *a1, const TPauliSpinor &v2)
{
   TPauliSpinor result(v2);
   result.fSpinor[0] -= *(a1++);
   result.fSpinor[1] -= *a1;
   return result;
}

inline TPauliSpinor operator*(const TPauliSpinor &vec, const LDouble_t &factor)
{
   TPauliSpinor result(vec);
   result.fSpinor[0] *= factor;
   result.fSpinor[1] *= factor;
   return result;
}

inline TPauliSpinor operator*(const LDouble_t &factor, const TPauliSpinor &vec)
{
   TPauliSpinor result(vec);
   result.fSpinor[0] *= factor;
   result.fSpinor[1] *= factor;
   return result;
}

inline TPauliSpinor operator*(const TPauliSpinor &vec, const Complex_t &factor)
{
   TPauliSpinor result(vec);
   result.fSpinor[0] *= factor;
   result.fSpinor[1] *= factor;
   return result;
}

inline TPauliSpinor operator*(const Complex_t &factor, const TPauliSpinor &vec)
{
   TPauliSpinor result(vec);
   result.fSpinor[0] *= factor;
   result.fSpinor[1] *= factor;
   return result;
}

inline TPauliSpinor operator/(const TPauliSpinor &vec, const Complex_t &factor)
{
   TPauliSpinor result(vec);
   result.fSpinor[0] /= factor;
   result.fSpinor[1] /= factor;
   return result;
}

inline TBuffer &operator>>(TBuffer &buf, TPauliSpinor *&obj)
{
   for (Int_t i=0; i<2; i++) {
      Double_t real,imag;
      buf >> real >> imag;
      obj->fSpinor[i] = Complex_t(real,imag);
   }
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TPauliSpinor *obj)
{
   for (Int_t i=0; i<2; i++) {
      Double_t real = obj->fSpinor[i].real();
      Double_t imag = obj->fSpinor[i].imag();
      buf << real << imag;
   }
   return buf;
}

#endif
