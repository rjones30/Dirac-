//
// TThreeRotation.h 
//
// This file is distributed as part of the Lorentz++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See LorentzPackage.h for details.

#ifndef ROOT_TThreeRotation
#define ROOT_TThreeRotation
 
#include "Double.h"
#include "TObject.h"
#include "TLorentzTransform.h"
#include "TError.h"
 
#include <math.h>

class TThreeVectorComplex;

 
class TThreeRotation : public TLorentzTransform {
 
public:
   TThreeRotation() : TLorentzTransform() { }
   explicit TThreeRotation(const TThreeVectorReal &axis);
   explicit TThreeRotation(const TUnitVector &ahat, const LDouble_t angle);
   explicit TThreeRotation(const LDouble_t phi, const LDouble_t theta, const LDouble_t psi);
   TThreeRotation(const TThreeRotation &another);
 
   virtual ~TThreeRotation() { }
 
   TThreeVectorReal Axis() const;
   void GetAxis(TUnitVector &ahat, LDouble_t &angle) const;
   void GetEuler(LDouble_t &phi, LDouble_t &theta, LDouble_t &psi) const;

   TThreeRotation &operator=(const TThreeRotation &source);
   TThreeRotation &operator*=(const TThreeRotation &source);
 
   TThreeRotation &Null();
   TThreeRotation &Transpose();
   TThreeRotation &Invert();

   TThreeRotation &SetAxis(const TThreeVectorReal &axis);
   TThreeRotation &SetAxis(const TUnitVector &ahat, const LDouble_t angle);
   TThreeRotation &SetEuler(const LDouble_t &phi,
                       const LDouble_t &theta,
                       const LDouble_t &psi);

   TThreeVectorReal operator*(const TThreeVectorReal &vec) const;
   TThreeVectorComplex operator*(const TThreeVectorComplex &vec) const;
   TThreeRotation operator*(const TThreeRotation &rotOp) const;

   void Print(Option_t *option="");
 
   ClassDef(TThreeRotation,1)  // Space-rotation operator class
};
 
//----- inlines ----------------------------------------------------------------

inline TThreeRotation::TThreeRotation(const TThreeVectorReal &axis)
 : TLorentzTransform()
{
   SetAxis(axis);
}

inline TThreeRotation::TThreeRotation(const TUnitVector &ahat, const LDouble_t angle)
 : TLorentzTransform()
{
   SetAxis(ahat,angle);
}

inline TThreeRotation::TThreeRotation(const LDouble_t phi,
                            const LDouble_t theta,
                            const LDouble_t psi)
 : TLorentzTransform()
{
   SetEuler(phi,theta,psi);
}

inline TThreeRotation::TThreeRotation(const TThreeRotation &another)
 : TLorentzTransform((TLorentzTransform)another)
{}
 
#endif
