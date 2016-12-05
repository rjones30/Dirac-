//
// TThreeRotation.h 
//
// This file is distributed as part of the Lorentz++ package,
// within the Dirac++ toolkit for computing the amplitudes for
// Feynman graphs. See LorentzPackage.h for details.

#ifndef ROOT_TThreeRotation
#define ROOT_TThreeRotation
 
#include "TObject.h"
#include "TLorentzTransform.h"
#include "TError.h"
 
#include <math.h>

class TThreeVectorComplex;

 
class TThreeRotation : public TLorentzTransform {
 
public:
   TThreeRotation() : TLorentzTransform() { }
   TThreeRotation(const TThreeVectorReal &axis);
   TThreeRotation(const TUnitVector &ahat, const Double_t angle);
   TThreeRotation(const Double_t phi, const Double_t theta, const Double_t psi);
   TThreeRotation(const TThreeRotation &another);
 
   virtual ~TThreeRotation() { }
 
   TThreeVectorReal Axis() const;
   void GetAxis(TUnitVector &ahat, Double_t &angle) const;
   void GetEuler(Double_t &phi, Double_t &theta, Double_t &psi) const;

   TThreeRotation &operator=(const TThreeRotation &source);
   TThreeRotation &operator*=(const TThreeRotation &source);
 
   TThreeRotation &Null();
   TThreeRotation &Transpose();
   TThreeRotation &Invert();

   TThreeRotation &SetAxis(const TThreeVectorReal &axis);
   TThreeRotation &SetAxis(const TUnitVector &ahat, const Double_t angle);
   TThreeRotation &SetEuler(const Double_t &phi,
                       const Double_t &theta,
                       const Double_t &psi);

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

inline TThreeRotation::TThreeRotation(const TUnitVector &ahat, const Double_t angle)
 : TLorentzTransform()
{
   SetAxis(ahat,angle);
}

inline TThreeRotation::TThreeRotation(const Double_t phi,
                            const Double_t theta,
                            const Double_t psi)
 : TLorentzTransform()
{
   SetEuler(phi,theta,psi);
}

inline TThreeRotation::TThreeRotation(const TThreeRotation &another)
 : TLorentzTransform((TLorentzTransform)another)
{}
 
#endif
