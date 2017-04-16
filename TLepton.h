//
// TLepton.h 
//
// This file is distributed as part of the Dirac++ package,
// a general toolkit for computing the amplitudes for Feynman
// graphs. See DiracPackage.h for details.

#ifndef ROOT_TLepton
#define ROOT_TLepton

#include "Double.h"
#include "TFourVectorComplex.h"
#include "TDiracMatrix.h"

const LDouble_t mElectron=0.51099907e-3;
const LDouble_t mMuon=0.105658389;


class TLepton {

protected:
   LDouble_t       fMass;         // mass of lepton
   TFourVectorReal fMomentum;     // four-momentum of lepton
   TPauliMatrix    fSpinDensity;  // spin density matrix in helicity basis

public:
   TLepton(const LDouble_t mass=0);
   explicit TLepton(const TFourVectorReal &p, const LDouble_t mass=0);
   explicit TLepton(const TThreeVectorReal &p, const LDouble_t mass=0);
   virtual ~TLepton() { }
   LDouble_t Mass() const;
   TFourVectorReal Mom() const;
   TThreeVectorReal Pol() const;
   TPauliMatrix &SDM() const;
   TLepton &SetMass(LDouble_t mass);
   TLepton &SetMom(const TFourVectorReal &p);
   TLepton &SetMom(const TThreeVectorReal &p);
   TLepton &SetPol(const TThreeVectorReal &polar);
   TLepton &AllPol();

   friend TBuffer &operator>>(TBuffer &buf, TLepton *&obj);
   friend TBuffer &operator<<(TBuffer &buf, const TLepton *obj);
   void Print(Option_t *option="");
 
   ClassDef(TLepton,1)  // Lepton state description class
};

inline TLepton::TLepton(const LDouble_t mass)
{
   fMass = mass;
}

inline TLepton::TLepton(const TFourVectorReal &p, const LDouble_t mass)
{
   fMass = mass;
   fMomentum = p;
   fSpinDensity = TPauliMatrix(0.5);
}

inline TLepton::TLepton(const TThreeVectorReal &p, const LDouble_t mass)
{
   fMass = mass;
   fMomentum = TFourVectorReal(sqrt(p.LengthSqr()+mass*mass),p);
   fSpinDensity = TPauliMatrix(0.5);
}

inline LDouble_t TLepton::Mass() const
{
   return fMass;
}

inline TFourVectorReal TLepton::Mom() const
{
   return fMomentum;
}

inline TPauliMatrix &TLepton::SDM() const
{
   return (TPauliMatrix &)fSpinDensity;
}

inline TLepton &TLepton::SetMass(LDouble_t mass)
{
   fMass = mass;
   return *this;
}

inline TLepton &TLepton::SetMom(const TFourVectorReal &p)
{
   fMomentum = p;
   return *this;
}

inline TLepton &TLepton::SetMom(const TThreeVectorReal &p)
{
   TFourVectorReal p4(p.Length(),p);
   p4[0] = sqrt(p4[0]*p4[0] + fMass*fMass);
   return SetMom(p4);
}

inline TLepton &TLepton::SetPol(const TThreeVectorReal &polar)
{
   fSpinDensity.SetDensity(polar);
   return *this;
}

inline TLepton &TLepton::AllPol()
{
   fSpinDensity = TPauliMatrix(1.0);
   return *this;
}

inline TBuffer &operator>>(TBuffer &buf, TLepton *&obj)
{
   TFourVectorReal *mom=(TFourVectorReal *)&obj->fMomentum;
   TPauliMatrix *sdm=(TPauliMatrix *)&obj->fSpinDensity;
   Double_t mass;
   buf >> mass;
   obj->fMass = mass;
   buf >> mom;
   buf >> sdm;
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TLepton *obj)
{
   TFourVectorReal *mom=(TFourVectorReal *)&obj->fMomentum;
   TPauliMatrix *sdm=(TPauliMatrix *)&obj->fSpinDensity;
   Double_t mass = obj->fMass;
   buf << mass;
   buf << mom;
   buf << sdm;
   return buf;
}

#endif
