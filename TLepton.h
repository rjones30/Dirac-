//
// TLepton.h 
//
// This file is distributed as part of the Dirac++ package,
// a general toolkit for computing the amplitudes for Feynman
// graphs. See DiracPackage.h for details.

#ifndef ROOT_TLepton
#define ROOT_TLepton

#include "TFourVectorComplex.h"
#include "TDiracMatrix.h"

const Double_t mElectron=0.51099907e-3;
const Double_t mMuon=0.105658389;


class TLepton {

protected:
   Double_t        fMass;	  // mass of lepton
   TFourVectorReal fMomentum;	  // four-momentum of lepton
   TPauliMatrix	   fSpinDensity;  // spin density matrix in helicity basis

public:
   TLepton(const Double_t mass=0);
   TLepton(const TFourVectorReal &p, const Double_t mass=0);
   TLepton(const TThreeVectorReal &p, const Double_t mass=0);
   virtual ~TLepton() { }
   Double_t Mass() const;
   TFourVectorReal Mom() const;
   TThreeVectorReal Pol() const;
   TPauliMatrix &SDM() const;
   TLepton &SetMass(Double_t mass);
   TLepton &SetMom(const TFourVectorReal &p);
   TLepton &SetMom(const TThreeVectorReal &p);
   TLepton &SetPol(const TThreeVectorReal &polar);
   TLepton &AllPol();

   friend TBuffer &operator>>(TBuffer &buf, TLepton *&obj);
   friend TBuffer &operator<<(TBuffer &buf, const TLepton *obj);
   void Print(Option_t *option="");
 
   ClassDef(TLepton,1)  // Lepton state description class
};

inline TLepton::TLepton(const Double_t mass)
{
   fMass = mass;
}

inline TLepton::TLepton(const TFourVectorReal &p, const Double_t mass)
{
   fMass = mass;
   fMomentum = p;
   fSpinDensity = 0.5;
}

inline TLepton::TLepton(const TThreeVectorReal &p, const Double_t mass)
{
   TFourVectorReal p4(p.Length(),p);
   p4[0] = sqrt(p4[0]*p4[0] + mass*mass);
   TLepton(p4,mass);
}

inline Double_t TLepton::Mass() const
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

inline TLepton &TLepton::SetMass(Double_t mass)
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
   fSpinDensity = 1.0;
   return *this;
}

inline TBuffer &operator>>(TBuffer &buf, TLepton *&obj)
{
   TFourVectorReal *mom=(TFourVectorReal *)&obj->fMomentum;
   TPauliMatrix *sdm=(TPauliMatrix *)&obj->fSpinDensity;
   buf >> obj->fMass;
   buf >> mom;
   buf >> sdm;
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TLepton *obj)
{
   TFourVectorReal *mom=(TFourVectorReal *)&obj->fMomentum;
   TPauliMatrix *sdm=(TPauliMatrix *)&obj->fSpinDensity;
   buf << obj->fMass;
   buf << mom;
   buf << sdm;
   return buf;
}

#endif
