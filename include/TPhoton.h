//
// TPhoton.h 
//
// This file is distributed as part of the Dirac++ package,
// a general toolkit for computing the amplitudes for Feynman
// graphs. See DiracPackage.h for details.

#ifndef ROOT_TPhoton
#define ROOT_TPhoton

#include "Double.h"
#include "TFourVectorComplex.h"
#include "TPauliMatrix.h"

const LDouble_t alphaQED=1/137.036;
const LDouble_t hbarcSqr=389.37966;    // in GeV*GeV*ub


class TPhoton {

protected:
   TFourVectorReal fMomentum;      // four-momentum of photon
   TPauliMatrix    fSpinDensity;   // spin density matrix in helicity basis

public:
   TPhoton() { }
   explicit TPhoton(const TFourVectorReal &p);
   explicit TPhoton(const TThreeVectorReal &p);
   virtual ~TPhoton() { }
   TFourVectorReal Mom() const;
   TThreeVectorReal Pol() const;
   TPauliMatrix &SDM() const;
   TFourVectorComplex Eps(const Int_t mode) const;
   TFourVectorComplex EpsStar(const Int_t mode) const;
   TPhoton &SetMom(const TFourVectorReal &p);
   TPhoton &SetMom(const TThreeVectorReal &p);
   TPhoton &SetPol(const TThreeVectorReal &polar);
   TPhoton &SetPlanePolarization(const TThreeVectorReal &phat,
                                 LDouble_t pol=1);
   TPhoton &SetEllipticalPolarization(const TThreeVectorReal &phat,
                                      LDouble_t circ,
                                      LDouble_t pol=1);
   TThreeVectorReal GetPolarizationPlane() const;
   TPhoton &AllPol();

   friend TBuffer &operator>>(TBuffer &buf, TPhoton *&obj);
   friend TBuffer &operator<<(TBuffer &buf, const TPhoton *obj);
   void Print(Option_t *option="");
 
   ClassDef(TPhoton,1)  // Photon state description class
};

class TGhoston : public TPhoton {

public:
   TGhoston() { }
   TGhoston(const TFourVectorReal &p) : TPhoton(p) { }
   TGhoston(const TThreeVectorReal &p) : TPhoton(p) { }
   TGhoston(const TPhoton &obj) : TPhoton(obj) { }
   ~TGhoston() { }

   TFourVectorComplex Eps(const Int_t mode) const;
   TFourVectorComplex EpsStar(const Int_t mode) const;

   void Print(Option_t *option="");
 
   ClassDef(TGhoston,1)  // Ghost photon for checking gauge invariance
};

inline TPhoton::TPhoton(const TFourVectorReal &p)
{
   fMomentum = p;
   fSpinDensity = TPauliMatrix(0.5);
}

inline TPhoton::TPhoton(const TThreeVectorReal &p)
{
   TFourVectorReal p4(p.Length(),p);
   fMomentum = p4;
   fSpinDensity = TPauliMatrix(0.5);
}

inline TFourVectorReal TPhoton::Mom() const
{
   return fMomentum;
}

inline TPauliMatrix &TPhoton::SDM() const
{
   return (TPauliMatrix &)fSpinDensity;
}

inline TPhoton &TPhoton::SetMom(const TFourVectorReal &p)
{
   fMomentum = p;
   return *this;
}

inline TPhoton &TPhoton::SetMom(const TThreeVectorReal &p)
{
   TFourVectorReal p4(p.Length(),p);
   return SetMom(p4);
}

inline TFourVectorComplex TPhoton::EpsStar(const Int_t mode) const
{
   return Eps(mode).Conj();
}

inline TPhoton &TPhoton::AllPol()
{
   fSpinDensity = TPauliMatrix(1.0);
   return *this;
}

inline TBuffer &operator>>(TBuffer &buf, TPhoton *&obj)
{
   TFourVectorReal *mom=(TFourVectorReal *)&obj->fMomentum;
   TPauliMatrix *sdm=(TPauliMatrix *)&obj->fSpinDensity;
   buf >> mom;
   buf >> sdm;
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const TPhoton *obj)
{
   TFourVectorReal *mom=(TFourVectorReal *)&obj->fMomentum;
   TPauliMatrix *sdm=(TPauliMatrix *)&obj->fSpinDensity;
   buf << mom;
   buf << sdm;
   return buf;
}

inline TFourVectorComplex TGhoston::EpsStar(const Int_t mode) const
{
   return Eps(mode).Conj();
}

#endif
