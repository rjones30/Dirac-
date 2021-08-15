//
// TCrossSection.h 
//
// This file is distributed as part of the Dirac++ package,
// a general toolkit for computing the amplitudes for Feynman
// graphs. See DiracPackage.h for details.

#ifndef ROOT_TCrossSection
#define ROOT_TCrossSection

#include "Double.h"
#include "TBuffer.h"

class TPhoton;
class TLepton;
class TThreeVectorReal;

class TCrossSection {

public:
   virtual ~TCrossSection() { }

   static LDouble_t Compton(const TPhoton &gIn, const TLepton &eIn,
                            const TPhoton &gOut, const TLepton &eOut);
   static LDouble_t Bremsstrahlung(const TLepton &eIn, const TLepton &eOut,
                                   const TPhoton &gOut);
   static LDouble_t PairProduction(const TPhoton &gIn,
                                   const TLepton &eOut, const TLepton &pOut);
   static LDouble_t TripletProduction(const TPhoton &gIn, const TLepton &eIn,
                                      const TLepton &pOut, const TLepton &eOut2,
                                      const TLepton &eOut3);
   static LDouble_t BetheHeitlerNucleon(const TPhoton &gIn,
                                        const TLepton &nIn,
                                        const TLepton &pOut,
                                        const TLepton &eOut,
                                        const TLepton &nOut,
                                        LDouble_t F1spacelike,
                                        LDouble_t F2spacelike,
                                        LDouble_t F1timelike,
                                        LDouble_t F2timelike);
   static LDouble_t eeBremsstrahlung(const TLepton &eIn0,
                                     const TLepton &eIn1,
                                     const TLepton &eOut2, 
                                     const TLepton &eOut3,
                                     const TPhoton &gOut);
   static LDouble_t ePairProduction(const TLepton &eIn,
                                    const TLepton &eOut,
                                    const TLepton &lpOut,
                                    const TLepton &lnOut);
   static LDouble_t eTripletProduction(const TLepton &eIn,
                                       const TLepton &eOut,
                                       const TLepton &lpOut,
                                       const TLepton &lnOut,
                                       const TLepton &teIn,
                                       const TLepton &teOut);

   void Print(Option_t *option="");

   ClassDef(TCrossSection,1)  // Several useful QED cross sections
};

#endif
