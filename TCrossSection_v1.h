//
// TCrossSection_v1.h 
//
// This file is distributed as part of the Dirac++ package,
// a general toolkit for computing the amplitudes for Feynman
// graphs. See DiracPackage.h for details.



// Warning: the algorithms in this class are now obsolete.
// New class TCrossSection implements the same interface,
// using improved algorithms!
// -- Richard Jones, April 18, 2017



#ifndef ROOT_TCrossSection_v1
#define ROOT_TCrossSection_v1

#include "Double.h"
#include "TBuffer.h"

class TPhoton;
class TLepton;
class TThreeVectorReal;

class TCrossSection_v1 {

public:
   virtual ~TCrossSection_v1() { }

   static LDouble_t Compton(const TPhoton &gIn, const TLepton &eIn,
                            const TPhoton &gOut, const TLepton &eOut);
   static LDouble_t Bremsstrahlung(const TLepton &eIn, const TLepton &eOut,
                                   const TPhoton &gOut);
   static LDouble_t PairProduction(const TPhoton &gIn,
                                   const TLepton &eOut, const TLepton &pOut);
   static LDouble_t TripletProduction(const TPhoton &gIn, const TLepton &eIn,
                                      const TLepton &pOut, const TLepton &eOut2,
                                      const TLepton &eOut3);
   static LDouble_t eeBremsstrahlung(const TLepton &eIn0,
                                     const TLepton &eIn1,
                                     const TLepton &eOut2, 
                                     const TLepton &eOut3,
                                     const TPhoton &gOut);

   void Print(Option_t *option="");

   ClassDef(TCrossSection_v1,1)  // Several useful QED cross sections
};

#endif
