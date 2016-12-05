//
// TLepton.cxx
//
// author:  Richard T. Jones  11/16/98
// version:  Dec. 12, 1998  v1.00
//
/*************************************************************************
 * Copyright(c) 1998, University of Connecticut, All rights reserved.    *
 * Author: Richard T. Jones, Asst. Prof. of Physics                      *
 *                                                                       *
 * Permission to use, copy, modify and distribute this software and its  *
 * documentation for non-commercial purposes is hereby granted without   *
 * fee, provided that the above copyright notice appears in all copies   *
 * and that both the copyright notice and this permission notice appear  *
 * in the supporting documentation. The author makes no claims about the *
 * suitability of this software for any purpose.                         *
 * It is provided "as is" without express or implied warranty.           *
 *************************************************************************/
//////////////////////////////////////////////////////////////////////////
//
// The TLepton class provides a representation for a lepton in terms
// of its four-momentum and polarization.  It can be off-shell.  The
// polarization is represented by its spin density matrix in the rest
// frame.  The axis of quantization is the direction of the momentum,
// of the particle, or the z axis if the particle is at rest.
//
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "TLepton.h"

ClassImp(TLepton)

TThreeVectorReal TLepton::Pol() const
{
   TThreeVectorReal pol;
   Double_t cons;
   fSpinDensity.Decompose(cons,pol);
   return (pol *= 2);
}

void TLepton::Streamer(TBuffer &buf)
{
   // Put/get a TLepton object to/from stream buffer buf.

   TLepton *me=this;
   if (buf.IsReading()) {
      buf >> me;
   } else {
      buf << me;
   }
}

void TLepton::Print(Option_t *option)
{
   // Output a TLepton class member in ascii form.

   std::cout << "TLepton is" << std::endl;
   std::cout << "Mass: " << fMass << std::endl;
   std::cout << "Momentum: (" << fMomentum[0] << ", " << fMomentum[1] << ", "
        << fMomentum[2] << ", " << fMomentum[3] << " )" << std::endl;
   std::cout << "Helicity density matrix:" << std::endl;
   std::cout << fSpinDensity[0][0] << "  " << fSpinDensity[0][1] << std::endl;
   std::cout << fSpinDensity[1][0] << "  " << fSpinDensity[1][1] << std::endl;
}
