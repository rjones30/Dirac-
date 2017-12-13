//
// TCrossSection.cxx
//
// author:  Richard T. Jones  11/16/98
// version:  Sep. 10, 2011  v1.01
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
// The TCrossSection class is simply a collection of static functions
// that return the differential cross sections for a small set of
// electromagnetic reactions.  They are calculated to first order in
// QED according to the following Feynman rules.
//
// 1. Cross section is defined as a transition rate density divided by
//    incident flux density.  This assumes a 2-body initial state.
//
// 2. The cross section is invariant with respect to boosts along the
//    axis of relative motion in the initial state.  Furthermore, the
//    states are normalized such that the initial flux density, final
//    density and matrix element factors are each individually Lorentz
//    scalars under boosts along the beam-target axis.
//
//                                            2
//                                    | M   |
//                                4          fi
//            d[sigma]  =  ( 2pi )   --------- d[rho(final)]
//                                     F(in)
//
// 3. F(in) is equal to the product [ 4 E(beam) E(target) ] times the
//    relative velocity between beam and target in the initial state.
//
// 4. rho(final) is calculated in whatever frame the user has chosen
//    to specify the kinematics.  It consists of a term of the form
//                 -1      -3    3
//             (2E)   (2pi)   d p
//    for each final-state fermion or photon, accompanied by a four-
//    dimensional delta function expressing momentum conservation.
//
// 5. M(fi) is calculated in the same frame as was used for F(in) and
//    rho(final).  For tree-level diagrams the following rules apply.
//
// 6. Each fermion line begins and ends with an initial-state or final-
//    state particle.  Each one consists of a series of directed line
//    segments with interaction vertices between them.  For each chain,
//    write down the following terms, starting at the end and working
//    back to the beginning:
//
//    * write down the factor Ubar(p,s) if the line terminates in a
//      final-state fermion, or Vbar(p,s) if it ends with an initial-
//      state antifermion, using the appropriate p,s for this state;
//    * at each vertex write down a factor gamma(mu) for the current;
//    * for each intermediate segment include a propagator of the form
//      1/(pSlash - m) where p is obtained by enforcing momentum
//      conservation at all vertices;
//    * write down the final factor U(p,s) if the line begins with an
//      initial-state ferminon, or V(p,s) if it begins with a final-
//      state antifermion, using the appropriate p,s for this state;
//    * include one power of the coupling e for each vertex;
//    * for each external photon line, write down the 4-vector eps(mu)
//      for initial-state photons, or epsStar(mu) for final-state;
//    * for each internal photon line, write down the photon propagator
//      g(mu,nu)/q2 where q2 is the invariant mass of the virtual photon;
//    * sum over the repeated Lorentz indices.
//
// 7. Counting the powers of 2pi, it turns out that at tree level they
//    cancel in M(fi).  To see this consider the following rules:
//
//    * count -4 for every fermion or photon propagator;
//    * count +4 for every vertex except one;
//
// 8. Compute separate amplitudes M(fi) for each possible combination
//    of initial and final state particle helicities. Then sum over
//    all possible combinations weighted by the spin-density matrix
//    elements for each particle. The spin-density matrix for final
//    state particles is properly considered to be an output to a
//    calculation, not an input. To understand how it is being used as
//    an input in this package, use the following rules.  These rules
//    amount to treating the SDM matrices for initial and final state
//    particle on exactly the same footing.
//
//    * if the SDM specified for a final-state particle corresponds
//      to a pure spin state for that particle, either a lepton or
//      a boson, then the cross section for production into that
//      specific polarization for that final state particle is
//      returned;
//    * if the SDM specified for a final-state particle is simply
//      the unit matrix, then the cross section is computed as a
//      sum over all final polarizations for that particle;
//    * if a more general SDM is specified for a final-state particle
//      then the cross section is calculated as if the SDM gives the
//      detection efficiency for the final-state particle as a
//      function of its polarization state, and the final cross 
//      section weighted by the detection efficiency for that
//      particle is returned.
//
// 9. Gather up all powers of e and rewrite them as sqrt(4 pi alpha)
//
// 10. Include however many powers of hbar and c and powers of 10 are
//     needed to get the result into the appropriate units,
//     eg. microbarns/sr or barns/GeV^2/sr/r.
//
//////////////////////////////////////////////////////////////////////////

#define DEBUGGING 1

#include <iostream>
#include "TCrossSection.h"

ClassImp(TCrossSection)

#include "TPhoton.h"
#include "TLepton.h"
#include "TLorentzBoost.h"

const LDouble_t PI_=2*atan2(1.,0.);

inline LDouble_t sqr(LDouble_t x) { return x*x; }
inline Complex_t sqr(Complex_t x) { return x*x; }

LDouble_t TCrossSection::Compton(const TPhoton &gIn, const TLepton &eIn,
                                 const TPhoton &gOut, const TLepton &eOut)
{
   // Calculates the Compton differential cross section for scattering of
   // a photon from a free lepton.  Units are microbarns per steradian
   // in solid angle of the scattered photon, where the solid angle is
   // that of the photon in the frame chosen by the user.

   TPhoton gIncoming(gIn),  *gI=&gIncoming;
   TLepton eIncoming(eIn),  *eI=&eIncoming;
   TPhoton gOutgoing(gOut), *gF=&gOutgoing;
   TLepton eOutgoing(eOut), *eF=&eOutgoing;

/*******************************************
   TThreeVector bhat(.333,.777,-.666);
   TLorentzBoost btest(bhat,0.0);
   gI->SetMom(gI->Mom().Boost(btest));
   eI->SetMom(eI->Mom().Boost(btest));
   gF->SetMom(gF->Mom().Boost(btest));
   eF->SetMom(eF->Mom().Boost(btest));
*******************************************/

   // Obtain the initial,final lepton state vectors
   TDiracSpinor uI[2];
   uI[0].SetStateU(eI->Mom(), +0.5);
   uI[1].SetStateU(eI->Mom(), -0.5);
   TDiracSpinor uF[2];
   uF[0].SetStateU(eF->Mom(), +0.5);
   uF[1].SetStateU(eF->Mom(), -0.5);

   // Assume without checking that initial,final leptons have same mass
   const LDouble_t mLepton = eI->Mass();

   // Obtain the electron propagators for the two diagrams
   TDiracMatrix dm;
   LDouble_t edenom1 = +2 * eI->Mom().ScalarProd(gI->Mom());
   LDouble_t edenom2 = -2 * eI->Mom().ScalarProd(gF->Mom());
   TDiracMatrix ePropagator1(dm.Slash(eI->Mom() + gI->Mom()) + mLepton);
   TDiracMatrix ePropagator2(dm.Slash(eI->Mom() - gF->Mom()) + mLepton);
   ePropagator1 /= edenom1;
   ePropagator2 /= edenom2;

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   Complex_t invAmp[2][2][2][2];
   for (Int_t gi=0; gi < 2; gi++) {
      for (Int_t gf=0; gf < 2; gf++) {
         TDiracMatrix D;
         TDiracMatrix epsI;
         epsI.Slash(gI->Eps(gi+1));
         TDiracMatrix epsF;
         epsF.Slash(gF->EpsStar(gf+1));
         D = epsF * ePropagator1 * epsI + epsI * ePropagator2 * epsF;
         for (Int_t hi=0; hi < 2; hi++) {
            for (Int_t hf=0; hf < 2; hf++) {
               invAmp[hi][hf][gi][gf] = uF[hf].ScalarProd(D * uI[hi]);
            }
         }
      }
   }

   // Average over initial and final spins
   Complex_t ampSquared(0);
   for (Int_t gi=0; gi < 2; gi++) {
    for (Int_t gibar=0; gibar < 2; gibar++) {
     for (Int_t gf=0; gf < 2; gf++) {
      for (Int_t gfbar=0; gfbar < 2; gfbar++) {
       for (Int_t hi=0; hi < 2; hi++) {
        for (Int_t hibar=0; hibar < 2; hibar++) {
         for (Int_t hf=0; hf < 2; hf++) {
          for (Int_t hfbar=0; hfbar < 2; hfbar++) {
             ampSquared += invAmp[hi][hf][gi][gf] *
                 std::conj(invAmp[hibar][hfbar][gibar][gfbar]) *
                           eI->SDM()[hi][hibar] * 
                           eF->SDM()[hfbar][hf] *
                           gI->SDM()[gi][gibar] *
                           gF->SDM()[gfbar][gf];
          }
         }
        }
       }
      }
     }
    }
   }

   // Obtain the kinematical factors:
   //    (1) 1/flux factor from initial state 1/(4*qin*rootS)
   //    (2) rho from density of final states factor
   // where the general relativistic expression for rho is
   //  rho = pow(2*PI_,4-3*N) delta4(Pin-Pout) [d4 P1] [d4 P2] ... [d4 PN]
   // using differential forms [d4 P] = d4P delta(P.P - m*m) where P.P is
   // the invariant norm of four-vector P, m is the known mass of the
   // corresponding particle, and N is the number of final state particles.
   //    (3) absorb two powers of 4*PI_ into sqr(alphaQED)

   const LDouble_t fluxIn = 4*gI->Mom()[0]*(eI->Mom().Length()+eI->Mom()[0]);
   const LDouble_t rhoFin = sqr(gF->Mom()[0])/eF->Mom().ScalarProd(gF->Mom())/4;
   const LDouble_t kinFactor = 4*rhoFin/fluxIn;

   LDouble_t diffXsect = hbarcSqr*sqr(alphaQED)*real(ampSquared)*kinFactor;

   return diffXsect;

   // The unpolarized Klein Nishina formula is here for comparison
   LDouble_t sinSqrTheta = 1 - sqr(gF->Mom()[3]/gF->Mom()[0]);
   LDouble_t KleinNishinaResult = sqr(alphaQED/mLepton)/2;
   KleinNishinaResult *= sqr(gF->Mom()[0]/gI->Mom()[0]);
   KleinNishinaResult *= (gF->Mom()[0]/gI->Mom()[0]) +
                         (gI->Mom()[0]/gF->Mom()[0]) - sinSqrTheta;
   KleinNishinaResult *= hbarcSqr;
   return KleinNishinaResult;
}

LDouble_t TCrossSection::Bremsstrahlung(const TLepton &eIn,
                                        const TLepton &eOut,
                                        const TPhoton &gOut)
{
   // Calculates the bremsstrahlung cross section for scattering of
   // a lepton from an atom at a particular recoil momentum vector q.
   // The cross section is returned as d(sigma)/(dk dphi d^3 q) where k is
   // the energy of the bremsstrahlung photon and phi is the azimuthal angle
   // of the photon.  The polar angle of the photon is fixed by kinematics.
   // It is assumed that eIn.Mom()[0] = eOut.Mom()[0]+gOut.Mom()[0], that
   // the energy carried away by the recoil is zero in the laboratory frame,
   // but it is not checked.  The calculation is performed in the lab frame.
   // This cross section is only a partial result, because it does not
   // include the integral d^3 q over the form factor of the target.  This
   // depends on the crystal structure of the target atom, and so is left to
   // be carried out by more specialized code.  Units are microbarns/GeV^4/r.

   TLepton eIncoming(eIn),  *eI=&eIncoming;
   TPhoton gOutgoing(gOut), *gF=&gOutgoing;
   TLepton eOutgoing(eOut), *eF=&eOutgoing;

   // Obtain the initial,final lepton state vectors
   TDiracSpinor uI[2];
   uI[0].SetStateU(eI->Mom(), +0.5);
   uI[1].SetStateU(eI->Mom(), -0.5);
   TDiracSpinor uF[2];
   uF[0].SetStateU(eF->Mom(), +0.5);
   uF[1].SetStateU(eF->Mom(), -0.5);

   // Assume without checking that initial, final leptons have same mass
   const LDouble_t mLepton = eI->Mass();

   TFourVectorReal qRecoil(eI->Mom() - eF->Mom() - gF->Mom());

   // Obtain the electron propagators for the two diagrams
   TDiracMatrix dm;
   LDouble_t edenom1 = qRecoil.InvariantSqr() - 2 * qRecoil.ScalarProd(eI->Mom());
   LDouble_t edenom2 = qRecoil.InvariantSqr() + 2 * qRecoil.ScalarProd(eF->Mom());
   TDiracMatrix ePropagator1 = dm.Slash(eI->Mom() - qRecoil) + mLepton;
   TDiracMatrix ePropagator2 = dm.Slash(eF->Mom() + qRecoil) + mLepton;
   ePropagator1 /= edenom1;
   ePropagator2 /= edenom2;

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   Complex_t invAmp[2][2][2];
   for (Int_t gf=0; gf < 2; gf++) {
      TDiracMatrix D;
      TDiracMatrix epsF;
      epsF.Slash(gF->EpsStar(gf+1));
      D = epsF * ePropagator1 * gamma0 + gamma0 * ePropagator2 * epsF;
      for (Int_t hi=0; hi < 2; hi++) {
         for (Int_t hf=0; hf < 2; hf++) {
            invAmp[hi][hf][gf] = uF[hf].ScalarProd(D * uI[hi]);
         }
      }
   }

   // Sum over spins
   Complex_t ampSquared(0);
   Complex_t AAbar[2][2] = {0};
   for (Int_t gf=0; gf < 2; gf++) {
    for (Int_t gfbar=0; gfbar < 2; gfbar++) {
     for (Int_t hi=0; hi < 2; hi++) {
      for (Int_t hibar=0; hibar < 2; hibar++) {
       for (Int_t hf=0; hf < 2; hf++) {
        for (Int_t hfbar=0; hfbar < 2; hfbar++) {
           Complex_t term = invAmp[hi][hf][gf] *
                  std::conj(invAmp[hibar][hfbar][gfbar]) *
                            eI->SDM()[hi][hibar] *
                            eF->SDM()[hfbar][hf] *
                            gF->SDM()[gfbar][gf];
           AAbar[gf][gfbar] += term;
           ampSquared += term;
        }
       }
      }
     }
    }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: bad Bremsstrahlung amplitudes:" << std::endl
                << "  These guys should be all real positive:" << std::endl
                << "    ampSquared = " << ampSquared << std::endl
                << "    AAbar[0][0] = " << AAbar[0][0] << std::endl
                << "    AAbar[1][1] = " << AAbar[1][1] << std::endl
                << "  The rest of these should be conjugate pairs:" << std::endl
                << "    AAbar[i][j]: " << AAbar[0][1] << ", "
                                       << AAbar[1][0] <<  std::endl
      ;
   }
#endif

   // Obtain the kinematical factors:
   //    (1) 1/flux factor from initial state 1/(2E)
   //    (2) rho from density of final states factor
   // where the general relativistic expression for rho is
   //  rho = pow(2*PI_,4-3*N) delta4(Pin-Pout) [d4 P1] [d4 P2] ... [d4 PN]
   // using differential forms [d4 P] = d4P delta(P.P - m*m) where P.P is
   // the invariant norm of four-vector P, m is the known mass of the
   // corresponding particle, and N is the number of final state particles.
   //    (3) 1/pow(qRecoil,4) from the virtual photon propagator
   //    (4) absorb three powers of 4*PI_ into pow(alphaQED,3)
   // To get a simple expression for the density of final states,
   // I redefined the solid angle for the outgoing photon around
   // the momentum axis of the final electron+photon, rather than
   // the incoming electron direction.

   LDouble_t kinFactor = 1/sqr(2*PI_*eI->Mom()[0]); // |qRecoil| << E/c
   LDouble_t diffXsect = hbarcSqr*pow(alphaQED,3)*real(ampSquared)
                         *kinFactor/sqr(qRecoil.InvariantSqr());
   return diffXsect;
}

LDouble_t TCrossSection::PairProduction(const TPhoton &gIn,
                                        const TLepton &eOut,
                                        const TLepton &pOut)
{
   // Calculates the e+e- pair production cross section for a
   // gamma ray off an atom at a particular recoil momentum vector q.
   // The cross section is returned as d(sigma)/(dE dphi d^3q) where E is
   // the energy of the final-state electron and phi is its azimuthal angle.
   // The polar angles of the pair are fixed by momentum conservation.
   // It is assumed that gIn.Mom()[0] = eOut.Mom()[0]+pOut.Mom()[0], that
   // the energy carried away by the recoil is zero in the laboratory frame,
   // but it is not checked.  The calculation is performed in the lab frame.
   // This cross section is only a partial result, because it does not
   // include the integral d^3 q over the form factor of the target.  This
   // depends on the crystal structure of the target atom, and so is left to
   // be carried out by more specialized code.  Units are microbarns/GeV^4/r.

   TPhoton gIncoming(gIn),  *gI=&gIncoming;
   TLepton eOutgoing(eOut), *eF=&eOutgoing;
   TLepton pOutgoing(pOut), *pF=&pOutgoing;

   // Obtain the two lepton state vectors
   TDiracSpinor uF[2];
   uF[0].SetStateU(eF->Mom(), +0.5);
   uF[1].SetStateU(eF->Mom(), -0.5);
   TDiracSpinor vF[2];
   vF[0].SetStateV(pF->Mom(), +0.5);
   vF[1].SetStateV(pF->Mom(), -0.5);

   // Assume without checking that eF, pF are particle, antiparticle
   const LDouble_t mLepton = eF->Mass();

   TFourVectorReal qRecoil(gI->Mom() - eF->Mom() - pF->Mom());

   // Obtain the electron propagators for the two diagrams
   TDiracMatrix dm;
   LDouble_t edenom1 = -2 * gI->Mom().ScalarProd(eF->Mom());
   LDouble_t edenom2 = -2 * gI->Mom().ScalarProd(pF->Mom());
   TDiracMatrix ePropagator1 = dm.Slash(eF->Mom() - gI->Mom()) + mLepton;
   TDiracMatrix ePropagator2 = dm.Slash(gI->Mom() - pF->Mom()) + mLepton;
   ePropagator1 /= edenom1;
   ePropagator2 /= edenom2;

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   Complex_t invAmp[2][2][2];
   for (Int_t gi=0; gi < 2; gi++) {
      TDiracMatrix D;
      TDiracMatrix epsI;
      epsI.Slash(gI->Eps(gi+1));
      D = epsI * ePropagator1 * gamma0 + gamma0 * ePropagator2 * epsI;
      for (Int_t hi=0; hi < 2; hi++) {
         for (Int_t hf=0; hf < 2; hf++) {
            invAmp[hi][hf][gi] = uF[hf].ScalarProd(D * vF[hi]);
         }
      }
   }

   // Sum over spins
   Complex_t ampSquared=0;
   for (Int_t gi=0; gi < 2; gi++) {
    for (Int_t gibar=0; gibar < 2; gibar++) {
     for (Int_t hi=0; hi < 2; hi++) {
      for (Int_t hibar=0; hibar < 2; hibar++) {
       for (Int_t hf=0; hf < 2; hf++) {
        for (Int_t hfbar=0; hfbar < 2; hfbar++) {
           ampSquared += invAmp[hi][hf][gi] *
               std::conj(invAmp[hibar][hfbar][gibar]) *
                         pF->SDM()[hi][hibar] *
                         eF->SDM()[hfbar][hf] *
                         gI->SDM()[gi][gibar];
        }
       }
      }
     }
    }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: bad PairProduction amplitudes:" << std::endl
                << "  These guys should be all real positive:" << std::endl
                << "    ampSquared = " << ampSquared << std::endl;
   }
#endif

   // Obtain the kinematical factors:
   //    (1) 1/flux factor from initial state 1/(2E)
   //    (2) rho from density of final states factor
   // where the general relativistic expression for rho is
   //  rho = pow(2*PI_,4-3*N) delta4(Pin-Pout) [d4 P1] [d4 P2] ... [d4 PN]
   // using differential forms [d4 P] = d4P delta(P.P - m*m) where P.P is
   // the invariant norm of four-vector P, m is the known mass of the
   // corresponding particle, and N is the number of final state particles.
   //    (3) 1/pow(qRecoil,4) from the virtual photon propagator
   //    (4) absorb three powers of 4*PI_ into pow(alphaQED,3)
   // To get a simple expression for the density of final states,
   // I redefined the solid angle for the outgoing electron around
   // the momentum axis of the pair, rather than the incoming photon.

   LDouble_t kinFactor = 1/sqr(2*PI_*gI->Mom()[0]);
   LDouble_t diffXsect = hbarcSqr*pow(alphaQED,3)*real(ampSquared)
                       *kinFactor/sqr(qRecoil.InvariantSqr());
   return diffXsect;
}

LDouble_t TCrossSection::TripletProduction(const TPhoton &gIn,
                                           const TLepton &eIn,
                                           const TLepton &pOut,
                                           const TLepton &eOut2,
                                           const TLepton &eOut3)
{
   // Calculates the e-e+e- triplet production cross section for a gamma
   // ray off a free electron at a particular recoil momentum vector qR.
   // The cross section is returned as d(sigma)/(dE+ dphi+ d^3q) where E+ is
   // the energy of the final-state positron and phi+ is its azimuthal angle
   // about the direction formed by the momentum pOut.Mom() + eOut2.Mom().
   // The polar angles of the pair are fixed by momentum conservation.
   // It is assumed that momentum conservation is respected by the momenta
   //     gIn.Mom() + eIn.Mom() = pOut.Mom() + eOut2.Mom() + eOut3.Mom()
   // but it is not checked.  The calculation is performed in whatever frame
   // the user specifies through the momenta passed in the argument objects.
   // This cross section is only a partial result, because it does not
   // include the integral d^3 q over the form factor of the target.  This
   // depends on the internal structure of the target atom, and so is left to
   // be carried out by more specialized code.  Units are microbarns/GeV^4/r.

   TPhoton gIncoming(gIn), *g0=&gIncoming;
   TLepton eIncoming(eIn), *e0=&eIncoming;
   TLepton pOutgoing(pOut), *e1=&pOutgoing;
   TLepton eOutgoing2(eOut2), *e2=&eOutgoing2;
   TLepton eOutgoing3(eOut3), *e3=&eOutgoing3;

   // Assume without checking that all leptons have the same mass;
   const LDouble_t mLepton = e0->Mass();

   // Obtain the four lepton state vectors
   TDiracSpinor u0[2]; // incoming electron
   u0[0].SetStateU(e0->Mom(), +0.5);
   u0[1].SetStateU(e0->Mom(), -0.5);
   TDiracSpinor v1[2]; // outgoing positron
   v1[0].SetStateV(e1->Mom(), +0.5);
   v1[1].SetStateV(e1->Mom(), -0.5);
   TDiracSpinor u2[2]; // outgoing electron 2
   u2[0].SetStateU(e2->Mom(), +0.5);
   u2[1].SetStateU(e2->Mom(), -0.5);
   TDiracSpinor u3[2]; // outgoing electron 3
   u3[0].SetStateU(e3->Mom(), +0.5);
   u3[1].SetStateU(e3->Mom(), -0.5);

   // There are 8 tree-level diagrams for triplet production.  They can be
   // organized into pairs that share a similar structure.  Two of them
   // resemble Compton scattering with e+e- (Dalitz) splitting of the final
   // gamma (CD), and two resemble gamma decay plus scattering from an electron
   // target (GD).  The next 2 are clones of the CD diagrams, with final-state
   // electrons swapped with each other.  The final 2 are clones of the GD
   // diagrams with final-state electrons swapped.  Each diagram amplitude
   // involves 2 Dirac matrix product chains, one beginning with the final-
   // state positron (1) and the other beginning with the initial-state
   // electron (0).  Each of these comes with one Lorentz index [mu=0..4]
   // and one photon spin index [j=0,1] which must be summed over at the end.
   // The following naming scheme will help to keep track of which amplitude
   // factor is being computed:
   //
   //    {dm}{diag}{swap}
   // where
   //    {dm} = dm or some other symbol for Dirac matrix
   //    {diag} = CD or GD, distinguishes type of diagram
   //    {swap} = 2 or 3, which final electron connects to the initial one
   // For example, dmGD3 refers to the Dirac matrix product coming from the
   // (diag=GD) pair of diagrams with final-state electron (swap=3) connected
   // to the initial-state electron.

   // Pre-compute the electron propagators (a,b suffix for 2 diagrams in pair)
   TDiracMatrix dm;
   LDouble_t edenomCD2a = +2 * g0->Mom().ScalarProd(e0->Mom());
   LDouble_t edenomCD2b = -2 * g0->Mom().ScalarProd(e2->Mom());
   LDouble_t edenomGD2a = -2 * g0->Mom().ScalarProd(e1->Mom());
   LDouble_t edenomGD2b = -2 * g0->Mom().ScalarProd(e3->Mom());
   TDiracMatrix epropCD2a = dm.Slash(g0->Mom() + e0->Mom()) + mLepton;
   TDiracMatrix epropCD2b = dm.Slash(e2->Mom() - g0->Mom()) + mLepton;
   TDiracMatrix epropGD2a = dm.Slash(g0->Mom() - e1->Mom()) + mLepton;
   TDiracMatrix epropGD2b = dm.Slash(e3->Mom() - g0->Mom()) + mLepton;
   epropCD2a /= edenomCD2a;
   epropCD2b /= edenomCD2b;
   epropGD2a /= edenomGD2a;
   epropGD2b /= edenomGD2b;
   TDiracMatrix epropCD3a(epropCD2a);
   TDiracMatrix epropCD3b(epropGD2b);
   TDiracMatrix epropGD3a(epropGD2a);
   TDiracMatrix epropGD3b(epropCD2b);

   // Pre-compute the photon propagators (no a,b suffix needed)
   LDouble_t gpropCD2 = 1 / (e1->Mom() + e3->Mom()).InvariantSqr();
   LDouble_t gpropGD2 = 1 / (e0->Mom() - e2->Mom()).InvariantSqr();
   LDouble_t gpropCD3 = 1 / (e1->Mom() + e2->Mom()).InvariantSqr();
   LDouble_t gpropGD3 = 1 / (e0->Mom() - e3->Mom()).InvariantSqr();

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   const TDiracMatrix gamma1(kDiracGamma1);
   const TDiracMatrix gamma2(kDiracGamma2);
   const TDiracMatrix gamma3(kDiracGamma3);
   TDiracMatrix gamma[4] = {gamma0, gamma1, gamma2, gamma3};

   // Compute the product chains of Dirac matrices
   Complex_t invAmp[2][2][2][2][2] = {0};
   for (Int_t gi=0; gi < 2; gi++) {
      for (Int_t mu=0; mu < 4; mu++) {
         TDiracMatrix epsI;
         epsI.Slash(g0->Eps(gi+1));
         TDiracMatrix CD2;
         CD2 = gamma[mu] * epropCD2a * epsI + epsI * epropCD2b * gamma[mu];
         CD2 *= gpropCD2;
         TDiracMatrix GD2;
         GD2 = gamma[mu] * epropGD2a * epsI + epsI * epropGD2b * gamma[mu];
         GD2 *= gpropGD2;
         TDiracMatrix CD3;
         CD3 = gamma[mu] * epropCD3a * epsI + epsI * epropCD3b * gamma[mu];
         CD3 *= gpropCD3;
         TDiracMatrix GD3;
         GD3 = gamma[mu] * epropGD3a * epsI + epsI * epropGD3b * gamma[mu];
         GD3 *= gpropGD3;
         for (Int_t h0=0; h0 < 2; h0++) {
            for (Int_t h1=0; h1 < 2; h1++) {
               for (Int_t h2=0; h2 < 2; h2++) {
                  for (Int_t h3=0; h3 < 2; h3++) {
                     invAmp[h0][h1][h2][h3][gi] += 
                           Complex_t(((mu == 0)? +1.L : -1.L) * (
                              u3[h3].ScalarProd(gamma[mu] * v1[h1]) *
                              u2[h2].ScalarProd(CD2 * u0[h0])
                            - u2[h2].ScalarProd(gamma[mu] * v1[h1]) *
                              u3[h3].ScalarProd(CD3 * u0[h0])
                            + u2[h2].ScalarProd(gamma[mu] * u0[h0]) *
                              u3[h3].ScalarProd(GD2 * v1[h1])
                            - u3[h3].ScalarProd(gamma[mu] * u0[h0]) *
                              u2[h2].ScalarProd(GD3 * v1[h1]))
                          );
                  }
               }
            }
         }
      }
   }

   // Sum over spins
   Complex_t ampSquared(0);
   for (Int_t gi=0; gi < 2; gi++) {
    for (Int_t gibar=0; gibar < 2; gibar++) {
     for (Int_t h0=0; h0 < 2; ++h0) {
      for (Int_t h0bar=0; h0bar < 2; ++h0bar) {
       for (Int_t h1=0; h1 < 2; ++h1) {
        for (Int_t h1bar=0; h1bar < 2; ++h1bar) {
         for (Int_t h2=0; h2 < 2; ++h2) {
          for (Int_t h2bar=0; h2bar < 2; ++h2bar) {
           for (Int_t h3=0; h3 < 2; ++h3) {
            for (Int_t h3bar=0; h3bar < 2; ++h3bar) {
               ampSquared += invAmp[h0][h1][h2][h3][gi] *
                   std::conj(invAmp[h0bar][h1bar][h2bar][h3bar][gibar]) *
                             e0->SDM()[h0][h0bar] *
                             e1->SDM()[h1][h1bar] *
                             e2->SDM()[h2bar][h2] *
                             e3->SDM()[h3bar][h3] *
                             g0->SDM()[gi][gibar];
            }
           }
          }
         }
        }
       }
      }
     }
    }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: bad triplets amplitude: " << std::endl
                << "  These guys should be all real positive:" << std::endl
                << "    ampSquared = " << ampSquared << std::endl;
   }
#endif

   // Obtain the kinematical factors:
   //    (1) 1/flux from initial state 1/(4 kin [p0 + E0])
   //    (2) rho from density of final states factor
   // where the general relativistic expression for rho is
   //  rho = pow(2*PI_,4-3*N) delta4(Pin-Pout) [d4 P1] [d4 P2] ... [d4 PN]
   // using differential forms [d4 P] = d4P delta(P.P - m*m) where P.P is
   // the invariant norm of four-vector P, m is the known mass of the
   // corresponding particle, and N is the number of final state particles.
   //    (3) absorb three powers of 4*PI_ into pow(alphaQED,3)

   LDouble_t fluxFactor = 4*g0->Mom()[0]*(e0->Mom().Length()+e0->Mom()[0]);
   LDouble_t rhoFactor = 1/(8*e3->Mom()[0]*(e1->Mom()+e2->Mom()).Length());
   LDouble_t piFactor = pow(2*PI_,4-9)*pow(4*PI_,3);
   LDouble_t diffXsect = hbarcSqr * pow(alphaQED,3) * real(ampSquared)
                        / fluxFactor * rhoFactor * piFactor;
   return diffXsect;
}

LDouble_t TCrossSection::BetheHeitlerNucleon(const TPhoton &gIn,
                                             const TLepton &nIn,
                                             const TLepton &pOut,
                                             const TLepton &eOut,
                                             const TLepton &nOut,
                                             LDouble_t F1spacelike,
                                             LDouble_t F2spacelike,
                                             LDouble_t F1timelike,
                                             LDouble_t F2timelike)
{
   // Calculates the e+e- Bethe Heitler production cross section for a gamma
   // ray off a free nucleon at a particular recoil momentum vector qR.
   // This is similar to TripletProduction, except that the nucleon has an
   // amomalous magnetic moment in addition to its Dirac magnetic moment,
   // and there are no identical particles to symmetrize in the final state.
   // The values of the F1 and F2 form factors that accompany these two terms
   // in the current must be evaluated by the user at the value of q2 for the
   // input kinematics and supplied as inputs through arguments F1 and F2.
   // The cross section is returned as d(sigma)/(dE+ dphi+ d^3q) where E+ is
   // the energy of the final-state positron and phi+ is its azimuthal angle
   // about the direction formed by the momentum pOut.Mom() + eOut2.Mom().
   // The polar angles of the pair are fixed by momentum conservation.
   // It is assumed that momentum conservation is respected by the momenta
   //     gIn.Mom() + pIn.Mom() = pOut.Mom() + eOut2.Mom() + eOut3.Mom()
   // but it is not checked.  The calculation is performed in whatever frame
   // the user specifies through the momenta passed in the argument objects.
   // Units are microbarns/GeV^4/r.

   TPhoton gIncoming(gIn), *g0=&gIncoming;
   TLepton nIncoming(nIn), *n0=&nIncoming;
   TLepton pOutgoing(pOut), *e1=&pOutgoing;
   TLepton eOutgoing(eOut), *e2=&eOutgoing;
   TLepton nOutgoing(nOut), *n3=&nOutgoing;

   const LDouble_t mLepton = e1->Mass();
   const LDouble_t mNucleon = n0->Mass();

   // Obtain the four fermion state vectors
   TDiracSpinor u0[2]; // incoming nucleon
   u0[0].SetStateU(n0->Mom(), +0.5);
   u0[1].SetStateU(n0->Mom(), -0.5);
   TDiracSpinor v1[2]; // outgoing positron
   v1[0].SetStateV(e1->Mom(), +0.5);
   v1[1].SetStateV(e1->Mom(), -0.5);
   TDiracSpinor u2[2]; // outgoing electron
   u2[0].SetStateU(e2->Mom(), +0.5);
   u2[1].SetStateU(e2->Mom(), -0.5);
   TDiracSpinor u3[2]; // outgoing nucleon
   u3[0].SetStateU(n3->Mom(), +0.5);
   u3[1].SetStateU(n3->Mom(), -0.5);

   // There are 4 tree-level diagrams for Bethe-Heitler production.  They can
   // be organized into pairs that share a similar structure.  Two of them
   // resemble Compton scattering with e+e- (Dalitz) splitting of the final
   // gamma (CD), and two resemble gamma decay plus rescattering from a nucleon
   // target (GD).  Each diagram amplitude involves 2 Dirac matrix product 
   // chains, one beginning with the final-state positron (1) and the other
   // beginning with the initial-state nucleon (0).  Each of these comes with
   // one Lorentz index [mu=0..4] and one photon spin index [j=0,1] which must
   // be summed over at the end.  The following naming scheme will help to keep
   // track of which amplitude factor is being computed:
   //
   //    {dm}{diag}
   // where
   //    {dm} = dm or some other symbol for Dirac matrix
   //    {diag} = CD or GD, distinguishes type of diagram
   // For example, dmGD refers to the Dirac matrix product coming from the
   // (diag=GD) pair.

   // Pre-compute the fermion propagators (a,b suffix for 2 diagrams in pair)
   TDiracMatrix dm;
   LDouble_t ndenomCDa = +2 * g0->Mom().ScalarProd(n0->Mom());
   LDouble_t ndenomCDb = -2 * g0->Mom().ScalarProd(n3->Mom());
   LDouble_t edenomGDa = -2 * g0->Mom().ScalarProd(e1->Mom());
   LDouble_t edenomGDb = -2 * g0->Mom().ScalarProd(e2->Mom());
   TDiracMatrix npropCDa = dm.Slash(g0->Mom() + n0->Mom()) + mNucleon;
   TDiracMatrix npropCDb = dm.Slash(n3->Mom() - g0->Mom()) + mNucleon;
   TDiracMatrix epropGDa = dm.Slash(g0->Mom() - e1->Mom()) + mLepton;
   TDiracMatrix epropGDb = dm.Slash(e2->Mom() - g0->Mom()) + mLepton;
   npropCDa /= ndenomCDa;
   npropCDb /= ndenomCDb;
   epropGDa /= edenomGDa;
   epropGDb /= edenomGDb;

   // Pre-compute the photon propagators (no a,b suffix needed)
   LDouble_t gpropCD = 1 / (e1->Mom() + e2->Mom()).InvariantSqr();
   LDouble_t gpropGD = 1 / (n0->Mom() - n3->Mom()).InvariantSqr();

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   const TDiracMatrix gamma1(kDiracGamma1);
   const TDiracMatrix gamma2(kDiracGamma2);
   const TDiracMatrix gamma3(kDiracGamma3);
   const TDiracMatrix sigma01(kDiracGamma0, kDiracGamma1);
   const TDiracMatrix sigma02(kDiracGamma0, kDiracGamma2);
   const TDiracMatrix sigma03(kDiracGamma0, kDiracGamma3);
   const TDiracMatrix sigma12(kDiracGamma1, kDiracGamma2);
   const TDiracMatrix sigma13(kDiracGamma1, kDiracGamma3);
   const TDiracMatrix sigma23(kDiracGamma2, kDiracGamma3);

   TDiracMatrix gamma[4] = {gamma0, gamma1, gamma2, gamma3};
   TFourVectorReal qCD(e1->Mom() + e2->Mom());
   TDiracMatrix JnucleonCD[4] =
            {
               gamma0 * F1timelike +
                Complex_t(0, F2timelike / (2 * mNucleon)) *
                 (-sigma01 * qCD[1] - sigma02 * qCD[2] - sigma03 * qCD[3]),
               gamma1 * F1timelike +
                Complex_t(0, F2timelike / (2 * mNucleon)) *
                 (-sigma01 * qCD[0] - sigma12 * qCD[2] - sigma13 * qCD[3]),
               gamma2 * F1timelike +
                Complex_t(0, F2timelike / (2 * mNucleon)) *
                 (-sigma02 * qCD[0] + sigma12 * qCD[1] - sigma23 * qCD[3]),
               gamma3 * F1timelike +
                Complex_t(0, F2timelike / (2 * mNucleon)) *
                 (-sigma03 * qCD[0] + sigma13 * qCD[1] + sigma23 * qCD[2])
            };
   TFourVectorReal qGD(n3->Mom() - n0->Mom());
   TDiracMatrix JnucleonGD[4] = 
            {
               gamma0 * F1spacelike +
                Complex_t(0, F2spacelike / (2 * mNucleon)) *
                 (-sigma01 * qGD[1] - sigma02 * qGD[2] - sigma03 * qGD[3]),
               gamma1 * F1spacelike +
                Complex_t(0, F2spacelike / (2 * mNucleon)) *
                 (-sigma01 * qGD[0] - sigma12 * qGD[2] - sigma13 * qGD[3]),
               gamma2 * F1spacelike +
                Complex_t(0, F2spacelike / (2 * mNucleon)) *
                 (-sigma02 * qGD[0] + sigma12 * qGD[1] - sigma23 * qGD[3]),
               gamma3 * F1spacelike +
                Complex_t(0, F2spacelike / (2 * mNucleon)) *
                 (-sigma03 * qGD[0] + sigma13 * qGD[1] + sigma23 * qGD[2])
            };

   // Compute the product chains of Dirac matrices
   Complex_t invAmp[2][2][2][2][2] = {0};
   for (Int_t gi=0; gi < 2; gi++) {
      for (Int_t mu=0; mu < 4; mu++) {
         TDiracMatrix epsI;
         epsI.Slash(g0->Eps(gi+1));
         TDiracMatrix CD;
         CD = JnucleonCD[mu] * npropCDa * epsI + epsI * npropCDb * JnucleonCD[mu];
         CD *= gpropCD;
         TDiracMatrix GD;
         GD = gamma[mu] * epropGDa * epsI + epsI * epropGDb * gamma[mu];
         GD *= gpropGD;
         for (Int_t h0=0; h0 < 2; h0++) {
            for (Int_t h1=0; h1 < 2; h1++) {
               for (Int_t h2=0; h2 < 2; h2++) {
                  for (Int_t h3=0; h3 < 2; h3++) {
                     invAmp[h0][h1][h2][h3][gi] += 
                           Complex_t(((mu == 0)? +1.L : -1.L) * (
                              u2[h2].ScalarProd(gamma[mu] * v1[h1]) *
                              u3[h3].ScalarProd(CD * u0[h0])
                            + u3[h3].ScalarProd(JnucleonGD[mu] * u0[h0]) *
                              u2[h2].ScalarProd(GD * v1[h1]))
                          );
                  }
               }
            }
         }
      }
   }

   // Sum over spins
   Complex_t ampSquared(0);
   for (Int_t gi=0; gi < 2; gi++) {
    for (Int_t gibar=0; gibar < 2; gibar++) {
     for (Int_t h0=0; h0 < 2; ++h0) {
      for (Int_t h0bar=0; h0bar < 2; ++h0bar) {
       for (Int_t h1=0; h1 < 2; ++h1) {
        for (Int_t h1bar=0; h1bar < 2; ++h1bar) {
         for (Int_t h2=0; h2 < 2; ++h2) {
          for (Int_t h2bar=0; h2bar < 2; ++h2bar) {
           for (Int_t h3=0; h3 < 2; ++h3) {
            for (Int_t h3bar=0; h3bar < 2; ++h3bar) {
               ampSquared += invAmp[h0][h1][h2][h3][gi] *
                   std::conj(invAmp[h0bar][h1bar][h2bar][h3bar][gibar]) *
                             n0->SDM()[h0][h0bar] *
                             e1->SDM()[h1][h1bar] *
                             e2->SDM()[h2bar][h2] *
                             n3->SDM()[h3bar][h3] *
                             g0->SDM()[gi][gibar];
            }
           }
          }
         }
        }
       }
      }
     }
    }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: bad Bethe-Heitler amplitude: " << std::endl
                << "  These guys should be all real positive:" << std::endl
                << "    ampSquared = " << ampSquared << std::endl;
   }
#endif

   // Obtain the kinematical factors:
   //    (1) 1/flux from initial state 1/(4 kin [p0 + E0])
   //    (2) rho from density of final states factor
   // where the general relativistic expression for rho is
   //  rho = pow(2*PI_,4-3*N) delta4(Pin-Pout) [d4 P1] [d4 P2] ... [d4 PN]
   // using differential forms [d4 P] = d4P delta(P.P - m*m) where P.P is
   // the invariant norm of four-vector P, m is the known mass of the
   // corresponding particle, and N is the number of final state particles.
   //    (3) absorb three powers of 4*PI_ into pow(alphaQED,3)

   LDouble_t fluxFactor = 4*g0->Mom()[0]*(n0->Mom().Length()+n0->Mom()[0]);
   LDouble_t rhoFactor = 1/(8*n3->Mom()[0]*(e1->Mom()+e2->Mom()).Length());
   LDouble_t piFactor = pow(2*PI_,4-9)*pow(4*PI_,3);
   LDouble_t diffXsect = hbarcSqr * pow(alphaQED,3) * real(ampSquared)
                        / fluxFactor * rhoFactor * piFactor;
   return diffXsect;
}

LDouble_t TCrossSection::eeBremsstrahlung(const TLepton &eIn0,
                                          const TLepton &eIn1,
                                          const TLepton &eOut2, 
                                          const TLepton &eOut3,
                                          const TPhoton &gOut)
{
   // Calculates the e-,e- bremsstrahlung cross section for the radiative
   // scattering of an energetic electron off a free electron in the target.
   // The cross section is returned as d(sigma)/(dk dphi d^3 q) where k is
   // the energy of the bremsstrahlung photon and phi is the azimuthal angle
   // of the photon.  The polar angle of the photon is fixed by kinematics,
   // given known values for its energy and the recoil momentum vector q.
   // It is assumed that the momenta specified in the input arguments eIn0,
   // eIn1, eOut2, eOut3, and gOut satisfy momentum conservation, that is
   // eIn0.Mom() + eIn1.Mom() == eOut2.Mom() + eOut3.Mom() + gOut.Mom(),
   // but it is not checked.  The calculation is performed in the lab frame.
   // This cross section is only a partial result, because it does not
   // include the integral d^3 q over the form factor of the target.  This
   // depends on the internal structure of the target atom, and so is left to
   // be carried out by more specialized code.  Units are microbarns/GeV^4/r.

   TLepton eIncoming0(eIn0), *e0=&eIncoming0;
   TLepton eIncoming1(eIn1), *e1=&eIncoming1;
   TLepton eOutgoing2(eOut2), *e2=&eOutgoing2;
   TLepton eOutgoing3(eOut3), *e3=&eOutgoing3;
   TPhoton gOutgoing(gOut), *g0=&gOutgoing;

   // The two leptons must be identical, not checked
   const LDouble_t mLepton=e0->Mass();

   // Obtain the four lepton state vectors
   TDiracSpinor u0[2];
   u0[0].SetStateU(e0->Mom(), +0.5);
   u0[1].SetStateU(e0->Mom(), -0.5);
   TDiracSpinor u1[2];
   u1[0].SetStateU(e1->Mom(), +0.5);
   u1[1].SetStateU(e1->Mom(), -0.5);
   TDiracSpinor u2[2];
   u2[0].SetStateU(e2->Mom(), +0.5);
   u2[1].SetStateU(e2->Mom(), -0.5);
   TDiracSpinor u3[2];
   u3[0].SetStateU(e3->Mom(), +0.5);
   u3[1].SetStateU(e3->Mom(), -0.5);

   // There are 8 tree-level diagrams for e,e bremstrahlung.  They can be
   // organized as follows. Diagram A1[A2] has initial[final] state radiation
   // from the lepton leg that connects from eIn0 to eOut2. Diagram B1[B2]
   // has initial[final] state radiation from the lepton leg that connects
   // from eIn1 to eOut3. Diagrams C1[C2] and D1[D2] are copies of A1[A2]
   // and B1[B2] respectively, with eOut2 and eOut3 labels swapped on the
   // outgoing legs. Each diagram amplitude involves 2 Dirac matrix product
   // chains, one connecting to eIn0 and the other to eIn1. Each of these
   // comes with one Lorentz index [mu=0..4] for the photon propagator, and
   // one photon spin index [j=0,1] for the external photon polarization, 
   // which must be summed over at the end.  The following naming scheme
   // will help to keep track of which amplitude factor is being computed:
   //
   //    {dm}{diag}{leg}
   // where
   //    {dm} = dm or some other symbol for Dirac matrix
   //    {diag} = A, B, C, or D to indicate which diagram pair
   //    {leg} = 0 or 1, chain with initial electron 0 or 1
   // For example, dmB1 refers to the Dirac matrix product coming from the
   // leg (leg=1) containing the initial-state electron eIn1, with radiation
   // from the right-hand leg connecting eIn1 to eOut3 (diag=B).

   // Pre-compute the electron propagators
   TDiracMatrix dm;
   LDouble_t edenomA1(-2 * g0->Mom().ScalarProd(e0->Mom()));
   LDouble_t edenomA2(+2 * g0->Mom().ScalarProd(e2->Mom()));
   LDouble_t edenomB1(-2 * g0->Mom().ScalarProd(e1->Mom()));
   LDouble_t edenomB2(+2 * g0->Mom().ScalarProd(e3->Mom()));
   TDiracMatrix epropA1 = dm.Slash(e0->Mom() - g0->Mom()) + mLepton;
   TDiracMatrix epropA2 = dm.Slash(e2->Mom() + g0->Mom()) + mLepton;
   TDiracMatrix epropB1 = dm.Slash(e1->Mom() - g0->Mom()) + mLepton;
   TDiracMatrix epropB2 = dm.Slash(e3->Mom() + g0->Mom()) + mLepton;
   epropA1 /= edenomA1;
   epropA2 /= edenomA2;
   epropB1 /= edenomB1;
   epropB2 /= edenomB2;
   TDiracMatrix epropC1(epropA1);
   TDiracMatrix epropC2(epropB2);
   TDiracMatrix epropD1(epropB1);
   TDiracMatrix epropD2(epropA2);

   // Pre-compute the photon propagators (no 1,2 suffix is needed)
   LDouble_t gpropA = 1 / (e1->Mom() - e3->Mom()).InvariantSqr();
   LDouble_t gpropB = 1 / (e0->Mom() - e2->Mom()).InvariantSqr();
   LDouble_t gpropC = 1 / (e1->Mom() - e2->Mom()).InvariantSqr();
   LDouble_t gpropD = 1 / (e0->Mom() - e3->Mom()).InvariantSqr();

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   const TDiracMatrix gamma1(kDiracGamma1);
   const TDiracMatrix gamma2(kDiracGamma2);
   const TDiracMatrix gamma3(kDiracGamma3);
   TDiracMatrix gamma[4] = {gamma0, gamma1, gamma2, gamma3};

   // Compute the product chains of Dirac matrices
   Complex_t invAmp[2][2][2][2][2] = {0};
   for (Int_t gf=0; gf < 2; gf++) {
      for (Int_t mu=0; mu < 4; mu++) {
         TDiracMatrix epsF;
         epsF.Slash(g0->EpsStar(gf+1));
         TDiracMatrix A;
         A = gamma[mu] * epropA1 * epsF + epsF * epropA2 * gamma[mu];
         A *= gpropA;
         TDiracMatrix B;
         B = gamma[mu] * epropB1 * epsF + epsF * epropB2 * gamma[mu];
         B *= gpropB;
         TDiracMatrix C;
         C = gamma[mu] * epropC1 * epsF + epsF * epropC2 * gamma[mu];
         C *= gpropC;
         TDiracMatrix D;
         D = gamma[mu] * epropD1 * epsF + epsF * epropD2 * gamma[mu];
         D *= gpropD;
         for (Int_t h0=0; h0 < 2; h0++) {
            for (Int_t h1=0; h1 < 2; h1++) {
               for (Int_t h2=0; h2 < 2; h2++) {
                  for (Int_t h3=0; h3 < 2; h3++) {
                     invAmp[h0][h1][h2][h3][gf] +=
                           Complex_t(((mu == 0)? +1.L : -1.L) * (
                              u3[h3].ScalarProd(gamma[mu] * u1[h1]) *
                              u2[h2].ScalarProd(A * u0[h0])
                            + u2[h2].ScalarProd(gamma[mu] * u0[h0]) *
                              u3[h3].ScalarProd(B * u1[h1])
                            - u2[h2].ScalarProd(gamma[mu] * u1[h1]) *
                              u3[h3].ScalarProd(C * u0[h0])
                            - u3[h3].ScalarProd(gamma[mu] * u0[h0]) *
                              u2[h2].ScalarProd(D * u1[h1]))
                          );
                  }
               }
            }
         }
      }
   }

   // Sum over spins
   Complex_t ampSquared(0);
   for (Int_t gf=0; gf < 2; gf++) {
    for (Int_t gfbar=0; gfbar < 2; gfbar++) {
     for (Int_t h0=0; h0 < 2; ++h0) {
      for (Int_t h0bar=0; h0bar < 2; ++h0bar) {
       for (Int_t h1=0; h1 < 2; ++h1) {
        for (Int_t h1bar=0; h1bar < 2; ++h1bar) {
         for (Int_t h2=0; h2 < 2; ++h2) {
          for (Int_t h2bar=0; h2bar < 2; ++h2bar) {
           for (Int_t h3=0; h3 < 2; ++h3) {
            for (Int_t h3bar=0; h3bar < 2; ++h3bar) {
               ampSquared += invAmp[h0][h1][h2][h3][gf] *
                   std::conj(invAmp[h0bar][h1bar][h2bar][h3bar][gfbar]) *
                             e0->SDM()[h0][h0bar] *
                             e1->SDM()[h1][h1bar] *
                             e2->SDM()[h2bar][h2] *
                             e3->SDM()[h3bar][h3] *
                             g0->SDM()[gfbar][gf];
            }
           }
          }
         }
        }
       }
      }
     }
    }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: bad eeBremsstrahlung amplitude: " << std::endl
                << "  These guys should all be real positive: " << std::endl
                << "    ampSquared = " << ampSquared << std::endl;
   }
#endif

   // Obtain the kinematical factors:
   //    (1) 1/flux factor from initial state 1/(4 E0 E1)
   //    (2) rho from density of final states factor
   // where the general relativistic expression for rho is
   //  rho = pow(2*PI_,4-3*N) delta4(Pin-Pout) [d4 P1] [d4 P2] ... [d4 PN]
   // using differential forms [d4 P] = d4P delta(P.P - m*m) where P.P is
   // the invariant norm of four-vector P, m is the known mass of the
   // corresponding particle, and N is the number of final state particles.
   //    (3) absorb three powers of 4*PI_ into pow(alphaQED,3)
   // To get a simple expression for the density of final states,
   // I redefined the solid angle for the outgoing photon around
   // the momentum axis of the final eOut2 + photon, rather than
   // the incoming electron direction.

   LDouble_t kinFactor = 1 / sqr(2*PI_ * e0->Mom()[0]);
   kinFactor /= 4 * e1->Mom()[0] * e3->Mom()[0];
   LDouble_t diffXsect = hbarcSqr * pow(alphaQED,3) *
                         real(ampSquared) * kinFactor;
   return diffXsect;
}

void TCrossSection::Streamer(TBuffer &buf)
{
   // All members are static; this function is a noop.
}

void TCrossSection::Print(Option_t *option)
{
   // All members are static; this function is a noop.
}
