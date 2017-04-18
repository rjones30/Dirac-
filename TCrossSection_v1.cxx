//
// TCrossSection_v1.cxx
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



// The algorithms for computing QED spin-dependent amplitudes implemented
// in this class turned out to be vulnerable to large rounding errors for
// a broad range of kinematics that are important in scattering processes
// at high energy. This is because they first square the amplitudes, then
// compute the sums over Dirac indices. Contracting over Dirac indices
// before squaring the amplitudes allows large fermion propagators that
// appear when the kinematics come out close to the pole to cancel with
// many orders of magnitude smaller rounding errors. The improved code
// is now found in class TCrossSection, and TCrossSection_v1 is obsolete
// except for purposes of comparison.
// -- Richard Jones, April 18, 2017



//////////////////////////////////////////////////////////////////////////
//
// The TCrossSection_v1 class is simply a collection of static functions
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
//    final-state fermion, or Vbar(p,s) if it ends with an initial-
//    state antifermion, using the appropriate p,s for this state;
//    * at each vertex write down a factor gamma(mu) for the current;
//    * for each intermediate segment include a propagator of the form
//      1/(pSlash - m) where p is obtained by enforcing momentum
//      conservation at all vertices;
//    * write down the final factor U(p,s) if the line begins with an
//      initial-state ferminon, of V(p,s) if it begins with a final-
//      state antifermion, using the appropriate p,s for this state;
//    * include one power of the coupling e for each vertex;
//    * for each external photon line, write down the 4-vector eps(mu)
//      and contract the index mu with the appropriate current index;
//    * for each internal photon line, write down the photon propagator
//      g(mu,nu)/q2 where q2 is the invariant mass of the virtual
//      photon, and contract the mu,nu with the appropriate currents.
//
// 7. Counting the powers of 2pi, it turns out that at tree level they
//    cancel in M(fi).  To see this consider the following rules:
//
//    * count -4 for every fermion or photon propagator;
//    * count +4 for every vertex except one;
//
// 8. Coming to square the factor M(fi), each fermion chain becomes a
//    ring by appending to each chain a second copy in adjoint order,
//    and taking the trace over dangling Dirac indices.  Such rings
//    contain two terms of the form U(p,s)Ubar(p,s) or V(p,s)Vbar(p,s),
//    one for the incoming and one for the outgoing leg of the chain.
//    Factor out the gamma(0) from the Ubar or Vbar, and then write
//    projector matrices P+(p,s) or P-(p,s) in to express the sum over
//    final spins or average over initial.  For example, to sum over
//    spins of a final-state electron, sum[P+(p,s)] = (pSlash + m).
//
// 9. Gather up all powers of e and rewrite them as sqrt(4 pi alpha)
//
// 10. Include however many powers of hbar and c and powers of 10 are
//    needed to get the result into the appropriate units,
//    eg. microbarns/sr or barns/GeV^2/sr/r.
//
//////////////////////////////////////////////////////////////////////////

#define DEBUGGING 1

#include <iostream>
#include "TCrossSection_v1.h"

ClassImp(TCrossSection_v1)

#include "TPhoton.h"
#include "TLepton.h"
#include "TLorentzBoost.h"

const LDouble_t PI_=2*atan2(1.,0.);

inline LDouble_t sqr(LDouble_t x) { return x*x; }
inline Complex_t sqr(Complex_t x) { return x*x; }

LDouble_t TCrossSection_v1::Compton(const TPhoton &gIn, const TLepton &eIn,
                                 const TPhoton &gOut, const TLepton &eOut)
{
   // Calculates the Compton differential cross section for scattering of
   // a photon from a free lepton.  Units are microbarns per steradian
   // in solid angle of the scattered photon, where the solid angle is
   // that of the photon in the frame chosen by the user.

   TPhoton gIncoming(gIn), *gI=&gIncoming;
   TLepton eIncoming(eIn), *eI=&eIncoming;
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

   // Obtain the initial,final lepton state matrices
   TDiracMatrix chiIn, chiOut;
   chiIn.SetUUbar(eI->Mom(),eI->SDM());
   chiOut.SetUUbar(eF->Mom(),eF->SDM());
   const LDouble_t mLepton=eI->Mass();

   // Obtain the electron propagators for the two diagrams
   TDiracMatrix ePropagator1, ePropagator2, dm;
   ePropagator1 = 1/(dm.Slash(eI->Mom() + gI->Mom()) - mLepton);
   ePropagator2 = 1/(dm.Slash(eI->Mom() - gF->Mom()) - mLepton);

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   TDiracMatrix invAmp[2][2];
   TDiracMatrix invAmpBar[2][2];
   for (Int_t i=0; i<2; i++) {
      for (Int_t j=0; j<2; j++) {
         TDiracMatrix current1, current2;
         current1.Slash(gF->EpsStar(j+1));
         current1 *= ePropagator1;
         current1 *= dm.Slash(gI->Eps(i+1));
         current2.Slash(gI->Eps(i+1));
         current2 *= ePropagator2;
         current2 *= dm.Slash(gF->EpsStar(j+1));
         invAmp[i][j] = invAmpBar[i][j] = current1 + current2;
         invAmp[i][j] *= chiIn;
         invAmpBar[i][j].Adjoint();
         invAmpBar[i][j].UniTransform(gamma0);
         invAmpBar[i][j] *= chiOut;
      }
   }

   // Average over initial and final spins
   Complex_t ampSquared=0;
   for (Int_t i=0; i<2; i++) {
      for (Int_t j=0; j<2; j++) {
         for (Int_t ii=0; ii<2; ii++) {
            for (Int_t jj=0; jj<2; jj++) {
               ampSquared += (invAmp[i][j]*invAmpBar[ii][jj]).Trace()
                             * gI->SDM()[i][ii] * gF->SDM()[jj][j];
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
   std::cout << "returning the Klein Nishina value..." << std::endl;
   return KleinNishinaResult;
}

LDouble_t TCrossSection_v1::Bremsstrahlung(const TLepton &eIn,
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

   TLepton eIncoming(eIn), *eI=&eIncoming;
   TPhoton gOutgoing(gOut), *gF=&gOutgoing;
   TLepton eOutgoing(eOut), *eF=&eOutgoing;

   // Obtain the initial,final lepton state matrices
   TDiracMatrix chiIn, chiOut;
   chiIn.SetUUbar(eI->Mom(),eI->SDM());
   chiOut.SetUUbar(eF->Mom(),eF->SDM());
   const LDouble_t mLepton=eI->Mass();
   TFourVectorReal qRecoil(eI->Mom() - eF->Mom() - gF->Mom());

   // Obtain the electron propagators for the two diagrams
   TDiracMatrix ePropagator1, ePropagator2, dm;
   ePropagator1 = (dm.Slash(eI->Mom() - qRecoil) + mLepton) /
                  (qRecoil.InvariantSqr() - 2*qRecoil.ScalarProd(eI->Mom()));
   ePropagator2 = (dm.Slash(eF->Mom() + qRecoil) + mLepton) /
                  (qRecoil.InvariantSqr() + 2*qRecoil.ScalarProd(eF->Mom()));

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   TDiracMatrix invAmp[2];
   TDiracMatrix invAmpBar[2];
   for (Int_t j=0; j<2; j++) {
      TDiracMatrix current1, current2;
      current1.Slash(gF->EpsStar(j+1));
      current1 *= ePropagator1;
      current1 *= gamma0;
      current2 = gamma0;
      current2 *= ePropagator2;
      current2 *= dm.Slash(gF->EpsStar(j+1));
      invAmp[j] = invAmpBar[j] = current1 + current2;
      invAmp[j] *= chiIn;
      invAmpBar[j].Adjoint();
      invAmpBar[j].UniTransform(gamma0);
      invAmpBar[j] *= chiOut;
   }

   // Sum over spins
   Complex_t ampSquared=0;
   for (Int_t j=0; j<2; j++) {
      for (Int_t jj=0; jj<2; jj++) {
         ampSquared += (invAmp[j]*invAmpBar[jj]).Trace() * gF->SDM()[jj][j];
      }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: problem with Bremsstrahlung amplitudes:" << std::endl
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
   // I redefined the solid angle for the outgoing photon around
   // the momentum axis of the final electron+photon, rather than
   // the incoming electron direction.

   LDouble_t kinFactor = 1/sqr(2*PI_*eI->Mom()[0]); // |qRecoil| << E/c
   LDouble_t diffXsect = hbarcSqr*pow(alphaQED,3)*real(ampSquared)
                         *kinFactor/sqr(qRecoil.InvariantSqr());
   return diffXsect;
}

LDouble_t TCrossSection_v1::PairProduction(const TPhoton &gIn,
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

   TPhoton gIncoming(gIn), *gI=&gIncoming;
   TLepton eOutgoing(eOut), *eF=&eOutgoing;
   TLepton pOutgoing(pOut), *pF=&pOutgoing;

   // Obtain the two lepton state matrices
   TDiracMatrix chiEle, chiPos;
   chiEle.SetUUbar(eF->Mom(),eF->SDM());
   chiPos.SetVVbar(pF->Mom(),pF->SDM());
   const LDouble_t mLepton=eF->Mass();
   TFourVectorReal qRecoil(gI->Mom() - eF->Mom() - pF->Mom());

   // Obtain the electron propagators for the two diagrams
   TDiracMatrix ePropagator1, ePropagator2, dm;
   ePropagator1 = (dm.Slash(eF->Mom() - gI->Mom()) + mLepton) /
                  (-2*gI->Mom().ScalarProd(eF->Mom()));
   ePropagator2 = (dm.Slash(gI->Mom() - pF->Mom()) + mLepton) /
                  (-2*gI->Mom().ScalarProd(pF->Mom()));

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   TDiracMatrix invAmp[2];
   TDiracMatrix invAmpBar[2];
   for (Int_t j=0; j<2; j++) {
      TDiracMatrix current1, current2;
      current1.Slash(gI->Eps(j+1));
      current1 *= ePropagator1;
      current1 *= gamma0;
      current2 = gamma0;
      current2 *= ePropagator2;
      current2 *= dm.Slash(gI->Eps(j+1));
      invAmp[j] = invAmpBar[j] = current1 + current2;
      invAmp[j] *= chiPos;
      invAmpBar[j].Adjoint();
      invAmpBar[j].UniTransform(gamma0);
      invAmpBar[j] *= chiEle;
   }

   // Sum over spins
   Complex_t ampSquared=0;
   for (Int_t j=0; j<2; j++) {
      for (Int_t jj=0; jj<2; jj++) {
         ampSquared += (invAmp[j]*invAmpBar[jj]).Trace() * gI->SDM()[j][jj];
      }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: problem with PairProduction amplitudes:" << std::endl
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

LDouble_t TCrossSection_v1::TripletProduction(
                        const TPhoton &gIn, const TLepton &eIn,
                        const TLepton &pOut, const TLepton &eOut2,
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

   const LDouble_t mLepton=e0->Mass();

   // Obtain the four lepton state matrices
   TDiracMatrix chi0,chi1,chi2,chi3;
   chi0.SetUUbar(e0->Mom(),e0->SDM());
   chi1.SetVVbar(e1->Mom(),e1->SDM());
   chi2.SetUUbar(e2->Mom(),e2->SDM());
   chi3.SetUUbar(e3->Mom(),e3->SDM());

   // There are 8 tree-level diagrams for triplet production.  They can be
   // organized into pairs that share a similar structure.  Two of them
   // resemble Compton scattering with e+e- (Dalitz) splitting of the final
   // gamma (CD), and two resemble Bethe-Heitler scattering from an electron
   // target (BH).  The next 2 are clones of the CD diagrams, with final-state
   // electrons swapped with each other.  The final 2 are clones of the BH
   // diagrams with final-state electrons swapped.  Each diagram amplitude
   // involves 2 Dirac matrix product chains, one beginning with the final-
   // state positron (1) and the other beginning with the initial-state
   // electron (0).  Each of these comes with one Lorentz index [mu=0..4]
   // and one photon spin index [j=0,1] which must be summed over at the end.
   // The following naming scheme will help to keep track of which amplitude
   // factor is being computed:
   //
   //    {dm}{diag}{swap}{leg}{Bar}[mu][j]
   // where
   //    {dm} = dm or some other symbol for Dirac matrix
   //    {diag} = CD or BH, distinguishes type of diagram
   //    {swap} = 2 or 3, which final electron connects to the initial one
   //    {leg} = 0 or 1, chain with initial electron (0) or final positron (1)
   //    {Bar} = "" or "Bar", representing the matrix or its adjoint partner
   //    [mu] = Lorentz index for amplitude factor
   //    [j] = initial photon spin index for amplitude factor
   // For example, dmBH31Bar[3][0] refers to the adjoint pair of the Dirac
   // matrix product (Bar) coming from the leg (leg=1) containing the final-
   // state positron of the Bethe-Heitler (diag=BH) pair of diagrams with
   // final-state electron (swap=3) connected to the initial-state electron,
   // Lorentz component (mu=3), photon spin component (j=0).
   //
   // The plan for computing the sums is as follows:
   //   1. compute all Dirac matrix chains for each leg of each diagram
   //   2. compute the adjoint pair for each of the above (xxxBar matrices)
   //   3. append chi matrices for the u(p) or v(p) spinor factors to each
   //   4. take traces of chains of two xxx[mu][j] and two xxxBar[nu][jj]
   //   5. sum above traces over diag,swap,diagBar,swapBar,mu,nu
   //   6. contract sum[j][jj] with the photon spin density matrix
   // The result of the final contraction is |Mfi|^2

   // Pre-compute the electron propagators (a,b suffix for 2 diagrams in pair)
   TDiracMatrix dm;
   TDiracMatrix epropCD2a = (dm.Slash(g0->Mom() + e0->Mom()) + mLepton) /
                            (2*g0->Mom().ScalarProd(e0->Mom()));
   TDiracMatrix epropCD2b = (dm.Slash(e2->Mom() - g0->Mom()) + mLepton) /
                            (-2*g0->Mom().ScalarProd(e2->Mom()));
   TDiracMatrix epropBH2a = (dm.Slash(g0->Mom() - e1->Mom()) + mLepton) /
                            (-2*g0->Mom().ScalarProd(e1->Mom()));
   TDiracMatrix epropBH2b = (dm.Slash(e3->Mom() - g0->Mom()) + mLepton) /
                            (-2*g0->Mom().ScalarProd(e3->Mom()));
   TDiracMatrix epropCD3a(epropCD2a);
   TDiracMatrix epropCD3b(epropBH2b);
   TDiracMatrix epropBH3a(epropBH2a);
   TDiracMatrix epropBH3b(epropCD2b);

   // Pre-compute the photon propagators (no a,b suffix needed)
   LDouble_t gpropCD2 = 1/(e1->Mom()+e3->Mom()).InvariantSqr();
   LDouble_t gpropBH2 = 1/(e0->Mom()-e2->Mom()).InvariantSqr();
   LDouble_t gpropCD3 = 1/(e1->Mom()+e2->Mom()).InvariantSqr();
   LDouble_t gpropBH3 = 1/(e0->Mom()-e3->Mom()).InvariantSqr();

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   const TDiracMatrix gamma1(kDiracGamma1);
   const TDiracMatrix gamma2(kDiracGamma2);
   const TDiracMatrix gamma3(kDiracGamma3);
   TDiracMatrix gamma[4];
   gamma[0] = gamma0;
   gamma[1] = gamma1;
   gamma[2] = gamma2;
   gamma[3] = gamma3;

   // Compute the product chains of Dirac matrices
   TDiracMatrix dmCD20[4][2],dmCD20Bar[4][2];
   TDiracMatrix dmCD21[4][2],dmCD21Bar[4][2];
   TDiracMatrix dmBH20[4][2],dmBH20Bar[4][2];
   TDiracMatrix dmBH21[4][2],dmBH21Bar[4][2];
   TDiracMatrix dmCD30[4][2],dmCD30Bar[4][2];
   TDiracMatrix dmCD31[4][2],dmCD31Bar[4][2];
   TDiracMatrix dmBH30[4][2],dmBH30Bar[4][2];
   TDiracMatrix dmBH31[4][2],dmBH31Bar[4][2];
   for (Int_t j=0; j<2; j++) {
      for (Int_t mu=0; mu<4; mu++) {
         TDiracMatrix currenta, currentb;

         currenta = gamma[mu];
         currenta *= epropCD2a;
         currenta *= dm.Slash(g0->Eps(j+1));
         currentb.Slash(g0->Eps(j+1));
         currentb *= epropCD2b;
         currentb *= gamma[mu];
         dmCD20[mu][j] = currenta + currentb;
         dmCD21[mu][j] = gamma[mu];

         currenta = gamma[mu];
         currenta *= epropBH2a;
         currenta *= dm.Slash(g0->Eps(j+1));
         currentb.Slash(g0->Eps(j+1));
         currentb *= epropBH2b;
         currentb *= gamma[mu];
         dmBH21[mu][j] = currenta + currentb;
         dmBH20[mu][j] = gamma[mu];

         currenta = gamma[mu];
         currenta *= epropCD3a;
         currenta *= dm.Slash(g0->Eps(j+1));
         currentb.Slash(g0->Eps(j+1));
         currentb *= epropCD3b;
         currentb *= gamma[mu];
         dmCD30[mu][j] = currenta + currentb;
         dmCD31[mu][j] = gamma[mu];

         currenta = gamma[mu];
         currenta *= epropBH3a;
         currenta *= dm.Slash(g0->Eps(j+1));
         currentb.Slash(g0->Eps(j+1));
         currentb *= epropBH3b;
         currentb *= gamma[mu];
         dmBH31[mu][j] = currenta + currentb;
         dmBH30[mu][j] = gamma[mu];
      }
   }

   // Compute adjoint pairs and append chi matrices
   for (Int_t j=0; j<2; j++) {
      for (Int_t mu=0; mu<4; mu++) {
         dmCD20Bar[mu][j] = dmCD20[mu][j];
         dmCD20Bar[mu][j].Adjoint();
         dmCD20Bar[mu][j].UniTransform(gamma0);
         dmCD20Bar[mu][j] *= chi2;
         dmCD20[mu][j] *= chi0;
         dmCD21Bar[mu][j] = dmCD21[mu][j];
         // dmCD21Bar[mu][j].Adjoint();
         // dmCD21Bar[mu][j].UniTransform(gamma0);
         dmCD21Bar[mu][j] *= chi3;
         dmCD21[mu][j] *= chi1;
 
         dmBH20Bar[mu][j] = dmBH20[mu][j];
         // dmBH20Bar[mu][j].Adjoint();
         // dmBH20Bar[mu][j].UniTransform(gamma0);
         dmBH20Bar[mu][j] *= chi2;
         dmBH20[mu][j] *= chi0;
         dmBH21Bar[mu][j] = dmBH21[mu][j];
         dmBH21Bar[mu][j].Adjoint();
         dmBH21Bar[mu][j].UniTransform(gamma0);
         dmBH21Bar[mu][j] *= chi3;
         dmBH21[mu][j] *= chi1;
 
         dmCD30Bar[mu][j] = dmCD30[mu][j];
         dmCD30Bar[mu][j].Adjoint();
         dmCD30Bar[mu][j].UniTransform(gamma0);
         dmCD30Bar[mu][j] *= chi3;
         dmCD30[mu][j] *= chi0;
         dmCD31Bar[mu][j] = dmCD31[mu][j];
         // dmCD31Bar[mu][j].Adjoint();
         // dmCD31Bar[mu][j].UniTransform(gamma0);
         dmCD31Bar[mu][j] *= chi2;
         dmCD31[mu][j] *= chi1;
 
         dmBH30Bar[mu][j] = dmBH30[mu][j];
         // dmBH30Bar[mu][j].Adjoint();
         // dmBH30Bar[mu][j].UniTransform(gamma0);
         dmBH30Bar[mu][j] *= chi3;
         dmBH30[mu][j] *= chi0;
         dmBH31Bar[mu][j] = dmBH31[mu][j];
         dmBH31Bar[mu][j].Adjoint();
         dmBH31Bar[mu][j].UniTransform(gamma0);
         dmBH31Bar[mu][j] *= chi2;
         dmBH31[mu][j] *= chi1;
      }
   }

   // Finally, the sums over traces
   Complex_t Mfi2[2][2]; 
   Complex_t CD2CD2bar[2][2] = {0};
   Complex_t CD2BH2bar[2][2] = {0};
   Complex_t CD2CD3bar[2][2] = {0};
   Complex_t CD2BH3bar[2][2] = {0};
   Complex_t BH2CD2bar[2][2] = {0};
   Complex_t BH2BH2bar[2][2] = {0};
   Complex_t BH2CD3bar[2][2] = {0};
   Complex_t BH2BH3bar[2][2] = {0};
   Complex_t CD3CD2bar[2][2] = {0};
   Complex_t CD3BH2bar[2][2] = {0};
   Complex_t CD3CD3bar[2][2] = {0};
   Complex_t CD3BH3bar[2][2] = {0};
   Complex_t BH3CD2bar[2][2] = {0};
   Complex_t BH3BH2bar[2][2] = {0};
   Complex_t BH3CD3bar[2][2] = {0};
   Complex_t BH3BH3bar[2][2] = {0};
   for (Int_t j=0; j<2; j++) {
      for (Int_t jj=0; jj<2; jj++) {
         for (Int_t mu=0; mu<4; mu++) {
            for (Int_t nu=0; nu<4; nu++) {
               LDouble_t sign;
               sign = (mu == 0)? 1 : -1;
               sign *= (nu == 0)? 1 : -1;
               TDiracMatrix dm0,dm1;
               // CD2 * CD2Bar
               dm0 = dmCD20[mu][j];
               dm0 *= dmCD20Bar[nu][jj];
               dm1 = dmCD21[mu][j];
               dm1 *= dmCD21Bar[nu][jj];
               CD2CD2bar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD2*gpropCD2;
               // CD2 * BH2Bar
               dm0 = dmCD20[mu][j];
               dm0 *= dmBH20Bar[nu][jj];
               dm1 = dmCD21[mu][j];
               dm1 *= dmBH21Bar[nu][jj];
               CD2BH2bar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD2*gpropBH2;
               // CD2 * CD3Bar
               dm0 = dmCD20[mu][j];
               dm0 *= dmCD30Bar[nu][jj];
               dm0 *= dmCD21[mu][j];
               dm0 *= dmCD31Bar[nu][jj];
               CD2CD3bar[j][jj] -= sign*dm0.Trace()*gpropCD2*gpropCD3;
               // CD2 * BH3Bar
               dm0 = dmCD20[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm0 *= dmCD21[mu][j];
               dm0 *= dmBH31Bar[nu][jj];
               CD2BH3bar[j][jj] -= sign*dm0.Trace()*gpropCD2*gpropBH3;
               // BH2 * CD2Bar
               dm0 = dmBH20[mu][j];
               dm0 *= dmCD20Bar[nu][jj];
               dm1 = dmBH21[mu][j];
               dm1 *= dmCD21Bar[nu][jj];
               BH2CD2bar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD2*gpropBH2;
               // BH2 * BH2Bar
               dm0 = dmBH20[mu][j];
               dm0 *= dmBH20Bar[nu][jj];
               dm1 = dmBH21[mu][j];
               dm1 *= dmBH21Bar[nu][jj];
               BH2BH2bar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropBH2*gpropBH2;
               // BH2 * CD3Bar
               dm0 = dmBH20[mu][j];
               dm0 *= dmCD30Bar[nu][jj];
               dm0 *= dmBH21[mu][j];
               dm0 *= dmCD31Bar[nu][jj];
               BH2CD3bar[j][jj] -= sign*dm0.Trace()*gpropBH2*gpropCD3;
               // BH2 * BH3Bar
               dm0 = dmBH20[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm0 *= dmBH21[mu][j];
               dm0 *= dmBH31Bar[nu][jj];
               BH2BH3bar[j][jj] -= sign*dm0.Trace()*gpropBH2*gpropBH3;
               // CD3 * CD2Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmCD20Bar[nu][jj];
               dm0 *= dmCD31[mu][j];
               dm0 *= dmCD21Bar[nu][jj];
               CD3CD2bar[j][jj] -= sign*dm0.Trace()*gpropCD2*gpropCD3;
               // CD3 * BH2Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmBH20Bar[nu][jj];
               dm0 *= dmCD31[mu][j];
               dm0 *= dmBH21Bar[nu][jj];
               CD3BH2bar[j][jj] -= sign*dm0.Trace()*gpropBH2*gpropCD3;
               // CD3 * CD3Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmCD30Bar[nu][jj];
               dm1 = dmCD31[mu][j];
               dm1 *= dmCD31Bar[nu][jj];
               CD3CD3bar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD3*gpropCD3;
               // CD3 * BH3Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm1 = dmCD31[mu][j];
               dm1 *= dmBH31Bar[nu][jj];
               CD3BH3bar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD3*gpropBH3;
               // BH3 * CD2Bar
               dm0 = dmBH30[mu][j];
               dm0 *= dmCD20Bar[nu][jj];
               dm0 *= dmBH31[mu][j];
               dm0 *= dmCD21Bar[nu][jj];
               BH3CD2bar[j][jj] -= sign*dm0.Trace()*gpropBH3*gpropCD2;
               // BH3 * BH2Bar
               dm0 = dmBH30[mu][j];
               dm0 *= dmBH20Bar[nu][jj];
               dm0 *= dmBH31[mu][j];
               dm0 *= dmBH21Bar[nu][jj];
               BH3BH2bar[j][jj] -= sign*dm0.Trace()*gpropBH3*gpropBH2;
               // BH3 * CD3Bar
               dm0 = dmBH30[mu][j];
               dm0 *= dmCD30Bar[nu][jj];
               dm1 = dmBH31[mu][j];
               dm1 *= dmCD31Bar[nu][jj];
               BH3CD3bar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD3*gpropBH3;
               // BH3 * BH3Bar
               dm0 = dmBH30[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm1 = dmBH31[mu][j];
               dm1 *= dmBH31Bar[nu][jj];
               BH3BH3bar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropBH3*gpropBH3;
            }
         }
         Mfi2[j][jj] = CD2CD2bar[j][jj] + CD2BH2bar[j][jj] + 
                       CD2CD3bar[j][jj] + CD2BH3bar[j][jj] +
                       BH2CD2bar[j][jj] + BH2BH2bar[j][jj] + 
                       BH2CD3bar[j][jj] + BH2BH3bar[j][jj] +
                       CD3CD2bar[j][jj] + CD3BH2bar[j][jj] + 
                       CD3CD3bar[j][jj] + CD3BH3bar[j][jj] +
                       BH3CD2bar[j][jj] + BH3BH2bar[j][jj] + 
                       BH3CD3bar[j][jj] + BH3BH3bar[j][jj];
      }
   }

   // Contract with the incident photon spin density matrix
   Complex_t ampSquared=0;
   for (Int_t j=0; j<2; j++) {
      for (Int_t jj=0; jj<2; jj++) {
         ampSquared += Mfi2[j][jj] * g0->SDM()[j][jj];
      }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: problem with triplets amplitudes:" << std::endl
                << "  These guys should be all real positive:" << std::endl
                << "    ampSquared = " << ampSquared << std::endl
                ;
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

LDouble_t TCrossSection_v1::eeBremsstrahlung(const TLepton &eIn0,
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

   // Obtain the four lepton state matrices
   TDiracMatrix chi0,chi1,chi2,chi3;
   chi0.SetUUbar(e0->Mom(),e0->SDM());
   chi1.SetUUbar(e1->Mom(),e1->SDM());
   chi2.SetUUbar(e2->Mom(),e2->SDM());
   chi3.SetUUbar(e3->Mom(),e3->SDM());

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
   //    {dm}{diag}{leg}{Bar}[mu][j]
   // where
   //    {dm} = dm or some other symbol for Dirac matrix
   //    {diag} = A, B, C, or D to indicate which diagram pair
   //    {leg} = 0 or 1, chain with initial electron 0 or 1
   //    {Bar} = "" or "Bar", representing the matrix or its adjoint partner
   //    [mu] = Lorentz index for photon propagator
   //    [j] = outgoing photon spin index
   // For example, dmB1Bar[3][0] refers to the adjoint variant of the Dirac
   // matrix product (Bar) coming from the leg (leg=1) containing the initial-
   // state electron eIn1, with radiation from the right-hand leg connecting
   // eIn1 to eOut3 (diag=B), Lorentz component mu=3, photon spin component
   // j=0.
   //
   // The plan for computing the sums is as follows:
   //   1. compute all Dirac matrix chains for each leg of each diagram
   //   2. compute the adjoint pair for each of the above (xxxBar matrices)
   //   3. append chi matrices for the u(p) or v(p) spinor factors to each
   //   4. take traces of chains of two xxx[mu][j] and two xxxBar[nu][jj]
   //   5. sum above traces over diag,rad,diagBar,radBar,mu,nu
   //   6. contract sum[j][jj] with the photon spin density matrix
   // The result of the final contraction is |Mfi|^2

   // Pre-compute the electron propagators
   TDiracMatrix dm;
   TDiracMatrix epropA1 = (dm.Slash(e0->Mom() - g0->Mom()) + mLepton) /
                          (-2*g0->Mom().ScalarProd(e0->Mom()));
   TDiracMatrix epropA2 = (dm.Slash(e2->Mom() + g0->Mom()) + mLepton) /
                          (2*g0->Mom().ScalarProd(e2->Mom()));
   TDiracMatrix epropB1 = (dm.Slash(e1->Mom() - g0->Mom()) + mLepton) /
                          (-2*g0->Mom().ScalarProd(e1->Mom()));
   TDiracMatrix epropB2 = (dm.Slash(e3->Mom() + g0->Mom()) + mLepton) /
                          (2*g0->Mom().ScalarProd(e3->Mom()));
   TDiracMatrix epropC1(epropA1);
   TDiracMatrix epropC2(epropB2);
   TDiracMatrix epropD1(epropB1);
   TDiracMatrix epropD2(epropA2);

   // Pre-compute the photon propagators (no 1,2 suffix is needed)
   LDouble_t gpropA = 1/(e1->Mom()-e3->Mom()).InvariantSqr();
   LDouble_t gpropB = 1/(e0->Mom()-e2->Mom()).InvariantSqr();
   LDouble_t gpropC = 1/(e1->Mom()-e2->Mom()).InvariantSqr();
   LDouble_t gpropD = 1/(e0->Mom()-e3->Mom()).InvariantSqr();

   // Evaluate the leading order Feynman amplitude
   const TDiracMatrix gamma0(kDiracGamma0);
   const TDiracMatrix gamma1(kDiracGamma1);
   const TDiracMatrix gamma2(kDiracGamma2);
   const TDiracMatrix gamma3(kDiracGamma3);
   TDiracMatrix gamma[4];
   gamma[0] = gamma0;
   gamma[1] = gamma1;
   gamma[2] = gamma2;
   gamma[3] = gamma3;

   // Compute the product chains of Dirac matrices
   TDiracMatrix dmA0[4][2],dmA0Bar[4][2];
   TDiracMatrix dmA1[4][2],dmA1Bar[4][2];
   TDiracMatrix dmB0[4][2],dmB0Bar[4][2];
   TDiracMatrix dmB1[4][2],dmB1Bar[4][2];
   TDiracMatrix dmC0[4][2],dmC0Bar[4][2];
   TDiracMatrix dmC1[4][2],dmC1Bar[4][2];
   TDiracMatrix dmD0[4][2],dmD0Bar[4][2];
   TDiracMatrix dmD1[4][2],dmD1Bar[4][2];
   for (Int_t j=0; j<2; j++) {
      for (Int_t mu=0; mu<4; mu++) {
         TDiracMatrix current1, current2;

         current1 = gamma[mu];
         current1 *= epropA1;
         current1 *= dm.Slash(g0->EpsStar(j+1));
         current2.Slash(g0->EpsStar(j+1));
         current2 *= epropA2;
         current2 *= gamma[mu];
         dmA0[mu][j] = current1 + current2;
         dmA1[mu][j] = gamma[mu];
 
         current1 = gamma[mu];
         current1 *= epropB1;
         current1 *= dm.Slash(g0->EpsStar(j+1));
         current2.Slash(g0->EpsStar(j+1));
         current2 *= epropB2;
         current2 *= gamma[mu];
         dmB1[mu][j] = current1 + current2;
         dmB0[mu][j] = gamma[mu];

         current1 = gamma[mu];
         current1 *= epropC1;
         current1 *= dm.Slash(g0->EpsStar(j+1));
         current2.Slash(g0->EpsStar(j+1));
         current2 *= epropC2;
         current2 *= gamma[mu];
         dmC0[mu][j] = current1 + current2;
         dmC1[mu][j] = gamma[mu];

         current1 = gamma[mu];
         current1 *= epropD1;
         current1 *= dm.Slash(g0->EpsStar(j+1));
         current2.Slash(g0->EpsStar(j+1));
         current2 *= epropD2;
         current2 *= gamma[mu];
         dmD1[mu][j] = current1 + current2;
         dmD0[mu][j] = gamma[mu];
      }
   }

   // Compute adjoint pairs and append chi matrices
   for (Int_t j=0; j<2; j++) {
      for (Int_t mu=0; mu<4; mu++) {
         dmA0Bar[mu][j] = dmA0[mu][j];
         dmA0Bar[mu][j].Adjoint();
         dmA0Bar[mu][j].UniTransform(gamma0);
         dmA0Bar[mu][j] *= chi2;
         dmA0[mu][j] *= chi0;
         dmA1Bar[mu][j] = dmA1[mu][j];
         // dmA1Bar[mu][j].Adjoint();
         // dmA1Bar[mu][j].UniTransform(gamma0);
         dmA1Bar[mu][j] *= chi3;
         dmA1[mu][j] *= chi1;
 
         dmB0Bar[mu][j] = dmB0[mu][j];
         // dmB0Bar[mu][j].Adjoint();
         // dmB0Bar[mu][j].UniTransform(gamma0);
         dmB0Bar[mu][j] *= chi2;
         dmB0[mu][j] *= chi0;
         dmB1Bar[mu][j] = dmB1[mu][j];
         dmB1Bar[mu][j].Adjoint();
         dmB1Bar[mu][j].UniTransform(gamma0);
         dmB1Bar[mu][j] *= chi3;
         dmB1[mu][j] *= chi1;
 
         dmC0Bar[mu][j] = dmC0[mu][j];
         dmC0Bar[mu][j].Adjoint();
         dmC0Bar[mu][j].UniTransform(gamma0);
         dmC0Bar[mu][j] *= chi3;
         dmC0[mu][j] *= chi0;
         dmC1Bar[mu][j] = dmC1[mu][j];
         // dmC1Bar[mu][j].Adjoint();
         // dmC1Bar[mu][j].UniTransform(gamma0);
         dmC1Bar[mu][j] *= chi2;
         dmC1[mu][j] *= chi1;
 
         dmD0Bar[mu][j] = dmD0[mu][j];
         // dmD0Bar[mu][j].Adjoint();
         // dmD0Bar[mu][j].UniTransform(gamma0);
         dmD0Bar[mu][j] *= chi3;
         dmD0[mu][j] *= chi0;
         dmD1Bar[mu][j] = dmD1[mu][j];
         dmD1Bar[mu][j].Adjoint();
         dmD1Bar[mu][j].UniTransform(gamma0);
         dmD1Bar[mu][j] *= chi2;
         dmD1[mu][j] *= chi1;
      }
   }

   // Finally, the sums over traces
   Complex_t Mfi2[2][2]; 
   Complex_t AAbar[2][2] = {0};
   Complex_t ABbar[2][2] = {0};
   Complex_t ACbar[2][2] = {0};
   Complex_t ADbar[2][2] = {0};
   Complex_t BAbar[2][2] = {0};
   Complex_t BBbar[2][2] = {0};
   Complex_t BCbar[2][2] = {0};
   Complex_t BDbar[2][2] = {0};
   Complex_t CAbar[2][2] = {0};
   Complex_t CBbar[2][2] = {0};
   Complex_t CCbar[2][2] = {0};
   Complex_t CDbar[2][2] = {0};
   Complex_t DAbar[2][2] = {0};
   Complex_t DBbar[2][2] = {0};
   Complex_t DCbar[2][2] = {0};
   Complex_t DDbar[2][2] = {0};
   for (Int_t j=0; j<2; j++) {
      for (Int_t jj=0; jj<2; jj++) {
         for (Int_t mu=0; mu<4; mu++) {
            for (Int_t nu=0; nu<4; nu++) {
               LDouble_t sign;
               sign = (mu == 0)? 1 : -1;
               sign *= (nu == 0)? 1 : -1;
               TDiracMatrix dm0,dm1;
               // A * ABar
               dm0 = dmA0[mu][j];
               dm0 *= dmA0Bar[nu][jj];
               dm1 = dmA1[mu][j];
               dm1 *= dmA1Bar[nu][jj];
               AAbar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropA*gpropA;
               // A * BBar
               dm0 = dmA0[mu][j];
               dm0 *= dmB0Bar[nu][jj];
               dm1 = dmA1[mu][j];
               dm1 *= dmB1Bar[nu][jj];
               ABbar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropA*gpropB;
               // A * CBar
               dm0 = dmA0[mu][j];
               dm0 *= dmC1Bar[nu][jj];
               dm0 *= dmA1[mu][j];
               dm0 *= dmC0Bar[nu][jj];
               ACbar[j][jj] -= sign*dm0.Trace()*gpropA*gpropC;
               // A * DBar
               dm0 = dmA0[mu][j];
               dm0 *= dmD1Bar[nu][jj];
               dm0 *= dmA1[mu][j];
               dm0 *= dmD0Bar[nu][jj];
               ADbar[j][jj] -= sign*dm0.Trace()*gpropA*gpropD;
               // B * ABar
               dm0 = dmB0[mu][j];
               dm0 *= dmA0Bar[nu][jj];
               dm1 = dmB1[mu][j];
               dm1 *= dmA1Bar[nu][jj];
               BAbar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropB*gpropA;
               // B * BBar
               dm0 = dmB0[mu][j];
               dm0 *= dmB0Bar[nu][jj];
               dm1 = dmB1[mu][j];
               dm1 *= dmB1Bar[nu][jj];
               BBbar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropB*gpropB;
               // B * CBar
               dm0 = dmB0[mu][j];
               dm0 *= dmC1Bar[nu][jj];
               dm0 *= dmB1[mu][j];
               dm0 *= dmC0Bar[nu][jj];
               BCbar[j][jj] -= sign*dm0.Trace()*gpropB*gpropC;
               // B * DBar
               dm0 = dmB0[mu][j];
               dm0 *= dmD1Bar[nu][jj];
               dm0 *= dmB1[mu][j];
               dm0 *= dmD0Bar[nu][jj];
               BDbar[j][jj] -= sign*dm0.Trace()*gpropB*gpropD;
               // C * ABar
               dm0 = dmC0[mu][j];
               dm0 *= dmA1Bar[nu][jj];
               dm0 *= dmC1[mu][j];
               dm0 *= dmA0Bar[nu][jj];
               CAbar[j][jj] -= sign*dm0.Trace()*gpropC*gpropA;
               // C * BBar
               dm0 = dmC0[mu][j];
               dm0 *= dmB1Bar[nu][jj];
               dm0 *= dmC1[mu][j];
               dm0 *= dmB0Bar[nu][jj];
               CBbar[j][jj] -= sign*dm0.Trace()*gpropC*gpropB;
               // C * CBar
               dm0 = dmC0[mu][j];
               dm0 *= dmC0Bar[nu][jj];
               dm1 = dmC1[mu][j];
               dm1 *= dmC1Bar[nu][jj];
               CCbar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropC*gpropC;
               // C * DBar
               dm0 = dmC0[mu][j];
               dm0 *= dmD0Bar[nu][jj];
               dm1 = dmC1[mu][j];
               dm1 *= dmD1Bar[nu][jj];
               CDbar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropC*gpropD;
               // D * ABar
               dm0 = dmD0[mu][j];
               dm0 *= dmA1Bar[nu][jj];
               dm0 *= dmD1[mu][j];
               dm0 *= dmA0Bar[nu][jj];
               DAbar[j][jj] -= sign*dm0.Trace()*gpropD*gpropA;
               // D * BBar
               dm0 = dmD0[mu][j];
               dm0 *= dmB1Bar[nu][jj];
               dm0 *= dmD1[mu][j];
               dm0 *= dmB0Bar[nu][jj];
               DBbar[j][jj] -= sign*dm0.Trace()*gpropD*gpropB;
               // D * CBar
               dm0 = dmD0[mu][j];
               dm0 *= dmC0Bar[nu][jj];
               dm1 = dmD1[mu][j];
               dm1 *= dmC1Bar[nu][jj];
               DCbar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropD*gpropC;
               // D * DBar
               dm0 = dmD0[mu][j];
               dm0 *= dmD0Bar[nu][jj];
               dm1 = dmD1[mu][j];
               dm1 *= dmD1Bar[nu][jj];
               DDbar[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropD*gpropD;
            }
         }
         Mfi2[j][jj] = AAbar[j][jj] + ABbar[j][jj] + 
                       ACbar[j][jj] + ADbar[j][jj] +
                       BAbar[j][jj] + BBbar[j][jj] + 
                       BCbar[j][jj] + BDbar[j][jj] +
                       CAbar[j][jj] + CBbar[j][jj] + 
                       CCbar[j][jj] + CDbar[j][jj] +
                       DAbar[j][jj] + DBbar[j][jj] + 
                       DCbar[j][jj] + DDbar[j][jj];
      }
   }

   // Contract with the incident photon spin density matrix
   Complex_t ampSquared=0;
   for (Int_t j=0; j<2; j++) {
      for (Int_t jj=0; jj<2; jj++) {
         ampSquared += Mfi2[j][jj] * g0->SDM()[jj][j];
      }
   }

#if DEBUGGING
   if (real(ampSquared) < 0 || fabs(ampSquared.imag()) > fabs(ampSquared / 1e8L))
   {
      std::cout << "Warning: problem with eeBremsstrahlung amplitudes:" << std::endl
                << "  These guys should be all real positive:" << std::endl
                << "    ampSquared = " << ampSquared << std::endl
                ;
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

void TCrossSection_v1::Streamer(TBuffer &buf)
{
   // All members are static; this function is a noop.
}

void TCrossSection_v1::Print(Option_t *option)
{
   // All members are static; this function is a noop.
}
