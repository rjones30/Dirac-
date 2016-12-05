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
//    			                       2
//    				        | M   |
//    				    4      fi
//    		d[sigma]  =  ( 2pi )   --------- d[rho(final)]
//    		                         F(in)
//
// 3. F(in) is equal to the product [ 4 E(beam) E(target) ] times the
//    relative velocity between beam and target in the initial state.
//
// 4. rho(final) is calculated in whatever frame the user has chosen
//    to specify the kinematics.  It consists of a term of the form
//		   -1     -3   3
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
//	final-state fermion, or Vbar(p,s) if it ends with an initial-
//	state antifermion, using the appropriate p,s for this state;
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

#include <iostream>
#include "TCrossSection.h"

ClassImp(TCrossSection)

#include "TPhoton.h"
#include "TLepton.h"
#include "TLorentzBoost.h"

const Double_t PI_=2*atan2(1.,0.);

inline Double_t sqr(Double_t x) { return x*x; }
inline Complex_t sqr(Complex_t x) { return x*x; }

Double_t TCrossSection::Compton(const TPhoton &gIn, const TLepton &eIn,
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
   const Double_t mLepton=eI->Mass();

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

   const Double_t fluxIn = 4*gI->Mom()[0]*(eI->Mom().Length()+eI->Mom()[0]);
   const Double_t rhoFin = sqr(gF->Mom()[0])/eF->Mom().ScalarProd(gF->Mom())/4;
   const Double_t kinFactor = 4*rhoFin/fluxIn;

   Double_t diffXsect = hbarcSqr*sqr(alphaQED)*real(ampSquared)*kinFactor;

   return diffXsect;

   // The unpolarized Klein Nishina formula is here for comparison
   Double_t sinSqrTheta = 1 - sqr(gF->Mom()[3]/gF->Mom()[0]);
   Double_t KleinNishinaResult = sqr(alphaQED/mLepton)/2;
   KleinNishinaResult *= sqr(gF->Mom()[0]/gI->Mom()[0]);
   KleinNishinaResult *= (gF->Mom()[0]/gI->Mom()[0]) +
                         (gI->Mom()[0]/gF->Mom()[0]) - sinSqrTheta;
   KleinNishinaResult *= hbarcSqr;
   return KleinNishinaResult;
}

Double_t TCrossSection::Bremsstrahlung(
                        const TLepton &eIn, const TLepton &eOut,
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
   const Double_t mLepton=eI->Mass();
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

   Double_t kinFactor = 1/sqr(2*PI_*eI->Mom()[0]); // |qRecoil| << mElectron
   Double_t diffXsect = hbarcSqr*pow(alphaQED,3)*real(ampSquared)
                       *kinFactor/sqr(qRecoil.InvariantSqr());
   return diffXsect;
}

Double_t TCrossSection::PairProduction(
                        const TPhoton &gIn,
                        const TLepton &eOut, const TLepton &pOut)
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
   const Double_t mLepton=eF->Mass();
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
         ampSquared += (invAmp[j]*invAmpBar[jj]).Trace() * gI->SDM()[jj][j];
      }
   }

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

   Double_t kinFactor = 1/sqr(2*PI_*gI->Mom()[0]);
   Double_t diffXsect = hbarcSqr*pow(alphaQED,3)*real(ampSquared)
                       *kinFactor/sqr(qRecoil.InvariantSqr());
   return diffXsect;
}

Double_t TCrossSection::TripletProduction(
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
   // Units are microbarns/GeV^4/r.

   TPhoton gIncoming(gIn), *g0=&gIncoming;
   TLepton eIncoming(eIn), *e0=&eIncoming;
   TLepton pOutgoing(pOut), *e1=&pOutgoing;
   TLepton eOutgoing2(eOut2), *e2=&eOutgoing2;
   TLepton eOutgoing3(eOut3), *e3=&eOutgoing3;

   const Double_t mLepton=e0->Mass();

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
   //	1. compute all Dirac matrix chains for each leg of each diagram
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
   TDiracMatrix epropCD2b = (dm.Slash(e3->Mom() - g0->Mom()) + mLepton) /
                            (-2*g0->Mom().ScalarProd(e3->Mom()));
   TDiracMatrix epropBH2a = (dm.Slash(g0->Mom() - e1->Mom()) + mLepton) /
                            (-2*g0->Mom().ScalarProd(e1->Mom()));
   TDiracMatrix epropBH2b = (dm.Slash(e2->Mom() - g0->Mom()) + mLepton) /
                            (-2*g0->Mom().ScalarProd(e2->Mom()));
   TDiracMatrix epropCD3a(epropCD2a);
   TDiracMatrix epropCD3b(epropBH2b);
   TDiracMatrix epropBH3a(epropBH2a);
   TDiracMatrix epropBH3b(epropCD2b);

   // Pre-compute the photon propagators (no a,b suffix needed)
   Double_t gpropCD2 = 1/(e1->Mom()+e2->Mom()).InvariantSqr();
   Double_t gpropBH2 = 1/(e0->Mom()-e3->Mom()).InvariantSqr();
   Double_t gpropCD3 = 1/(e1->Mom()+e3->Mom()).InvariantSqr();
   Double_t gpropBH3 = 1/(e0->Mom()-e2->Mom()).InvariantSqr();

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
         dmCD20Bar[mu][j] *= chi3;
         dmCD20[mu][j] *= chi0;
         dmCD21Bar[mu][j] = dmCD21[mu][j];
         // dmCD21Bar[mu][j].Adjoint();
         // dmCD21Bar[mu][j].UniTransform(gamma0);
         dmCD21Bar[mu][j] *= chi2;
         dmCD21[mu][j] *= chi1;
 
         dmBH20Bar[mu][j] = dmBH20[mu][j];
         // dmBH20Bar[mu][j].Adjoint();
         // dmBH20Bar[mu][j].UniTransform(gamma0);
         dmBH20Bar[mu][j] *= chi3;
         dmBH20[mu][j] *= chi0;
         dmBH21Bar[mu][j] = dmBH21[mu][j];
         dmBH21Bar[mu][j].Adjoint();
         dmBH21Bar[mu][j].UniTransform(gamma0);
         dmBH21Bar[mu][j] *= chi2;
         dmBH21[mu][j] *= chi1;
 
         dmCD30Bar[mu][j] = dmCD30[mu][j];
         dmCD30Bar[mu][j].Adjoint();
         dmCD30Bar[mu][j].UniTransform(gamma0);
         dmCD30Bar[mu][j] *= chi2;
         dmCD30[mu][j] *= chi0;
         dmCD31Bar[mu][j] = dmCD31[mu][j];
         // dmCD31Bar[mu][j].Adjoint();
         // dmCD31Bar[mu][j].UniTransform(gamma0);
         dmCD31Bar[mu][j] *= chi3;
         dmCD31[mu][j] *= chi1;
 
         dmBH30Bar[mu][j] = dmBH30[mu][j];
         // dmBH30Bar[mu][j].Adjoint();
         // dmBH30Bar[mu][j].UniTransform(gamma0);
         dmBH30Bar[mu][j] *= chi2;
         dmBH30[mu][j] *= chi0;
         dmBH31Bar[mu][j] = dmBH31[mu][j];
         dmBH31Bar[mu][j].Adjoint();
         dmBH31Bar[mu][j].UniTransform(gamma0);
         dmBH31Bar[mu][j] *= chi3;
         dmBH31[mu][j] *= chi1;
      }
   }

   // Finally, the sums over traces
   Complex_t Mfi2[2][2]; 
   Mfi2[0][0] = Mfi2[0][1] = Mfi2[1][0] = Mfi2[1][1] = Complex_t(0,0);
   for (Int_t j=0; j<2; j++) {
      for (Int_t jj=0; jj<2; jj++) {
         for (Int_t mu=0; mu<4; mu++) {
            for (Int_t nu=0; nu<4; nu++) {
               Double_t sign;
               sign = (mu == 0)? 1 : -1;
               sign *= (nu == 0)? 1 : -1;
               TDiracMatrix dm0,dm1;
               // CD2 * CD2Bar
               dm0 = dmCD20[mu][j];
               dm0 *= dmCD20Bar[nu][jj];
               dm1 = dmCD21[mu][j];
               dm1 *= dmCD21Bar[nu][jj];
               Mfi2[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD2*gpropCD2;
               // CD2 * BH2Bar
               dm0 = dmCD20[mu][j];
               dm0 *= dmBH20Bar[nu][jj];
               dm1 = dmCD21[mu][j];
               dm1 *= dmBH21Bar[nu][jj];
               Mfi2[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD2*gpropBH2;
               // CD2 * CD3Bar
               dm0 = dmCD20[mu][j];
               dm0 *= dmCD30Bar[nu][jj];
               dm0 *= dmCD21[mu][j];
               dm0 *= dmCD31Bar[nu][jj];
               Mfi2[j][jj] -= sign*dm0.Trace()*gpropCD2*gpropCD3;
               // CD2 * BH3Bar
               dm0 = dmCD20[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm0 *= dmCD21[mu][j];
               dm0 *= dmBH31Bar[nu][jj];
               Mfi2[j][jj] -= sign*dm0.Trace()*gpropCD2*gpropBH3;
               // BH2 * CD2Bar
               dm0 = dmBH20[mu][j];
               dm0 *= dmCD20Bar[nu][jj];
               dm1 = dmBH21[mu][j];
               dm1 *= dmCD21Bar[nu][jj];
               Mfi2[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD2*gpropBH2;
               // BH2 * BH2Bar
               dm0 = dmBH20[mu][j];
               dm0 *= dmBH20Bar[nu][jj];
               dm1 = dmBH21[mu][j];
               dm1 *= dmBH21Bar[nu][jj];
               Mfi2[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropBH2*gpropBH2;
               // BH2 * CD3Bar
               dm0 = dmBH20[mu][j];
               dm0 *= dmCD30Bar[nu][jj];
               dm0 *= dmBH21[mu][j];
               dm0 *= dmCD31Bar[nu][jj];
               Mfi2[j][jj] -= sign*dm0.Trace()*gpropBH2*gpropCD3;
               // BH2 * BH3Bar
               dm0 = dmBH20[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm0 *= dmBH21[mu][j];
               dm0 *= dmBH31Bar[nu][jj];
               Mfi2[j][jj] -= sign*dm0.Trace()*gpropBH2*gpropBH3;
               // CD3 * CD2Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmCD20Bar[nu][jj];
               dm0 *= dmCD31[mu][j];
               dm0 *= dmCD21Bar[nu][jj];
               Mfi2[j][jj] -= sign*dm0.Trace()*gpropCD2*gpropCD3;
               // CD3 * BH2Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmBH20Bar[nu][jj];
               dm0 *= dmCD31[mu][j];
               dm0 *= dmBH21Bar[nu][jj];
               Mfi2[j][jj] -= sign*dm0.Trace()*gpropBH2*gpropCD3;
               // CD3 * CD3Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmCD30Bar[nu][jj];
               dm1 = dmCD31[mu][j];
               dm1 *= dmCD31Bar[nu][jj];
               Mfi2[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD3*gpropCD3;
               // CD3 * BH3Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm1 = dmCD31[mu][j];
               dm1 *= dmBH31Bar[nu][jj];
               Mfi2[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD3*gpropBH3;
               // BH3 * CD2Bar
               dm0 = dmBH30[mu][j];
               dm0 *= dmCD20Bar[nu][jj];
               dm0 *= dmBH31[mu][j];
               dm0 *= dmCD21Bar[nu][jj];
               Mfi2[j][jj] -= sign*dm0.Trace()*gpropBH3*gpropCD2;
               // BH3 * BH2Bar
               dm0 = dmBH30[mu][j];
               dm0 *= dmBH20Bar[nu][jj];
               dm0 *= dmBH31[mu][j];
               dm0 *= dmBH21Bar[nu][jj];
               Mfi2[j][jj] -= sign*dm0.Trace()*gpropBH3*gpropBH2;
               // BH3 * CD3Bar
               dm0 = dmCD30[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm1 = dmCD31[mu][j];
               dm1 *= dmBH31Bar[nu][jj];
               Mfi2[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropCD3*gpropBH3;
               // BH3 * BH3Bar
               dm0 = dmBH30[mu][j];
               dm0 *= dmBH30Bar[nu][jj];
               dm1 = dmBH31[mu][j];
               dm1 *= dmBH31Bar[nu][jj];
               Mfi2[j][jj] += sign*dm0.Trace()*dm1.Trace()*gpropBH3*gpropBH3;
            }
         }
      }
   }

   // Contract with the incident photon spin density matrix
   Complex_t ampSquared=0;
   for (Int_t j=0; j<2; j++) {
      for (Int_t jj=0; jj<2; jj++) {
         ampSquared += Mfi2[j][jj] * g0->SDM()[jj][j];
      }
   }

   // Obtain the kinematical factors:
   //    (1) 1/flux from initial state 1/(4 kin [p0 + E0])
   //    (2) rho from density of final states factor
   // where the general relativistic expression for rho is
   //  rho = pow(2*PI_,4-3*N) delta4(Pin-Pout) [d4 P1] [d4 P2] ... [d4 PN]
   // using differential forms [d4 P] = d4P delta(P.P - m*m) where P.P is
   // the invariant norm of four-vector P, m is the known mass of the
   // corresponding particle, and N is the number of final state particles.
   //    (3) absorb three powers of 4*PI_ into pow(alphaQED,3)

   Double_t fluxFactor = 4*g0->Mom()[0]*(e0->Mom().Length()+e0->Mom()[0]);
   Double_t rhoFactor = 1/(8*e3->Mom()[0]*(e1->Mom()+e2->Mom()).Length());
   Double_t piFactor = pow(2*PI_,4-9)*pow(4*PI_,3);
   Double_t diffXsect = hbarcSqr * pow(alphaQED,3) * real(ampSquared)
                        / fluxFactor * rhoFactor * piFactor;
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
