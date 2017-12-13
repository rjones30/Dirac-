//
// Pairs.C
//
// Calculates the e+e- pair production rate including polarization
// effects for a given kinematics on an atomic target, where the atomic
// structure is represented by a scalar form factor.  For the case of
// fermionic targets, see Triplets.C for scattering from free electrons
// or BetheHeitler.C for scattering from a free nucleon   The kinematics
// are set by the initial photon energy kin, the mass of the e+e- pair M,
// the recoil momentum squared qR^2, the azimuthal angle of the recoil
// momentum phiR, the azimuthal angle of the outgoing positron about the
// pair momentum axis phi+, and energy of the outgoing positron E+.  The
// returned value is the differential cross section measured in
// microbarns/GeV^4/r, differential in (d^3 qR  dphi- dE-). Another
// useful expression for the differential measure is
//    (d^3 qR dphi- dE-) = (M / 2 kin) (dM dqR^2 dphiR dphi- dE-)
// The cross section contains the form factor squared that is supposed 
// to represent the carbon atom.
//
// author: richard.t.jones at uconn.edu
// version: january 1, 2000

#include <iomanip>

#include "Complex.h"
#include "TPhoton.h"
#include "TLepton.h"
#include "TLorentzBoost.h"
#include "TCrossSection.h"
#include "constants.h"

#include <TRandom2.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

#ifndef DEFINE_SQR_ON_STANDARD_TYPES
#define DEFINE_SQR_ON_STANDARD_TYPES
inline unsigned int sqr(unsigned int x) { return x*x; }
inline Int_t sqr(Int_t x) { return x*x; }
inline Float_t sqr(Float_t x) { return x*x; }
inline Double_t sqr(Double_t x) { return x*x; }
inline LDouble_t sqr(LDouble_t x) { return x*x; }
inline Complex_t sqr(Complex_t x) { return x*x; }
#endif

#ifndef STANDARD_VECTOR_CONSTANTS
#define STANDARD_VECTOR_CONSTANTS
const TThreeVectorReal zeroVector(0,0,0);
const TThreeVectorReal posXhat(1,0,0);
const TThreeVectorReal negXhat(-1,0,0);
const TThreeVectorReal posYhat(0,1,0);
const TThreeVectorReal negYhat(0,-1,0);
const TThreeVectorReal posZhat(0,0,1);
const TThreeVectorReal negZhat(0,0,-1);
#endif

#ifndef DEFINE_GLOBAL_RANDOM_GEN
#define DEFINE_GLOBAL_RANDOM_GEN
TRandom2 random_gen(0);
#endif

Double_t Pairs(Double_t *var, Double_t *par)
{
   LDouble_t kin=par[0];
   LDouble_t Epos=par[1]=var[0];
   LDouble_t phi12=par[2];
   LDouble_t Mpair=par[3];
   LDouble_t qR2=par[4];
   LDouble_t phiR=par[5];

   TPhoton gIn;
   TLepton eOut(mElectron), pOut(mElectron);

   // Solve for the rest of the kinematics, limit of high mass target
   LDouble_t qR=sqrt(qR2);
   LDouble_t Eele=kin-Epos;
   LDouble_t qmin=kin-sqrt(sqr(kin)-sqr(Mpair));
   LDouble_t costhetaR=(qR2+sqr(Mpair))/(2*kin*qR);
   if (costhetaR > 1) {
      // std::cout << "no kinematic solution because costhetaR > 1" << std::endl;
      return 0;
   }
   LDouble_t sinthetaR=sqrt(1-sqr(costhetaR));
   TThreeVectorReal qRecoil(qR*sinthetaR*cos(phiR),
                            qR*sinthetaR*sin(phiR),
                            qR*costhetaR);

   if (qRecoil[3] < qmin) {
      // std::cout << "no kinematic solution because qRecoil[3] < qmin" << std::endl;
      return 0;
   }

   TThreeVectorReal p;
   gIn.SetMom(p.SetPolar(kin,0,0));
   LDouble_t pStar2=sqr(Mpair/2)-sqr(mElectron);
   if (pStar2 < 0) {
      // std::cout << "no kinematic solution because pStar2 < 0" << std::endl;
      return 0;
   }
   LDouble_t pStar=sqrt(pStar2);
   LDouble_t p12mag=sqrt(sqr(kin)-sqr(Mpair));
   LDouble_t costhetastar=(Epos-kin/2)*Mpair/(pStar*p12mag);
   if (fabs(costhetastar) > 1) {
      // std::cout << "no kinematic solution because costhetastar < 1" << std::endl;
      return 0;
   }
   LDouble_t sinthetastar=sqrt(1-sqr(costhetastar));
   TThreeVectorReal k12(pStar*sinthetastar*cos(phi12),
                        pStar*sinthetastar*sin(phi12),
                        pStar*costhetastar);
   TFourVectorReal p1(Mpair/2,k12);
   TLorentzBoost toLab(qRecoil[1]/kin,qRecoil[2]/kin,(qRecoil[3]-kin)/kin);
   p1.Boost(toLab);
   pOut.SetMom(p1);
   TThreeVectorReal p2(gIn.Mom()-qRecoil-p1);
   eOut.SetMom(p2);

   // Set the initial,final polarizations
   gIn.SetPol(posXhat);
   eOut.AllPol();
   pOut.AllPol();

   // Multiply the basic cross section by the carbon atomic form factor
   // Here I chose a simple parameterization for the form factor
   const LDouble_t Z=6;
   const LDouble_t beta2=111*pow(Z,-1/3.)/mElectron;    // ff cutoff in /GeV
   const LDouble_t Fff=1/(1+sqr(beta2)*qRecoil.LengthSqr());
   LDouble_t result = TCrossSection::PairProduction(gIn,eOut,pOut);
   result *= sqr(Z*(1-Fff));
   return result;

   // The unpolarized Bethe-Heitler cross section is given here for comparison
   LDouble_t delta=136*mElectron/pow(Z,0.33333) * kin/(Eele*Epos);
   LDouble_t aCoul=sqr(alphaQED*Z);
   LDouble_t fCoul=aCoul*(1/(1+aCoul)+0.20206-0.0369*aCoul
                       +0.0083*pow(aCoul,2)-0.002*pow(aCoul,3));
   LDouble_t xsi=log(1440/pow(Z,0.66667))/(log(183/pow(Z,0.33333)-fCoul));
   LDouble_t FofZ=(8./3.)*log(Z) + ((kin < 0.05)? 0 : 8*fCoul);
   LDouble_t Phi1=20.867-3.242*delta+0.625*sqr(delta);
   LDouble_t Phi2=20.209-1.930*delta-0.086*sqr(delta);
   LDouble_t Phi0=21.12-4.184*log(delta+0.952);
   if (delta > 1) {
      Phi1=Phi2=Phi0;
   }
   result = hbarcSqr/sqr(mElectron)*pow(alphaQED,3)/kin
            * Z*(Z+xsi)
            * (
                 (sqr(Eele)+sqr(Epos))/sqr(kin)*(Phi1-FofZ/2)
                 + (2./3.)*(Eele*Epos)/sqr(kin)*(Phi2-FofZ/2)
              );
}

Int_t demoPairs(Double_t E0=9.,
                Double_t Epos=4.5,
                Double_t phi12=0,
                Double_t Mpair=2e-3,
                Double_t qR2=1e-6,
                Double_t phiR=0.)
{
   TCanvas *c1 = new TCanvas("c1","Pair Production Rate",200,10,700,500);
   TF1 *f1 = new TF1("f1",Pairs,0,Epos,6);
   Double_t params[6];
   params[0] = E0;
   params[1] = Epos;
   params[2] = phi12;
   params[3] = Mpair;
   params[4] = qR2;
   params[5] = phiR;
   f1->SetParameter(0,params[0]);
   f1->SetParameter(1,params[1]);
   f1->SetParameter(2,params[2]);
   f1->SetParameter(3,params[3]);
   f1->SetParameter(4,params[4]);
   f1->SetParameter(5,params[5]);
   //f1->DrawCopy("same");
   f1->Draw();
   c1->Update();
   return 0;
}

Int_t genPairs(Int_t N, Double_t kin=9., TFile *hfile=0, TTree *tree=0, Int_t prescale=1000)
{
   struct event_t {
      Double_t E0;
      Double_t Epos;
      Double_t phi12;
      Double_t Mpair;
      Double_t qR2;
      Double_t phiR;
      Double_t diffXS;
      Double_t weight;
      Double_t weightedXS;
   } event;
   TString leaflist("E0/D:Epos/D:phi12/D:Mpair/D:qR2/D:phiR/D:diffXS/D:weight/D:weightedXS");
   event.E0 = kin;

   if (hfile != 0) {
      if (tree == 0) {
         TString title;
         title.Form("e+e- pair production data, Egamma=%f",event.E0);
         tree = new TTree("epairXS",title);
      }
      tree->Branch("event",&event,leaflist,65536);
   }

   LDouble_t sum=0;
   LDouble_t sum2=0;
   for (int n=1; n<=N; n++) { 
      event.weight = 1;

      // generate Epos uniform on [0,E0]
      event.Epos = random_gen.Uniform(event.E0);
      event.weight *= event.E0;
   
      // generate phi12 uniform on [0,2pi]
      event.phi12 = random_gen.Uniform(2*PI_);
      event.weight *= 2*PI_;

      // generate phiR uniform on [0,2pi]
      event.phiR = random_gen.Uniform(2*PI_);
      event.weight *= 2*PI_;
   
#ifdef OLD_WEIGHTING

      // generate Mpair with weight (M0/M)^3
      LDouble_t M0=2*mElectron;
      event.Mpair = M0/sqrt(random_gen.Uniform(1.));
      event.weight *= pow(event.Mpair,3)/(2*M0*M0);

      // generate qR2 with weight 1/(q02 + qR2)^2
      LDouble_t q02=sqr(5e-5); // 50 keV/c cutoff parameter
      LDouble_t u=random_gen.Uniform(1);
      event.qR2 = q02*(1-u)/(u+1e-50);
      event.weight *= sqr(event.qR2+q02)/q02;

#else

      // generate Mpair with weight (1/M) / (Mcut^2 + M^2)
      LDouble_t Mmin=2*mElectron;
      LDouble_t Mcut=5e-3; // 5 MeV cutoff parameter
      LDouble_t um0 = 1+sqr(Mcut/Mmin);
      LDouble_t um = pow(um0,random_gen.Uniform(1));
      event.Mpair = Mcut/sqrt(um-1);
      event.weight *= event.Mpair*(sqr(Mcut)+sqr(event.Mpair))
                      *log(um0)/(2*sqr(Mcut));

      // generate qR^2 with weight (1/qR^2) / sqrt(qRcut^2 + qR^2)
      LDouble_t qRmin = sqr(event.Mpair)/(2*event.E0);
      LDouble_t qRcut = 1e-3; // 1 MeV/c cutoff parameter
      LDouble_t uq0 = qRmin/(qRcut+sqrt(sqr(qRcut)+sqr(qRmin)));
      LDouble_t uq = pow(uq0,random_gen.Uniform(1));
      event.qR2 = sqr(2*qRcut*uq/(1-sqr(uq)));
      event.weight *= event.qR2*sqrt(1+event.qR2/sqr(qRcut))
                      *(-2*log(uq0));

#endif

      // overall measure Jacobian factor
      event.weight *= event.Mpair/(2*event.E0);

      Double_t *par=&event.E0;
      Double_t *var=&event.Epos;
      event.diffXS = Pairs(var,par);
      event.weightedXS = event.diffXS*event.weight;

      if (tree != 0 && event.weightedXS > 0) {
         tree->Fill();
      }

      sum += event.weightedXS;
      sum2 += sqr(event.weightedXS);
      if (n/prescale*prescale == n) {
         std::cout << "est. total cross section after " << n << " events : "
              << sum/n << " +/- " << sqrt(sum2-sqr(sum)/n)/n << " ub" << std::endl;
      }
   }

   if (tree != 0) {
      tree->FlushBaskets();
   }
   if (hfile != 0) {
      hfile->Write();
   }
   return 0;
}
