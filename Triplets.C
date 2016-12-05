//
// Triplets.C
//
// Calculates the e+e- pair production rate on a free electron target,
// including incident photon polarization effects, for a given set of
// kinematics.  The kinematics are specified by the initial photon
// energy kin, the mass of the e+e- pair M, the recoil momentum vector
// qR (2 transverse components only), the azimuthal angle of the plane
// containing the e+e- pair phiR, and the energy of the pair positron
// E+.  The returned value is the differential cross section measured
// microbarns/GeV^4/r, differential in (d^3 qR  dphi+ dE+).  Another 
// useful expression for the differential measure is
//    (d^3 qR dphi- dE-) = (M / 2 kin) (dM dqR^2 dphiR dphi+ dE+)
//
// author: richard.t.jones at uconn.edu
// version: january 1, 2000

#include <iomanip>

#include "Complex.h"
#include "TPhoton.h"
#include "TLepton.h"
#include "TCrossSection.h"
#include "TLorentzBoost.h"
#include "constants.h"

#include <TRandom2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>

inline Double_t sqr(Double_t x) { return x*x; }
inline Complex_t sqr(Complex_t x) { return x*x; }

const TThreeVectorReal zeroVector(0,0,0);
const TThreeVectorReal posXhat(1,0,0);
const TThreeVectorReal negXhat(-1,0,0);
const TThreeVectorReal posYhat(0,1,0);
const TThreeVectorReal negYhat(0,-1,0);
const TThreeVectorReal posZhat(0,0,1);
const TThreeVectorReal negZhat(0,0,-1);

TRandom2 random_gen(0);

Double_t Triplets(Double_t *var, Double_t *par)
{
   Double_t kin=par[0];
   Double_t Epos=par[1];
   Double_t phi12=par[2];
   Double_t Mpair=par[3];
   Double_t qR2=par[4];
   Double_t phiR=par[5]=var[0];

   // Solve for the 4-vector qR
   Double_t qR=sqrt(qR2);
   Double_t E3=sqrt(qR2+sqr(mElectron));
   Double_t costhetaR=(sqr(Mpair)/2 + (kin+mElectron)*(E3-mElectron))/(kin*qR);
   if (fabs(costhetaR) > 1) {
      // cout << "no kinematic solution because |costhetaR| > 1" << endl;
      return 0;
   }
   Double_t qRperp=qR*sqrt(1-sqr(costhetaR));
   Double_t qRlong=qR*costhetaR;
   TFourVectorReal q3(E3,qRperp*cos(phiR),qRperp*sin(phiR),qRlong);

   // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
   Double_t k12star2=sqr(Mpair/2)-sqr(mElectron);
   if (k12star2 < 0) {
      // cout << "no kinematic solution because k12star2 < 0" << endl;
      return 0;
   }
   Double_t k12star=sqrt(k12star2);
   Double_t E12=kin+mElectron-E3;
   Double_t q12mag=sqrt(sqr(E12)-sqr(Mpair));
   Double_t costhetastar=(Epos-E12/2)*Mpair/(k12star*q12mag);
   if (fabs(costhetastar) > 1) {
      // cout << "no kinematic solution because |costhetastar| > 1" << endl;
      return 0;
   }
   Double_t sinthetastar=sqrt(1-sqr(costhetastar));
   TThreeVectorReal k12(k12star*sinthetastar*cos(phi12),
                        k12star*sinthetastar*sin(phi12),
                        k12star*costhetastar);
   TFourVectorReal q1(Mpair/2,-k12);
   TFourVectorReal q2(Mpair/2,k12);
   TLorentzBoost pairCMtolab(q3[1]/E12,q3[2]/E12,(q3[3]-kin)/E12);
   q1.Boost(pairCMtolab);
   q2.Boost(pairCMtolab);

   // Define the particle objects
   TPhoton g0;
   TLepton e0(mElectron),e1(mElectron),e2(mElectron),e3(mElectron);
   TThreeVectorReal p;
   g0.SetMom(p.SetPolar(kin,0,0));
   e0.SetMom(zeroVector);
   e1.SetMom(q1);
   e2.SetMom(q2);
   e3.SetMom(q3);

   // Set the initial, final polarizations
   g0.SetPol(posXhat);
   e0.SetPol(zeroVector);
   e1.AllPol();
   e2.AllPol();
   e3.AllPol();

   Double_t result = TCrossSection::TripletProduction(g0,e0,e1,e2,e3);
   return result;
}

Int_t demoTriplets(Double_t E0=9.,
                   Double_t Epos=4.5,
                   Double_t phi12=PI_/2,
                   Double_t Mpair=2e-3,
                   Double_t qR2=1e-6,
                   Double_t phiR=0.)
{
   TCanvas *c1 = new TCanvas("c1","Triplet Production Cross Section",200,10,700,500);
   TF1 *f1 = new TF1("f1",Triplets,0,2*PI_,6);
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

Int_t genTriplets(Int_t N, Double_t kin=9., TFile *hfile=0, TTree *tree=0, Int_t prescale=1000)
{
   struct event_t {
      Double_t E0;
      Double_t Epos;
      Double_t phi12;
      Double_t Mpair;
      Double_t qR2;
      Double_t phiR;
      Double_t thetaR;
      Double_t diffXS;
      Double_t weight;
      Double_t weightedXS;
   } event;
   TString leaflist("E0/D:Epos/D:phi12/D:Mpair/D:qR2/D:phiR/D:thetaR/D:diffXS/D:weight/D:weightedXS");
   event.E0 = kin;

   if (hfile != 0) {
      if (tree == 0) {
         TString title;
         title.Form("e-e+e- triplet production data, Egamma=%f",event.E0);
         tree = new TTree("epairXS",title);
      }
      tree->Branch("event",&event,leaflist,65536);
   }

   Double_t sum=0;
   Double_t sum2=0;
   for (int n=1; n<=N; n++) { 
      event.weight = 1;

      // generate E+ uniform on [0,E0]
      event.Epos = random_gen.Uniform(event.E0);
      event.weight *= event.E0;
   
      // generate phi12 uniform on [0,2pi]
      event.phi12 = random_gen.Uniform(2*PI_);
      event.weight *= 2*PI_;

      // generate phiR uniform on [0,2pi]
      event.phiR = random_gen.Uniform(2*PI_);
      event.weight *= 2*PI_;
   
      // generate Mpair with weight (1/M) / (Mcut^2 + M^2)
      Double_t Mmin=2*mElectron;
      Double_t Mcut=5e-3; // 5 MeV cutoff parameter
      Double_t um0 = 1+sqr(Mcut/Mmin);
      Double_t um = pow(um0,random_gen.Uniform(1));
      event.Mpair = Mcut/sqrt(um-1);
      event.weight *= event.Mpair*(sqr(Mcut)+sqr(event.Mpair))
                      *log(um0)/(2*sqr(Mcut));

      // generate qR^2 with weight (1/qR^2) / sqrt(qRcut^2 + qR^2)
      Double_t qRmin = sqr(event.Mpair)/(2*event.E0);
      Double_t qRcut = 1e-3; // 1 MeV/c cutoff parameter
      Double_t uq0 = qRmin/(qRcut+sqrt(sqr(qRcut)+sqr(qRmin)));
      Double_t uq = pow(uq0,random_gen.Uniform(1));
      event.qR2 = sqr(2*qRcut*uq/sqr(1-sqr(uq)));
      event.weight *= event.qR2*sqrt(1+event.qR2/sqr(qRcut))
                      *(-2*log(uq0));

      // overall measure Jacobian factor
      event.weight *= event.Mpair/(2*event.E0);

      // compute recoil polar angle thetaR
      Double_t E3 = sqrt(event.qR2+sqr(mElectron));
      Double_t costhetaR = (sqr(event.Mpair)/2 + (kin+mElectron)*(E3-mElectron)
                           )/(kin*sqrt(event.qR2));
      if (fabs(costhetaR) > 1) {
         // cout << "no kinematic solution because |costhetaR| > 1" << endl;
         event.thetaR = 99;
      }
      else {
         event.thetaR = acos(costhetaR);
      }

      Double_t *par=&event.E0;
      Double_t *var=&event.phiR;
      event.diffXS = Triplets(var,par);
      event.weightedXS = event.diffXS*event.weight;

      if (tree != 0) {
         tree->Fill();
      }

      sum += event.weightedXS;
      sum2 += sqr(event.weightedXS);
      if (n/prescale*prescale == n) {
         cout << "est. total cross section after " << n << " events : "
              << sum/n << " +/- " << sqrt(sum2-sqr(sum)/n)/n << " ub" << endl;
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
