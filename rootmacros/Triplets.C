//
// Triplets.C
//
// Calculates the e+e- pair production rate on a free electron target,
// including incident photon polarization effects, for a given set of
// kinematics.  The kinematics are specified by the initial photon
// energy kin, the mass of the e+e- pair M, the recoil momentum vector
// qR (2 transverse components only), the azimuthal angle of the plane
// containing the e+e- pair phi+, and the energy of the pair positron
// E+. The returned value is the differential cross section measured
// microbarns/GeV^4/r per electron, differential in (d^3 qR  dphi+ dE+).
// Another useful expression for the differential measure is
//    (d^3 qR dphi+ dE+) = (M / 2 kin) (dM dqR^2 dphiR dphi+ dE+)
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
#include "sqr.h"

#include <TRandom2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TF1.h>

TRandom2 Triplets_random_gen(0);
TH2D *Triplets_random_bias2D_u0u1 = 0;

//#define H_DIPOLE_FORM_FACTOR 1
LDouble_t FFatomic(LDouble_t qR);

Double_t Triplets(Double_t *var, Double_t *par)
{
   LDouble_t kin=par[0];
   LDouble_t Epos=par[1]=var[0];
   LDouble_t phi12=par[2];
   LDouble_t Mpair=par[3];
   LDouble_t qR2=par[4];
   LDouble_t phiR=par[5];

   // Solve for the 4-vector qR
   if (kin < 0 || Epos < mElectron || Mpair < 2 * mElectron || qR2 < 0) {
      // cout << "no kinematic solution, input parameters out of range"
      //      << endl;
      return 0;
   }
   LDouble_t qR=sqrt(qR2);
   LDouble_t E3=sqrt(qR2+sqr(mElectron));
   LDouble_t costhetaR=(sqr(Mpair)/2 + (kin+mElectron)*(E3-mElectron))/(kin*qR);
   if (fabs(costhetaR) > 1) {
      // cout << "no kinematic solution because |costhetaR| > 1" << endl;
      return 0;
   }
   LDouble_t qRperp=qR*sqrt(1-sqr(costhetaR));
   LDouble_t qRlong=qR*costhetaR;
   TFourVectorReal q3(E3,qRperp*cos(phiR),qRperp*sin(phiR),qRlong);

   // Solve for the c.m. momentum of e+ in the pair 1,2 rest frame
   LDouble_t k12star2=sqr(Mpair/2)-sqr(mElectron);
   if (k12star2 < 0) {
      // cout << "no kinematic solution because k12star2 < 0" << endl;
      return 0;
   }
   LDouble_t k12star=sqrt(k12star2);
   LDouble_t E12=kin+mElectron-E3;
   if (E12 < Mpair) {
      // cout << "no kinematic solution because E12 < Mpair" << endl;
      return 0;
   }
   LDouble_t q12mag=sqrt(sqr(E12)-sqr(Mpair));
   LDouble_t costhetastar=(Epos-E12/2)*Mpair/(k12star*q12mag);
   if (Epos < mElectron) {
      // cout << "no kinematic solution because Epos < mElectron" << endl;
      return 0;
   }
   else if (Epos > E12 - mElectron) {
      // cout << "no kinematic solution because Epos > E12 - mElectron"
      //      << endl;
      return 0;
   }
   else if (fabs(costhetastar) > 1) {
      // cout << "no kinematic solution because |costhetastar| > 1" << endl;
      return 0;
   }
   LDouble_t sinthetastar=sqrt(1-sqr(costhetastar));
   TThreeVectorReal k12(k12star*sinthetastar*cos(phi12),
                        k12star*sinthetastar*sin(phi12),
                        k12star*costhetastar);
   TFourVectorReal q1(Mpair/2,-k12);
   TFourVectorReal q2(Mpair/2,k12);
   TLorentzBoost pairCMtolab(q3[1]/E12,q3[2]/E12,(q3[3]-kin)/E12);
   q1.Boost(pairCMtolab);
   q2.Boost(pairCMtolab);

   // To avoid double-counting, return zero if recoil electron
   // momentum is greater than the momentum of the pair electron.
   if (q2.Length() < qR) {
      // cout << "recoil/pair electrons switched, returning 0" << endl;
      return 0;
   }

   // Define the particle objects
   TPhoton g0;
   TLepton e0(mElectron),e1(mElectron),e2(mElectron),e3(mElectron);
   g0.SetMom(TThreeVectorReal(0,0,kin));
   e0.SetMom(TThreeVectorReal(0,0,0));
   e1.SetMom(q1);
   e2.SetMom(q2);
   e3.SetMom(q3);

   // Set the initial, final polarizations
   g0.SetPol(TThreeVectorReal(1,0,0));
   e0.SetPol(TThreeVectorReal(0,0,0));
   e1.AllPol();
   e2.AllPol();
   e3.AllPol();

   LDouble_t result = TCrossSection::TripletProduction(g0,e0,e1,e2,e3);
   LDouble_t FF = FFatomic(e3.Mom().Length());
   return result * (1 - FF*FF);
}

void set_bias2D_u0u1(TH2D *bias2D)
{
   Triplets_random_bias2D_u0u1 = bias2D;
}

TH2D *get_bias2D_u0u1(TH2D *bias2D)
{
   return Triplets_random_bias2D_u0u1;
}

Int_t demoTriplets(Double_t E0=9.,
                   Double_t Epos=4.5,
                   Double_t phi12=PI_/2,
                   Double_t Mpair=2e-3,
                   Double_t qR2=1e-6,
                   Double_t phiR=0.)
{
   TCanvas *c1 = new TCanvas("c1","Triplet Production Cross Section",200,10,700,500);
   TF1 *f1 = new TF1("f1",Triplets,0,E0,6);
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
      Double_t urand[5];
   } event;
   TString leaflist("E0/D:Epos/D:phi12/D:Mpair/D:qR2/D:phiR/D:thetaR/D:"
                    "diffXS/D:weight/D:weightedXS/D:urand[5]/D");
   event.E0 = kin;

   if (hfile != 0) {
      if (tree == 0) {
         TString title;
         title.Form("e-e+e- triplet production data, Egamma=%f",event.E0);
         tree = new TTree("epairXS",title);
      }
      tree->Branch("event",&event,leaflist,65536);
   }

   LDouble_t sum0=0;
   LDouble_t sum1=0;
   LDouble_t sum2=0;
   for (int n=1; n<=N; n++) { 
      event.weight = 1;
      Triplets_random_gen.RndmArray(5, event.urand);
      if (Triplets_random_bias2D_u0u1) {
         Triplets_random_bias2D_u0u1->GetRandom2(event.urand[0], event.urand[1]);
         static LDouble_t fmean = 0;
         if (fmean == 0) {
            fmean = Triplets_random_bias2D_u0u1->Integral() /
                    Triplets_random_bias2D_u0u1->GetNbinsX() /
                    Triplets_random_bias2D_u0u1->GetNbinsY();
         }
         int i0 = Triplets_random_bias2D_u0u1->GetXaxis()->FindBin(event.urand[0]);
         int i1 = Triplets_random_bias2D_u0u1->GetYaxis()->FindBin(event.urand[1]);
         event.weight = fmean / Triplets_random_bias2D_u0u1->GetBinContent(i0,i1);
      }

      // generate E+ uniform on [0,E0]
      event.Epos = event.urand[2] * event.E0;
      event.weight *= event.E0;
   
      // generate phi12 uniform on [0,2pi]
      event.phi12 = event.urand[3] * 2*PI_;
      event.weight *= 2*PI_;

      // generate phiR uniform on [0,2pi]
      event.phiR = event.urand[4] * 2*PI_;
      event.weight *= 2*PI_;
   
      // generate Mpair with weight (1/M) / (Mcut^2 + M^2)
      LDouble_t Mmin=2*mElectron;
      LDouble_t Mcut=5e-3; // 5 MeV cutoff parameter
      LDouble_t um0 = 1+sqr(Mcut/Mmin);
      LDouble_t um = pow(um0, event.urand[0]);
      event.Mpair = Mcut/sqrt(um-1);
      event.weight *= event.Mpair*(sqr(Mcut)+sqr(event.Mpair))
                      *log(um0)/(2*sqr(Mcut));

      // generate qR^2 with weight (1/qR^2) / sqrt(qRcut^2 + qR^2)
      LDouble_t qRmin = sqr(event.Mpair)/(2*event.E0);
      LDouble_t qRcut = 1e-3; // 1 MeV/c cutoff parameter
      LDouble_t uq0 = qRmin/(qRcut+sqrt(sqr(qRcut)+sqr(qRmin)));
      LDouble_t uq = pow(uq0, event.urand[1]);
      event.qR2 = sqr(2*qRcut*uq/(1-sqr(uq)));
      event.weight *= event.qR2*sqrt(1+event.qR2/sqr(qRcut))
                      *(-2*log(uq0));

      // overall measure Jacobian factor
      event.weight *= event.Mpair/(2*event.E0);

      // compute recoil polar angle thetaR
      LDouble_t E3 = sqrt(event.qR2+sqr(mElectron));
      LDouble_t costhetaR = (sqr(event.Mpair)/2 + (kin+mElectron)*(E3-mElectron)
                           )/(kin*sqrt(event.qR2));
      if (fabs(costhetaR) > 1) {
         // cout << "no kinematic solution because |costhetaR| > 1" << endl;
         continue;
      }
      else {
         event.thetaR = acos(costhetaR);
      }

      Double_t *par=&event.E0;
      Double_t *var=&event.Epos;
      event.diffXS = Triplets(var,par);
      event.weightedXS = event.diffXS*event.weight;
      if (event.weight <= 0) {
         cout << "non-positive event weight " << event.weight << endl;
         continue;
      }
      else if (event.diffXS < 0) {
         cout << "negative differential cross section " << event.diffXS 
              << ", weight = " << event.weight << endl;
         continue;
      }

      if (tree != 0 && event.weightedXS > 0) {
         tree->Fill();
      }

      sum0 += 1;
      sum1 += event.weightedXS;
      sum2 += sqr(event.weightedXS);
      if (n/prescale*prescale == n) {
         cout << "est. total cross section after " << n << " events : "
              << sum1/sum0 << " +/- " << sqrt(sum2-sqr(sum1)/sum0)/sum0 
              << " ub" << endl;
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

LDouble_t FFatomic(LDouble_t qR)
{
   // return the atomic form factor of 4Be normalized to unity
   // at zero momentum transfer qR (GeV/c). Length is in Angstroms.

#if H_DIPOLE_FORM_FACTOR

   LDouble_t a0Bohr = 0.529177 / 1.97327e-6;
   LDouble_t ff = 1 / pow(1 + pow(a0Bohr * qR, 2), 2);

#else

   double Z=4;

   // parameterization given by online database at
   // http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction\
   //      /atomicformfactors/formfactors.php

   LDouble_t acoeff[] = {1.5919, 1.1278, 0.5391, 0.7029};
   LDouble_t bcoeff[] = {43.6427, 1.8623, 103.483, 0.5420};
   LDouble_t ccoeff[] = {0.0385};

   LDouble_t q_invA = qR / 1.97327e-6;
   LDouble_t ff = ccoeff[0];
   for (int i=0; i < 4; ++i) {
      ff += acoeff[i] * exp(-bcoeff[i] * pow(q_invA / (4 * M_PI), 2));
   }
   ff /= Z;

#endif

   return ff;
}
