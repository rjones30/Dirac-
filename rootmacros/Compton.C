//
// Compton.C
//
// Calculates the Compton differential cross section in the rest frame
// of the initial electron, assumed to be free.  Units are microbarns
// per unit solid angle of the scattered photon.
//
// author: richard.t.jones at uconn.edu
// version: january 1, 2000

#include "Complex.h"
#include "TPhoton.h"
#include "TLepton.h"
#include "TCrossSection.h"
#include "constants.h"
#include "sqr.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TAxis.h>

Double_t Compton(Double_t *var, Double_t *par)
{
   LDouble_t theta=var[0];
   LDouble_t kin=par[0];
   LDouble_t phi=par[1];
   Int_t siggIn=par[3];
   Int_t sigeIn=par[4];

   TLepton eIn(mElectron), eOut(mElectron);
   TPhoton gIn, gOut;

   // Solve for the rest of the kinematics
   LDouble_t kout = kin/(1+(kin/mElectron)*(1-cos(theta)));
   TThreeVectorReal p;
   gIn.SetMom(TThreeVectorReal(0,0,kin));
   eIn.SetMom(TThreeVectorReal(0,0,0));
   gOut.SetMom(p.SetPolar(kout,theta,phi));
   eOut.SetMom(gIn.Mom()+eIn.Mom()-gOut.Mom());

   // Set the initial,final polarizations
   switch (siggIn) {
    case -1:
      gIn.SetPol(TThreeVectorReal(0,0,-1));
      break;
    case  0:
      gIn.SetPol(TThreeVectorReal(0,0,0));
      break;
    case +1:
      gIn.SetPol(TThreeVectorReal(0,0,1));
      break;
    case +2:
      gIn.SetPol(TThreeVectorReal(0,1,0));
      break;
    default:
      std::cout << "Compton.C :"
           << "bad photon helicity specified, default to zero!"
           << std::endl;
      gIn.SetPol(TThreeVectorReal(0,0,0));
   }
   switch (sigeIn) {
    case -1:
      eIn.SetPol(TThreeVectorReal(0,0,-1));
      break;
    case  0:
      eIn.SetPol(TThreeVectorReal(0,0,0));
      break;
    case +1:
      eIn.SetPol(TThreeVectorReal(0,0,1));
      break;
    default:
      std::cout << "Compton.C :"
           << "bad electron helicity specified, default to zero!"
           << std::endl;
      eIn.SetPol(TThreeVectorReal(0,0,0));
   }
   gOut.AllPol();
   eOut.AllPol();

   LDouble_t result=TCrossSection::Compton(gIn,eIn,gOut,eOut);
   return result;
}

Double_t ComptonAsym(Double_t *var, Double_t *par)
{
   LDouble_t sigpp,sigpm,sigmp,sigmm;
   LDouble_t sigp0,sigm0,sig0p,sig0m;
   par[3] = +1;
   par[4] = +1;
   sigpp = Compton(var,par);
   par[3] = +1;
   par[4] = -1;
   sigpm = Compton(var,par);
   par[3] = -1;
   par[4] = +1;
   sigmp = Compton(var,par);
   par[3] = -1;
   par[4] = -1;
   sigmm = Compton(var,par);
   return (sigpm - sigmm) / (sigpm + sigmm);
}

Int_t demoCompton(Double_t Ephot, Double_t phi)
{
   TCanvas *c1 = new TCanvas("c1","Compton Cross Section",200,10,700,500);
   TF1 *comp = new TF1("comp",Compton,0,3.1416,5);
   comp->SetParameter(0,Ephot);
   comp->SetParameter(1,phi);
   comp->SetParameter(3,+2);
   comp->SetParameter(4,0);
   comp->GetXaxis()->SetTitle("#theta (radians)");
   comp->GetYaxis()->SetTitle("d#sigma/d#Omega (#mub)");
   comp->GetYaxis()->SetTitleOffset(1.5);
   comp->Draw();
   c1->Update();
   return 0;
}

Int_t demoComptonAsym()
{
   TCanvas *c1 = new TCanvas("c1","Compton Cross Section",200,10,700,500);
   TF1 *asym = new TF1("asym",ComptonAsym,0,3.1416,5);
   asym->SetParameter(0,2e-5);
   asym->SetParameter(1,0);
   asym->SetParameter(2,1.2);
   asym->Draw();
   c1->Update();
   return 0;
}

Int_t demoComptonAsym(Double_t Ephot)
{
   TF1 *asym2 = new TF1("asym2",ComptonAsym,0,3.1416,5);
   asym2->SetParameter(0,Ephot);
   asym2->SetParameter(1,0);
   asym2->SetParameter(2,1.2);
   asym2->DrawCopy("SAME");
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   c1->Update();
   return 0;
}

Double_t ComptonBackScatter(Double_t *var, Double_t *par)
{
   LDouble_t kf=var[0];
   LDouble_t ki=par[0];
   LDouble_t phi=par[1];
   LDouble_t E0=par[2];

   TLepton eIn(mElectron), eOut(mElectron);
   TPhoton gIn, gOut;

   // Solve for the rest of the kinematics - in rest frame of reaction
   LDouble_t P0 = sqrt(E0*E0-mElectron*mElectron);
   LDouble_t S = mElectron*mElectron+2*ki*(E0+P0);
   LDouble_t rootS = sqrt(S);
   LDouble_t kstar = (S-mElectron*mElectron)/(2*rootS);
   TThreeVectorReal k(0,0,-kstar);
   gIn.SetMom(k);
   TThreeVectorReal p(0,0,kstar);
   eIn.SetMom(p);
   LDouble_t gamma = (E0+ki)/rootS;
   LDouble_t eta = (P0-ki)/rootS;
   LDouble_t costheta = (1-kf/(gamma*kstar))*gamma/eta;
   if (fabs(costheta) > 1.0) { return 0; }
   LDouble_t theta = acos(costheta);
   gOut.SetMom(k.SetPolar(kstar,theta+PI_,phi));
   eOut.SetMom(p.SetPolar(kstar,theta,phi));

   // Set the initial,final polarizations
   gIn.SetPol(TThreeVectorReal(0,0,0));
   eIn.SetPol(TThreeVectorReal(0,0,0));
   gOut.AllPol();
   eOut.AllPol();

   LDouble_t result=TCrossSection::Compton(gIn,eIn,gOut,eOut);

   result*=(2*PI_)/(eta*kstar);    // convert cross section from d(sigma)/d(Omega*)
                                   // to d(sigma)/d(kf)

   LDouble_t P=600;      // laser power in Watts (peak times duty factor)
   LDouble_t G=250;      // laser cavity gain factor
   LDouble_t tau=2e-12;  // laser pulse length (s) times crossing factor
   LDouble_t rC=10e-6;   // neck radius of cavity beam (m)
   LDouble_t rB=10e-6;   // electron beam radius (m)
   LDouble_t I=1.0e-6;   // electron beam current (A)
   LDouble_t L=0.0450;   // effective length of cavity (m)
   Int_t    N=2;         // number of passes through beam
   Int_t    pulsed=1;    // indicate whether laser is pulsed

   result*=1e-34;        // convert from microbarns to m^2
   result*=N*P*G/(ki*1.6e-10);    // ki is converted from GeV to Joules
   result/=2*PI_*(rB*rB+rC*rC);   // area of overlap region
   result*=I/1.6e-19;    // rate of electrons in beam
   if (pulsed)
     result*=tau;        // duration of pulse crossing (lab frame)
   else
     result*=L/3e8;      // time spent in crossing region (lab frame)

   return result;
}

Int_t demoComptonBackScatter()
{
   TCanvas *c2 = new TCanvas("c2","Compton Backscatter Cross Section",200,10,700,500);
   TF1 *f2 = new TF1("f2",ComptonBackScatter,0,12,3);
   f2->SetParameter(0,4.8e-9);
   f2->SetParameter(1,0);
   f2->SetParameter(2,12);
   f2->Draw();
   c2->Update();
   return 0;
}
