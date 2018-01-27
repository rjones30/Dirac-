//
// Brems.C
//
// Calculates the bremsstrahlung production rate including polarization
// effects for a given kinematics.  The kinematics are set by the value
// of the recoil momentum and the momentum of the final photon.  The
// results are reported as probability per GeV photon energy per beam
// electron (Brems) or as linear polarization (BremsPolarization).
//
// author: richard.t.jones at uconn.edu
// version: january 1, 2000

#include "Complex.h"
#include "TPhoton.h"
#include "TLepton.h"
#include "TCrossSection.h"
#include "constants.h"
#include "sqr.h"

#include <TCanvas.h>
#include <TF1.h>

Double_t Brems(Double_t *var, Double_t *par)
{
   LDouble_t kout=var[0];
   TThreeVectorReal qRecoil(9.83425e-6,0.,par[0]);
   LDouble_t phi=par[1];
   LDouble_t pin=par[2];

   TLepton eIn(mElectron), eOut(mElectron);
   TPhoton gOut;

   // Solve for the rest of the kinematics
   TThreeVectorReal p;
   LDouble_t pout = pin-kout;
   eIn.SetMom(TThreeVectorReal(0,0,pin));
   LDouble_t A = kout * (pin - qRecoil[3]);
   LDouble_t B = 2 * qRecoil[1] * kout * cos(phi);
   LDouble_t C = qRecoil.LengthSqr() + 
                 2 * kout * (sqrt(sqr(pin) + sqr(mElectron)) - pin) -
                 2 * qRecoil[3] * pout;
   LDouble_t discrim = B*B - 4*A*C;
   LDouble_t theta1, theta2;
   if (discrim >= 0) {
      theta1 = (-B - sqrt(discrim)) / (2*A);
      theta2 = (-B + sqrt(discrim)) / (2*A);
   }
   else {
      return 0;
   }
   LDouble_t theta = (theta1 > 0)? theta1 : theta2;
   if (theta <= 0)
      return 0;

   gOut.SetMom(p.SetPolar(kout,theta,phi));
   p = eIn.Mom()-qRecoil-p;
   eOut.SetMom(p);

   // Set the initial,final polarizations
   eIn.SetPol(TThreeVectorReal(0,0,0));
   gOut.AllPol();
   eOut.AllPol();

   // Multiply the basic cross section by the form factors
   LDouble_t result=TCrossSection::Bremsstrahlung(eIn,eOut,gOut);
   const LDouble_t Z=6;
   const LDouble_t Sff=8*Z;
   const LDouble_t Aphonon=0.5e9;                      // in /GeV**2
   const LDouble_t Gff=exp(-Aphonon*qRecoil.LengthSqr()/2);
   const LDouble_t beta=111*pow(Z,-1/3.)/mElectron;    // ff cutoff in /GeV
   const LDouble_t Fff=1/(1+sqr(beta)*qRecoil.LengthSqr());
   const LDouble_t hbarc=1.97327e-6;                   // in GeV.Angstroms
   const LDouble_t Vcell=45.5;                         // in Angstroms**3
   const LDouble_t XffSqr_d3q=pow(2*PI_*hbarc,3)/Vcell;
   result *= sqr(Sff) * sqr(Gff) * sqr(1-Fff) * XffSqr_d3q;

   // Multiply the cross section by target thickness
   result *= 1e-34;                    // from ub to m**2
   const LDouble_t t=20e-6;            // in m
   result *= t/(Vcell*1e-30);
result *= 2.2e-6/1.6e-19;
result *= 2*PI_;
   return result;
}

Int_t demoBrems(Double_t phi=0)
{
   TCanvas *c1 = new TCanvas("c1","Bremsstrahlung Production Rate",200,10,700,500);
   TF1 *f1 = new TF1("f1",Brems,7,10.0,3);
   Double_t params[3];
   params[0] = 33.e-9;//(0,55.1e-9);
   params[1] = phi;
   params[2] = 12;
   f1->SetParameter(0,params[0]);
   f1->SetParameter(1,params[1]);
   f1->SetParameter(2,params[2]);
   f1->Draw();
   c1->Update();
   std::cout << "Integral from 8.4 to 9.0 is " << f1->Integral(8.4,9.0)
             << std::endl;
   return 0;
}

Double_t BremsPolarization(Double_t *var, Double_t *par)
{
   TThreeVectorReal qRecoil(9.83425e-6,0.,par[0]);
   LDouble_t kout=par[1];
   LDouble_t phi=var[0];
   LDouble_t pin=par[2];

   TLepton eIn(mElectron), eOut(mElectron);
   TPhoton gOut;

   // Solve for the rest of the kinematics
   TThreeVectorReal p;
   LDouble_t pout = pin-kout;
   eIn.SetMom(TThreeVectorReal(0,0,pin));
   LDouble_t theta=0, thetaSqr;
   for (Int_t i=0; i<5; i++) {
      thetaSqr = (pout/kout)*(2*pin*qRecoil[3])/sqr(mElectron) -1
                -(pin/kout)*(qRecoil.LengthSqr())/sqr(mElectron)
                -(2*qRecoil[1]/mElectron)*theta*cos(phi);
      if (thetaSqr < 0) { return 0; }
      theta = sqrt(thetaSqr);
   }
   theta *= (mElectron/pin);
   gOut.SetMom(p.SetPolar(kout,theta,phi));
   p = eIn.Mom()-qRecoil-p;
   eOut.SetMom(p);

   // Measure the transverse polarization in the plane of qRecoil
   eIn.SetPol(TThreeVectorReal(0,0,1));
   eOut.AllPol();
   gOut.SetPol(TThreeVectorReal(1,0,0));
   LDouble_t Xrate=TCrossSection::Bremsstrahlung(eIn,eOut,gOut);
   gOut.SetPol(TThreeVectorReal(0,1,0));
   LDouble_t Yrate=TCrossSection::Bremsstrahlung(eIn,eOut,gOut);

   return (Xrate-Yrate)/(Xrate+Yrate);
}

Int_t demoBremsPolarization()
{
   TCanvas *c2 = new TCanvas("c1","Bremsstrahlung Polarization",200,10,700,500);
   TF1 *f2 = new TF1("f2",BremsPolarization,0,6.3,3);
   Double_t params[3];
   params[0] = 33.e-9;//(0,55.1e-9);
   params[1] = 8.70;
   params[2] = 12;
   f2->SetParameter(0,params[0]);
   f2->SetParameter(1,params[1]);
   f2->SetParameter(2,params[2]);
   f2->Draw();
   c2->Update();
   return 0;
}
