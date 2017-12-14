//
// Qpolar.C
//
// Calculates the quantum amplitude for the scattering of a polarized 
// virtual photon from a polarized quark.  The final state consists of
// only the scattered quark.  The photon energy is adjusted to keep the
// final quark on-shell. The mass of the quark is adjustable.
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

Complex_t Qpolar(Double_t *var, Double_t *par)
{
   LDouble_t thetai=var[0];  // angle of incident quark wrt the photon axis (r)
   LDouble_t pin=par[0];     // momentum of the incident quark [GeV/c]
   LDouble_t Q2=par[1];      // Q^2 of the virtual photon [GeV^2]
   LDouble_t sigeIn=par[2];  // polarization state of initial quark

   const LDouble_t mQuark(0);

   TLepton qIn(mQuark), qFi(mQuark);
   TPhoton gIn;

   // Solve for the rest of the kinematics
   LDouble_t Ein=sqrt(pin*pin+mQuark*mQuark);
   TFourVectorReal pIn(Ein,pin*sin(thetai),0,-pin*cos(thetai));
   LDouble_t A=pIn[0]*pIn[0]-pIn[3]*pIn[3];
   LDouble_t B=-Q2*pIn[3];
   LDouble_t C=-Q2*(pIn[0]*pIn[0]+Q2/4);
   LDouble_t D=B*B-4*A*C;
   LDouble_t qin1=(D > 0)? (-B-sqrt(D))/(2*A) : 0;
   LDouble_t qin2=(D > 0)? (-B+sqrt(D))/(2*A) : 0;
   TFourVectorReal phIn(0,0,0,0);
   if ((qin1 > 1e-6) && (qin1*qin1 > Q2)) {
     phIn = TFourVectorReal(sqrt(qin1*qin1-Q2),0,0,qin1);
   }
   else if ((qin2 > 1e-6) && (qin2*qin2 > Q2)) {
     phIn = TFourVectorReal(sqrt(qin2*qin2-Q2),0,0,qin2);
   }
   TFourVectorReal pFi(pIn+phIn);
   pFi[0] = sqrt(mQuark*mQuark+pFi[1]*pFi[1]+pFi[2]*pFi[2]+pFi[3]*pFi[3]);
   qIn.SetMom(pIn);
   qFi.SetMom(pFi);
   gIn.SetMom(pFi-pIn);

   // Set the initial,final polarizations
   switch ((int)sigeIn) {
    case -1:
      qIn.SetPol(TThreeVectorReal(0,0,-1));
      qFi.SetPol(TThreeVectorReal(0,0,-1));
      break;
    case  0:
      qIn.SetPol(TThreeVectorReal(0,0,0));
      qFi.SetPol(TThreeVectorReal(0,0,0));
      break;
    case +1:
      qIn.SetPol(TThreeVectorReal(0,0,1));
      qFi.SetPol(TThreeVectorReal(0,0,1));
      break;
    default:
      cout << "Qpolar.C :"
           << "bad quark helicity specified, default to zero!"
           << endl;
      qIn.SetPol(TThreeVectorReal(0,0,0));
      qFi.SetPol(TThreeVectorReal(0,0,0));
   }

   // Compute the current amplitude
// qIn.Print();
// qFi.Print();
// gIn.Print();
   TDiracMatrix qmat1,qmat2;
   qmat1.SetUUbar(qIn.Mom(),qIn.SDM());
   qmat2.SetUUbar(qFi.Mom(),qFi.SDM());
   TDiracMatrix Aslash,AslashStar;
   AslashStar.Slash(gIn.EpsStar(1));
   Aslash.Slash(gIn.Eps(1));
   TDiracMatrix cur;
   cur = AslashStar*qmat2*Aslash*qmat1;
// qmat1.Print();
// qmat2.Print();
// Aslash.Print();
// AslashStar.Print();
// cur.Print();

// LDouble_t ans;
// cout << "continue? ";
// cin >> ans;

   return cur.Trace();
}

Double_t reQpolar(Double_t *var, Double_t *par)
{
   Complex_t amp = Qpolar(var,par);
   return amp.real();
}

Double_t imQpolar(Double_t *var, Double_t *par)
{
   Complex_t amp = Qpolar(var,par);
   return amp.imag();
}

Double_t absQpolar(Double_t *var, Double_t *par)
{
   Complex_t amp = Qpolar(var,par);
   return abs(amp);
}

Int_t demoQpolar()
{
   TCanvas *c1 = new TCanvas("c1","quark scattering amplitude",200,10,700,500);
   TF1 *f1 = new TF1("f1",absQpolar,0,2,3);
   f1->SetParameter(0,1.0);
   f1->SetParameter(1,4.0);
   f1->SetParameter(2,+1);
   f1->Draw();
   c1->Update();
   TCanvas *c2 = new TCanvas("c2","quark scattering amplitude",200,10,700,500);
   TF1 *f2 = new TF1("f2",reQpolar,0,2,3);
   f2->SetParameter(0,1.0);
   f2->SetParameter(1,4.0);
   f2->SetParameter(2,+1);
   f2->SetLineColor(2);
   f2->Draw();
   c2->Update();
   TCanvas *c3 = new TCanvas("c3","quark scattering amplitude",200,10,700,500);
   TF1 *f3 = new TF1("f3",imQpolar,0,2,3);
   f3->SetParameter(0,1.0);
   f3->SetParameter(1,4.0);
   f3->SetParameter(2,+1);
   f3->SetLineColor(3);
   f3->Draw();
   c3->Update();
   TCanvas *c4 = new TCanvas("c4","quark scattering amplitude",200,10,700,500);
   TF1 *f4 = new TF1("f4",absQpolar,0,2,3);
   f4->SetParameter(0,1.0);
   f4->SetParameter(1,4.0);
   f4->SetParameter(2,-1);
   f4->Draw();
   c4->Update();
   TCanvas *c5 = new TCanvas("c5","quark scattering amplitude",200,10,700,500);
   TF1 *f5 = new TF1("f5",reQpolar,0,2,3);
   f5->SetParameter(0,1.0);
   f5->SetParameter(1,4.0);
   f5->SetParameter(2,-1);
   f5->SetLineColor(2);
   f5->Draw();
   c5->Update();
   TCanvas *c6 = new TCanvas("c6","quark scattering amplitude",200,10,700,500);
   TF1 *f6 = new TF1("f6",imQpolar,0,2,3);
   f6->SetParameter(0,1.0);
   f6->SetParameter(1,4.0);
   f6->SetParameter(2,-1);
   f6->SetLineColor(3);
   f6->Draw();
   c6->Update();
   return 0;
}
