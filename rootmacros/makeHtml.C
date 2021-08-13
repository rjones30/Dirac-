#include <THtml.h>

void makeHtml()
{
   THtml docs;
   docs.MakeClass("TThreeVectorReal");
   docs.MakeClass("TThreeVectorComplex");
   docs.MakeClass("TFourVectorReal");
   docs.MakeClass("TFourVectorComplex");
   docs.MakeClass("TLorentzTransform");
   docs.MakeClass("TLorentzBoost");
   docs.MakeClass("TThreeRotation");
   docs.MakeClass("TPauliSpinor");
   docs.MakeClass("TPauliMatrix");
   docs.MakeClass("TDiracSpinor");
   docs.MakeClass("TDiracMatrix");
   docs.MakeClass("TPhoton");
   docs.MakeClass("TLepton");
   docs.MakeClass("TCrossSection");
   docs.MakeClass("TDiracMatrix");
}
