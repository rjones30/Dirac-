#include "Double.h"

#pragma link C++ class TThreeVectorReal-;

#pragma link C++ function operator+(const TThreeVectorReal&,const TThreeVectorReal&);
#pragma link C++ function operator-(const TThreeVectorReal&,const TThreeVectorReal&);
#pragma link C++ function operator*(const TThreeVectorReal&,const LDouble_t);
#pragma link C++ function operator*(const LDouble_t,const TThreeVectorReal&);
#pragma link C++ function operator/(const TThreeVectorReal&,const LDouble_t);
#pragma link C++ function operator>>(TBuffer&,TThreeVectorReal*&);
#pragma link C++ function operator<<(TBuffer&,const TThreeVectorReal*);
