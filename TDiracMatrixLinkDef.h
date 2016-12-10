#pragma link C++ class TDiracMatrix-;

#pragma link C++ function operator+(const TDiracMatrix&,const TDiracMatrix&);
#pragma link C++ function operator+(const Complex_t&,const TDiracMatrix&);
#pragma link C++ function operator+(const TDiracMatrix&,const Complex_t&);
#pragma link C++ function operator+(const LDouble_t&,const TDiracMatrix&);
#pragma link C++ function operator+(const TDiracMatrix&,const LDouble_t&);
#pragma link C++ function operator-(const TDiracMatrix&,const TDiracMatrix&);
#pragma link C++ function operator-(const Complex_t&,const TDiracMatrix&);
#pragma link C++ function operator-(const TDiracMatrix&,const Complex_t&);
#pragma link C++ function operator-(const LDouble_t&,const TDiracMatrix&);
#pragma link C++ function operator-(const TDiracMatrix&,const LDouble_t&);
#pragma link C++ function operator*(const TDiracMatrix&,const TDiracMatrix&);
#pragma link C++ function operator*(const Complex_t&,const TDiracMatrix&);
#pragma link C++ function operator*(const TDiracMatrix&,const Complex_t&);
#pragma link C++ function operator*(const LDouble_t&,const TDiracMatrix&);
#pragma link C++ function operator*(const TDiracMatrix&,const LDouble_t&);
#pragma link C++ function operator/(const TDiracMatrix&,const TDiracMatrix&);
#pragma link C++ function operator/(const Complex_t&,const TDiracMatrix&);
#pragma link C++ function operator/(const TDiracMatrix&,const Complex_t&);
#pragma link C++ function operator/(const LDouble_t&,const TDiracMatrix&);
#pragma link C++ function operator/(const TDiracMatrix&,const LDouble_t&);
#pragma link C++ function operator>>(TBuffer&,TDiracMatrix*&);
#pragma link C++ function operator<<(TBuffer&,const TDiracMatrix*);

#pragma link C++ enum EDiracIndex;
