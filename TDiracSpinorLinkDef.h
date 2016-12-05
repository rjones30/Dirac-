#pragma link C++ class TDiracSpinor-;

#pragma link C++ function operator==(const Complex_t*,const TDiracSpinor&);
#pragma link C++ function operator!=(const Complex_t*,const TDiracSpinor&);
#pragma link C++ function operator+(const TDiracSpinor&,const TDiracSpinor&);
#pragma link C++ function operator+(const TDiracSpinor&,const Complex_t*);
#pragma link C++ function operator+(const Complex_t*,const TDiracSpinor&);
#pragma link C++ function operator-(const TDiracSpinor&,const TDiracSpinor&);
#pragma link C++ function operator-(const TDiracSpinor&,const Complex_t*);
#pragma link C++ function operator-(const Complex_t*,const TDiracSpinor&);
#pragma link C++ function operator*(const TDiracSpinor&,const Double_t&);
#pragma link C++ function operator*(const Double_t&,const TDiracSpinor&);
#pragma link C++ function operator*(const TDiracSpinor&,const Complex_t&);
#pragma link C++ function operator*(const Complex_t&,const TDiracSpinor&);
#pragma link C++ function operator*(const TDiracMatrix&,const TDiracSpinor&);
#pragma link C++ function operator/(const TDiracSpinor&,const Complex_t&);
#pragma link C++ function operator>>(TBuffer&,TDiracSpinor*&);
#pragma link C++ function operator<<(TBuffer&,const TDiracSpinor*);
