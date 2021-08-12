#pragma link C++ class TPauliMatrix-;

#pragma link C++ function operator+(const TPauliMatrix&,const TPauliMatrix&);
#pragma link C++ function operator+(const Complex_t&,const TPauliMatrix&);
#pragma link C++ function operator+(const TPauliMatrix&,const Complex_t&);
#pragma link C++ function operator-(const TPauliMatrix&,const TPauliMatrix&);
#pragma link C++ function operator-(const Complex_t&,const TPauliMatrix&);
#pragma link C++ function operator-(const TPauliMatrix&,const Complex_t&);
#pragma link C++ function operator*(const TPauliMatrix&,const TPauliMatrix&);
#pragma link C++ function operator*(const Complex_t&,const TPauliMatrix&);
#pragma link C++ function operator*(const TPauliMatrix&,const Complex_t&);
#pragma link C++ function operator/(const TPauliMatrix&,const TPauliMatrix&);
#pragma link C++ function operator/(const Complex_t&,const TPauliMatrix&);
#pragma link C++ function operator/(const TPauliMatrix&,const Complex_t&);
#pragma link C++ function operator>>(TBuffer&,TPauliMatrix*&);
#pragma link C++ function operator<<(TBuffer&,const TPauliMatrix*);

#pragma link C++ enum EPauliIndex;
