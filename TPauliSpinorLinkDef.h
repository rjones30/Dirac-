#pragma link C++ class TPauliSpinor-;

#pragma link C++ function operator==(const Complex_t*,const TPauliSpinor&);
#pragma link C++ function operator!=(const Complex_t*,const TPauliSpinor&);
#pragma link C++ function operator+(const TPauliSpinor&,const TPauliSpinor&);
#pragma link C++ function operator+(const TPauliSpinor&,const Complex_t*);
#pragma link C++ function operator+(const Complex_t*,const TPauliSpinor&);
#pragma link C++ function operator-(const TPauliSpinor&,const TPauliSpinor&);
#pragma link C++ function operator-(const TPauliSpinor&,const Complex_t*);
#pragma link C++ function operator-(const Complex_t*,const TPauliSpinor&);
#pragma link C++ function operator*(const TPauliSpinor&,const LDouble_t&);
#pragma link C++ function operator*(const LDouble_t&,const TPauliSpinor&);
#pragma link C++ function operator*(const TPauliSpinor&,const Complex_t&);
#pragma link C++ function operator*(const Complex_t&,const TPauliSpinor&);
#pragma link C++ function operator*(const TPauliMatrix&,const TPauliSpinor&);
#pragma link C++ function operator/(const TPauliSpinor&,const Complex_t&);
#pragma link C++ function operator>>(TBuffer&,TPauliSpinor*&);
#pragma link C++ function operator<<(TBuffer&,const TPauliSpinor*);
