//
//  Test.C
//
// Implements a suite of simple tests to verify that the basic classes
// of the Dirac++ package have been built successfully.
//
// author: richard.t.jones at uconn.edu
// version: january 1, 2000

#include <iostream>

#include "Complex.h"
#include "TFourVectorReal.h"
#include "TFourVectorComplex.h"
#include "TThreeRotation.h"
#include "TLorentzBoost.h"
#include "TDiracSpinor.h"
#include "TDiracMatrix.h"

int tests()
{
   TThreeVectorReal v3v;
   v3v.Zero();
   std::cout << "tests.C is now loaded." << std::endl;
   std::cout << "To run a test of a class, type Test<class-name>()" << std::endl;
   std::cout << "Example: TestFourVectorReal()" << std::endl;
   return 0;
}

void TestThreeVectorReal()
{
   TThreeVectorReal v3v;
   v3v.Zero();
   std::cout << "Zero(): "; v3v.Print();

   v3v[1] = 1; v3v[2] = 2; v3v[3] = 3;
   std::cout << "Set to 1,2,3: "; v3v.Print();

   TThreeVectorReal v3w(4,5,6);
   std::cout << "Initialized to 4,5,6: "; v3w.Print();

   Float_t aarray[]={7,8,9,10};
   TThreeVectorReal v3x(aarray);
   std::cout << "Initialized to (Float_t *) 7,8,9: "; v3x.Print();

   LDouble_t darray[]={10,11,12,13};
   TThreeVectorReal v3y(darray);
   std::cout << "Initialized to (LDouble_t *) 10,11,12: "; v3y.Print();

   TThreeVectorReal v3z(v3v);
   std::cout << "Initialized to (TThreeVectorReal) 1,2,3: "; v3z.Print();

   v3z.SpaceInv();
   std::cout << "SpaceInv(): "; v3z.Print();

   v3z.Normalize(1);
   std::cout << "Normalize(1): "; v3z.Print();

   v3z.SetPolar(10,2.5,-1.2);
   std::cout << "SetPolar(10,2.5,-1.2) : GetPolar(";
   std::cout << v3z.Length() << "," << v3z.Theta() << "," << v3z.Phi() << ")" << std::endl;

   std::cout << "z == z ? " << ((v3z == v3z) ? "yes!" : "no!") << std::endl;
   std::cout << "z != w ? " << ((v3z != v3w) ? "yes!" : "no!") << std::endl;

   v3z.Rotate(-1.2,2.5,0);
   LDouble_t v3zhat[]={0,0,10};
   std::cout << "Rotate(phi,theta,0) gets back z axis? ";
   std::cout << ((v3z == TThreeVectorReal(v3zhat)) ? "yes!" : "no!" ) << std::endl;

   std::cout << "(w cross x) != 0 ? ";
   std::cout << ((v3z.Cross(v3w,v3x).Length() >= v3w.Resolution()) ? "yes!" : "no!" ) << std::endl;

   std::cout << "w dot (w cross x) == 0 ? ";
   std::cout << ((v3w.Dot(v3z.Cross(v3w,v3x)) < v3w.Resolution()) ? "yes!" : "no!" ) << std::endl;

   v3y.Cross(v3w,v3w-v3x);
   v3z.Cross(v3x,v3w);
   std::cout << "w cross (w - x) == x cross w ? ";
   std::cout << ((v3y.Cross(v3w,v3w-v3x) == v3z.Cross(v3x,v3w)) ? "yes!" : "no!" ) << std::endl;
}

void TestThreeVectorComplex()
{
   TThreeVectorComplex v3v;
   v3v.Zero();
   std::cout << "Zero(): "; v3v.Print();

   v3v[1] = 1; v3v[2] = 2; v3v[3] = 3;
   std::cout << "Set to 1,2,3: "; v3v.Print();

   TThreeVectorComplex v3w(4,5,6);
   std::cout << "Initialized to 4,5,6: "; v3w.Print();

   Float_t farray[]={7,8,9,10};
   TThreeVectorComplex v3x(farray);
   std::cout << "Initialized to (Float_t *) 7,8,9: "; v3x.Print();

   LDouble_t darray[]={10,11,12,13};
   TThreeVectorComplex v3y(darray);
   std::cout << "Initialized to (LDouble_t *) 10,11,12: "; v3y.Print();

   TThreeVectorComplex v3z(v3v);
   std::cout << "Initialized to (TThreeVectorComplex) 1,2,3: "; v3z.Print();

   v3z.SpaceInv();
   std::cout << "SpaceInv(): "; v3z.Print();

   v3z.Normalize(1);
   std::cout << "Normalize(1): "; v3z.Print();

   std::cout << "z == z ? " << ((v3z == v3z) ? "yes!" : "no!") << std::endl;
   std::cout << "z != w ? " << ((v3z != v3w) ? "yes!" : "no!") << std::endl;

   TThreeVectorReal v3r(242.99,-485,0.6666);
   v3z = v3r.Normalize(10);
   v3z.Rotate(v3r.Phi(),v3r.Theta(),0);
   LDouble_t v3zhat[]={0,0,10};
   std::cout << "Rotate(phi,theta,0) gets back z axis? ";
   std::cout << ((v3z == TThreeVectorReal(v3zhat)) ? "yes!" : "no!" ) << std::endl;

   std::cout << "(w cross x) != 0 ? ";
   std::cout << ((v3z.Cross(v3w,v3x).Length() >= v3w.Resolution()) ? "yes!" : "no!" ) << std::endl;

   std::cout << "w dot (w cross x) == 0 ? ";
   std::cout << ((abs(v3w.Dot(v3z.Cross(v3w,v3x))) < v3w.Resolution()) ? "yes!" : "no!" ) << std::endl;

   v3y.Cross(v3w,v3w-v3x);
   v3z.Cross(v3x,v3w);
   std::cout << "w cross (w - x) == x cross w ? ";
   std::cout << ((v3y.Cross(v3w,v3w-v3x) == v3z.Cross(v3x,v3w)) ? "yes!" : "no!" ) << std::endl;
}

void TestFourVectorReal()
{
   TFourVectorReal v4v;
   v4v.Zero();
   std::cout << "Zero(): "; v4v.Print();

   v4v[0] = 10; v4v[1] = 1; v4v[2] = 2; v4v[3] = 3;
   std::cout << "Set to 10,1,2,3: "; v4v.Print();

   TFourVectorReal v4w(-4,4,5,6);
   std::cout << "Initialized to -4,4,5,6: "; v4w.Print();

   Float_t aarray[]={7,8,9,10};
   TFourVectorReal v4x(aarray);
   std::cout << "Initialized to (Float_t *) 7,8,9,10: "; v4x.Print();

   LDouble_t darray[]={10,11,12,13};
   TFourVectorReal v4y(darray);
   std::cout << "Initialized to (LDouble_t *) 10,11,12,13: "; v4y.Print();

   TFourVectorReal v4z(v4v);
   std::cout << "Initialized to (TFourVectorReal) 10,1,2,3: "; v4z.Print();

   v4z.SpaceInv();
   std::cout << "SpaceInv(): "; v4z.Print();

   v4z /= v4z.Invariant();
   std::cout << "v4z /= v4z.Invariant(): "; v4z.Print();

   TFourVectorReal v4u(v4z);
   v4z.Boost(0,0,.99);
   std::cout << "v4z.Boost(0,0,.99): "; v4z.Print();

   std::cout << "Lorentz invariant? ";
   std::cout << ((abs(v4z.Invariant()-v4u.Invariant()) < v4z.Resolution()) ?
           "yes!" : "no!") << std::endl;

   v4z.Boost(0,0,-.99);
   std::cout << "Boost back, recover original vector? ";
   std::cout << ((v4z == v4u) ? "yes!" : "no!") << std::endl;
}

void TestFourVectorComplex()
{
   TFourVectorComplex v4v;
   v4v.Zero();
   std::cout << "Zero(): "; v4v.Print();

   v4v[0] = 10; v4v[1] = 1; v4v[2] = 2; v4v[3] = 3;
   std::cout << "Set to 10,1,2,3: "; v4v.Print();

   TFourVectorComplex v4w(-4,4,5,6);
   std::cout << "Initialized to -4,4,5,6: "; v4w.Print();

   Float_t aarray[]={7,8,9,10};
   TFourVectorComplex v4x(aarray);
   std::cout << "Initialized to (Float_t *) 7,8,9,10: "; v4x.Print();

   LDouble_t darray[]={10,11,12,13};
   TFourVectorComplex v4y(darray);
   std::cout << "Initialized to (LDouble_t *) 10,11,12,13: "; v4y.Print();

   TFourVectorComplex v4z(v4v);
   std::cout << "Initialized to (TFourVectorComplex) 10,1,2,3: "; v4z.Print();

   v4z.SpaceInv();
   std::cout << "SpaceInv(): "; v4z.Print();

   v4z /= v4z.Invariant();
   std::cout << "v4z /= v4z.Invariant(): "; v4z.Print();

   TFourVectorComplex v4u(v4z);
   v4z.Boost(0,0,.99);
   std::cout << "v4z.Boost(0,0,.99): "; v4z.Print();

   std::cout << "Lorentz invariant? ";
   std::cout << ((abs(v4z.Invariant()-v4u.Invariant()) < v4z.Resolution()) ?
           "yes!" : "no!") << std::endl;

   v4z.Boost(0,0,-.99);
   std::cout << "Boost back, recover original vector? ";
   std::cout << ((v4z == v4u) ? "yes!" : "no!") << std::endl;
}

void TestRotation()
{
   TThreeVectorReal w1(0,0,1.2);
   TThreeVectorReal w2(0,0.887,0);
   TThreeVectorReal w3(0,0,-1.557);
   TThreeRotation r1(w1), r2(w2), r3(w3);
   TThreeRotation r123=r3*r2*r1;
   LDouble_t phi=0,theta=0,psi=0;
   r123.GetEuler(phi,theta,psi);
   std::cout << "First rotate 1.2 radians about z axis" <<std::endl;
   std::cout << "Then rotate 0.887 radians about y' axis" <<std::endl;
   std::cout << "Last rotate -1.557 radians about z'' axis" <<std::endl;
   std::cout << "GetEuler returns phi=" << phi
        << ", theta=" << theta << ", psi=" << psi << std::endl;

   TThreeRotation r4;
   r123 *= r4.SetEuler(1.557,-0.887,-1.2);
   std::cout << "Take the product with SetEuler(1.557,-0.887,-1.2)" << std::endl;
   std::cout << "Do we get back the null rotation? "
        << (r123.IsNull() ? "yes!" : "no!") << std::endl;

   r123 = r3*r2*r1;
   r123.Invert();
   std::cout << "Does Invert() after SetEuler(a,b,c) give SetEuler(-c,-b,-a)? "
        << (r123 == r4 ? "yes!" : "no!" ) << std::endl;

   std::cout << "Is a product of rotations a rotation? "
        << (r123.IsRotation() ? "yes!" : "no!") << std::endl;
   std::cout << "Is a rotation matrix orthogonal? "
        << (r123.IsOrthogonal() ? "yes!" : "no!") << std::endl;
   std::cout << "Is a rotation distinct from a boost? "
        << (r123.IsLorentzBoost() ? "no!" : "yes!") << std::endl;
}

void TestLorentzBoost()
{
   TLorentzBoost b1(0,0,0.6);
   TLorentzBoost b2(-.66,0.2,0.45);
   TLorentzBoost b3(0,0,-0.6);
   TLorentzTransform t13=b3*b1;
   std::cout << "Do colinear boosts commute? "
        << (t13.IsNull() ? "yes!" : "no!" ) << std::endl;

   TLorentzTransform t123=b3*b2*b1;
   std::cout << "But not if they are not colinear? "
        << ((t123!=b2) ? "yes!" : "no!" ) << std::endl;

   TLorentzBoost b4(.776,-.433,-0.06);
   TThreeVectorReal beta4 = b4.Beta();
   TLorentzBoost b4inv(-beta4);
   TLorentzTransform b5=b4*b4inv;
   std::cout << "Is a boost(beta) followed by boost(-beta) a null operation? "
        << (b5.IsNull() ? "yes!" : "no!") << std::endl;

   TFourVectorReal p(1,0,0,0);
   TLorentzTransform t12=b2*b1;
   TLorentzBoost b12,b123;
   TThreeRotation r12,r123;
   t12.Factorize(b12,r12);
   TThreeVectorReal q = b12.Beta();
   b3.SetBeta(q.SpaceInv());
   t123=b3*b12;
   std::cout << "Do B3*B2*B1 amount to a rotation if beta(B3) = -beta(B2*B1)? "
        << (t123.IsRotation() ? "yes!" : "no!" ) << std::endl;
}

void TestPauliSpinor()
{
   TUnitVector phat(28,-37,2.8);
   TUnitVector qhat(-4,.223,9);
   TPauliSpinor s1(phat),s2(qhat);
   TUnitVector pol=s1.Polar();
   std::cout << "Does Polar() give back what SetPolar put in? "
        << ((s1.Polar()==phat.Normalize(1)) ? "yes!" : "no!") << std::endl;

   TThreeRotation r1(.334,-21.554,1.28);
   TThreeVectorReal phatPrime=r1*phat;
   Complex_t dot12=s1.ScalarProd(s2);
   s1.Rotate(r1);
   s2.Rotate(r1);
   TThreeVectorReal polPrime = s1.Polar();
   std::cout << "Does polarization rotate like a three-vector? "
        << ((s1.Polar()==phatPrime) ? "yes!" : "no!") << std::endl;

   std::cout << "Is the scalar product invariant under rotations? "
        << ((dot12 == s1.ScalarProd(s2)) ? "yes!" : "no!") << std::endl;
}

void TestPauliMatrix()
{
   TThreeVectorReal axis(.2,1.11,-.44);
   TThreeRotation r1(axis),r2(10,2.33,-1),r3;
   r3.SetAxis(-axis,.333);
   TThreeRotation r123=r3*r2*r1;
   TPauliMatrix p1,p2,p3;
   p1.SetRotation(r1);
   p2.SetRotation(r2);
   p3.SetRotation(r3);
   TThreeVectorComplex b;
   Complex_t a(0);
   std::cout << "Is Pauli rotation matrix unitary, not hermetian? "
        << (p1.IsUnitary() ? "yes!" : "no!") << " and "
        << (p1.IsHermetian() ? "no!" : "yes!") << std::endl;

   p1.Decompose(a,b);
   b.Normalize(1);
   axis.Normalize(1);
   std::cout << "Does Decompose give back the rotation vector? "
        << ((b.ImagPart()==axis) ? "yes!" : "no!" ) << std::endl;
   TPauliMatrix p123=p3*p2*p1;
   p123.Decompose(a,b);
   b.Normalize(1);
   axis = r123.Axis();
   TThreeVectorReal ahat(axis);
   ahat.Normalize(1);
   std::cout << "Does the product of rotation equal the rotation product? "
        << ((b.ImagPart()==ahat) ? "yes" : "no")
        << " and "
        << ((abs(2.0L*acos(a)-axis.Length())<axis.Resolution())
           ? "yes!" : "no!") << std::endl;
}

void TestDiracSpinor()
{
   TFourVectorReal p(13,2.33,-5.8,-1.002);
   TDiracSpinor d1(p,-1);
   TDiracMatrix pslash;
   pslash.Slash(p);
   TDiracSpinor lhs=pslash*d1;
   TDiracSpinor rhs=p.Invariant()*d1;
   std::cout << "Does TDiracSpinor(p) satisfy the Dirac equation? "
        << ((lhs==rhs) ? "yes!" : "no!" ) << std::endl;

   TThreeVectorReal polar=d1.Upper().Polar();
   TThreeVectorReal phat(p);
   std::cout << "Does TDiracSpinor(p,-1) give correct polarization? "
        << ((polar==-phat.Normalize(1)) ? "yes!" : "no!") << std::endl;

   TThreeRotation r1(5,-2.34,.567);
   TDiracSpinor d1Prime(d1);
   TDiracMatrix dm1(r1);
   d1Prime.Operate(dm1);
   TThreeVectorReal polarPrime=d1Prime.Upper().Polar();
   TFourVectorReal pPrime=(TLorentzTransform)r1*p;
   pPrime.Normalize(1);
   std::cout << "Does the spinor polarization rotate like a vector? "
        << (((TThreeVectorReal)pPrime==-polarPrime) ? "yes!" : "no!") << std::endl;

   TLorentzBoost b1(.33,.44,-.55);
   dm1.SetBoost(b1);
   d1Prime=dm1*d1;
   pPrime=b1*p;
   pslash.Slash(pPrime);
   TDiracSpinor lhsPrime=pslash*d1Prime;
   TDiracSpinor rhsPrime=p.Invariant()*d1Prime;
   std::cout << "Does boosted spinor satisfy the Dirac equation with same m? "
        << ((lhsPrime==rhsPrime) ? "yes!" : "no!" ) << std::endl;

   std::cout << "Normalization still sqrt(2*E)? "
        << ((abs(d1Prime.NormSqr()-2*pPrime[0])<p.Resolution()) ?
             "yes!" : "no!") << std::endl;
}

void TestDiracMatrix()
{
   TDiracMatrix dm1(1.0L);
   std::cout << "Does TDiracMatrix(1.0) give identity matrix? "
        << (dm1.IsIdentity() ? "yes! " : " no!") << std::endl;

   TFourVectorReal p(-3.48,2.26,1.56,-0.96);
   dm1.SetUUbar(p);
   TDiracMatrix dm2;
   dm2.Slash(p);
   LDouble_t m=p.Invariant();
   dm2 += m;
   std::cout << "Is the UUbar-projector equal to pSlash+m ? "
        << (dm1==dm2 ? "yes!" : "no!") << std::endl;

   TDiracMatrix sigma3(kDiracSigma3);
   TDiracMatrix gamma1(kDiracGamma1),gamma2(kDiracGamma2);
   Complex_t i_(0,1);
   TDiracMatrix commutator12=(i_/2.L)*(gamma1*gamma2-gamma2*gamma1);
   std::cout << "Is Sigma_3 equal to (i/2)[gamma_1,gamma_2] ? "
        << (commutator12==sigma3 ? "yes!" : "no!") << std::endl;

   Complex_t minusOne(-1,0);
   std::cout << "Component of Sigma_3 belonging to commutator [2,1] = -1 ? "
        << ((sigma3.Component(kDiracGamma2,kDiracGamma1) == minusOne) ?
          "yes!" : "no!") << std::endl;
}
