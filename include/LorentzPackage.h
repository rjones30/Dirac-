//
// Lorentz Algebra Package
//
// The present package implements all the basic algorithms dealing
// with three-vectors and four-vectors, together with their transform-
// ations.  Four-vectors are derived from three-vectors and inherit
// all of their members.  Direct access to the components is provided
// through the subscript operator [i] which covers the range 1...3 for
// three-vectors and 0...3 for four-vectors.  Transformations are
// implemented as a friend class so that they can operate directly on
// the data members of the vector, which are otherwise hidden.  The
// special transformations (rotations and boosts) inherit from the
// general class LorentzTransform.  Products of rotations are other
// rotations, whereas the product of a boost with anything is simply
// a LorentzTransform.  The LorentzTransform objects can be tested
// for the property of being a pure rotation or boost.  They can also
// implement non-isochronous and improper transformations.
//
// Rotations may be specified either by Euler angles or by a rotation
// axis.  All angles are assumed to be in radians.  Vector classes are
// defined for both Double_t and Complex_t generic types.  For complex
// vectors there are several additional member functions to deal with
// operations that are specific to complex numbers.
//
// The classes comprising this package are:
//   TThreeVectorReal is a base class
//   TThreeVectorComplex is a base class
//   TFourVectorReal is a TThreeVectorReal
//   TFourVectorComplex is a TThreeVectorComplex
//   TLorentzTransform is a base class
//   TThreeRotation is a TLorentzTransform
//   TLorentzBoost is a TLorentzTransform
// The following aliases are defined for these classes:
//   TUnitVector is an alias for TThreeVectorReal
//   TThreeVector is an alias for TThreeVectorReal
//   TFourVector is an alias for TFourVectorReal
//
// This package depends on the ROOT framework (http://root.cern.ch).
//
// author: richard.t.jones at uconn.edu
// version: january 1, 2000
