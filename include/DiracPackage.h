//
// Dirac Spinor Algebra Package
//
// The present package implements all the basic algorithms dealing
// with Dirac spinors, which form a fundamental representation of the
// group SL(2,2).  The basic classes are DiracSpinor and DiracMatrix,
// which are 4x1 and 4x4 complex matrices, respectively. In this
// context any complex 4x4 matrix that operates on Dirac spinors is
// called a Dirac matrix, and not simply the four or five Dirac gamma
// matrices.  The standard representation of Dirac is used for the
// gamma matrices 0-5.  The generators of the Lorentz group are the
// Sigma (rotation generators) and Kappa (boost generators) matrices.
//
// The standard matrices are identified by a discrete index of enum
// type EDiracIndex.  A EDiracIndex can take on a value from the list
//    kDiracOne,    kDiracGamma1,   kDiracGamma2,   kDiracGamma3,
//    kDiracGamma4, kDiracGamma5,   kDiracSigma1,   kDiracSigma2,
//    kDiracSigma3, kDiracKappa1,   kDiracKappa2,   kDiracKappa3.
// The constructor invoked with two EDiracIndex values i,j returns
// i_/2 [TDiracMatrix(i),TDiracMatrix(j)] where [a,b] denotes the com-
// utator of matrices a and b, and i_ is the positive square root of
// -1.  In general Dirac matrices describe operators and Dirac spinors
// describe relativistic fermion states.  Dirac matrices are also used
// to describe mixed states, ensembles that contain mixtures of
// particles described by more than one Dirac spinor.
//
// Spinors and matrices can be transformed under rotations and boosts
// according to the commutation rules for the group.  The most general
// transformation combining rotations and boosts is described by the
// LorentzTransform group defined in TFourVector.h.  All angles are
// assumed to be in radians.
//
// This package depends on the ROOT framework (http://root.cern.ch).
//
// author: richard.t.jones at uconn.edu
// version: january 1, 2000
