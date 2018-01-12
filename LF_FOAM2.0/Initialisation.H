dimensionedScalar one = dimensionedScalar("one", dimensionSet(0,1,0,0,0,0,0), 1.0);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
	IOobject
	(
	    "transportProperties",
	     runTime.constant(),
	     mesh,
	     IOobject::MUST_READ,
	      IOobject::NO_WRITE
	 )
);

dictionary LevelSetParameter
(
	transportProperties.subDict("LevelSetParameter")
);

scalar EPS = readScalar(LevelSetParameter.lookup("interface_thickness"));

scalar CFL = readScalar(LevelSetParameter.lookup("cfl_number"));

word phaseName(wordList(transportProperties.lookup("phases"))[0]);

volScalarField alpha
(
       IOobject
       (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       ),
       mesh
);

/*volScalarField H
(
      IOobject
      (
          "H",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("H",dimensionSet(0,0,0,0,0,0,0),scalar(0))
);

volScalarField dirac
(
      IOobject
      (
          "dirac",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("diracf",dimensionSet(0,0,0,0,0,0,0),scalar(0))
);

surfaceScalarField theta
(
      IOobject
      (
          "theta",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("theta",dimensionSet(0,0,0,0,0,0,0),scalar(0))
);
*/

const surfaceScalarField deltaCoeff = alpha.mesh().deltaCoeffs();
const scalar DX = (1.0/min(deltaCoeff).value());
  
scalar Dtau  = CFL*DX;
scalar Gamma = 0.75*DX;
scalar NITER = EPS/Dtau;
scalar ERR (1.0);

volScalarField d
(
      IOobject
      (
          "d",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("d",dimensionSet(0,0,0,0,0,0,0),scalar(0)),
      alpha.boundaryField().types()
);

d  == Gamma*(2*alpha-1.0);

volScalarField d0
(
      IOobject
      (
          "d0",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("d0",dimensionSet(0,0,0,0,0,0,0),scalar(0)),
      d.boundaryField().types()
);

d0 == d;

for (int loop = 0; loop < NITER; loop++)
{

  // calculate distance function field
  
/*  fvScalarMatrix dEqn
  (
        fvm::Sp(scalar(1),d) + Dtau*sign(d0)*mag(one*fvc::grad(d)) == d + Dtau*sign(d0)
  );
  dEqn.solve();
*/
    Info<<"loop ="<<loop<<tab <<"min(d) = "<<min(d).value() <<tab<<"max(d)="<<max(d).value()<<endl;
  d += Dtau*sign(d0)*(scalar(1)-mag(one*fvc::grad(d)));

}

ERR = CALC_ERR(d, mag(fvc::grad(d)), EPS);

Info << "NITER=" << NITER << tab << "ERR=" << ERR << endl;
//Info << "ERR=" <<ERR<<endl;

Info
<< "  Min(d) = " << min(d).value()
<< "  Max(d) = " << max(d).value()
<< endl;

//#include "Reinitialisation.H"

//H == CALC_HEAVISIDE(d,EPS);
//dirac == CALC_DIRAC(d,EPS);

// Construct interface from U phi and d distribution
LSimmiscibleIncompressibleTwoPhaseMixture mixture(U, phi,d);
