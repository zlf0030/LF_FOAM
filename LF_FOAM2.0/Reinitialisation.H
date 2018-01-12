const surfaceScalarField deltaCoeff1 = alpha.mesh().deltaCoeffs();
const scalar dX = (1.0/max(deltaCoeff1).value());
scalar Dtau  = CFL*dX;
scalar Gamma = 0.75*dX;
scalar EPS = NITER*Dtau;

d  == Gamma*(2*alpha1-1.0);

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
/*	
	fvScalarMatrix dEqn
	(
    		fvm::Sp(scalar(1),d) + Dtau*sign(d0)*mag(one*fvc::grad(d)) == d + Dtau*sign(d0)
	);
	dEqn.solve();
*/	
	d += Dtau*sign(d0)*(scalar(1)-mag(one*fvc::grad(d)));

}

ERR = CALC_ERR(d, mag(fvc::grad(d)), EPS);

Info <<"EPS = "<< EPS << "NITER=" << NITER << tab << "ERR=" << ERR << endl;
//Info << "ERR=" <<ERR<<endl;

Info
<< "  Min(d) = " << min(d).value()
<< "  Max(d) = " << max(d).value()
<< endl;
