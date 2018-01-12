/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LStwoPhaseMixture.H"
#include "FUNCTIONS.H"
#include "fvc.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LStwoPhaseMixture::calculateH()
{
    H_ == CALC_HEAVISIDE(d_,EPS);
}


void Foam::LStwoPhaseMixture::reinitialisation()
{
    dimensionedScalar one = dimensionedScalar("one", dimensionSet(0,1,0,0,0,0,0), 1.0);
    const surfaceScalarField deltaCoeff = d_.mesh().deltaCoeffs();
    scalar DX = (1.0/max(deltaCoeff).value());
    scalar Dtau  = CFL*DX;
    int NITER = EPS/Dtau;
    scalar ERR (1.0);

    d_  == 1.5*DX*(alpha1_-0.5);

    volScalarField d0
    (
        IOobject
        (
            "d0",
            d_.instance(),
            d_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        d_.mesh(),
        dimensionedScalar("d0",dimensionSet(0,0,0,0,0,0,0),scalar(0)),
        d_.boundaryField().types()
    );

    d0 == d_;

    for (int loop = 0; loop <= NITER; loop++)
    {

        // calculate distance function field
        d_ += Dtau*sign(d0)*(scalar(1)-mag(one*fvc::grad(d_)));

    }

    ERR = CALC_ERR(d_, mag(fvc::grad(d_)), EPS);
    Info <<"EPS = "<< EPS << "NITER=" << NITER << tab << "ERR=" << ERR << endl;
    Info << "  Min(d) = " << min(d_).value() << "  Max(d) = " << max(d_).value()<< endl;

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LStwoPhaseMixture::LStwoPhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phase1Name_(wordList(dict.lookup("phases"))[0]),
    phase2Name_(wordList(dict.lookup("phases"))[1]),

    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - alpha1_
    ),
    
    d_
    (
        IOobject
        (
            "d",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("d",dimensionSet(0,0,0,0,0,0,0),scalar(0)),
        alpha1_.boundaryField().types()
    ),
    
    H_
    (
         IOobject
        (
            "H",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("H",dimensionSet(0,0,0,0,0,0,0),scalar(0))
    ),

    LevelSetParameter(dict.subDict("LevelSetParameter")),

    EPS(readScalar(LevelSetParameter.lookup("interface_thickness"))),

    CFL(readScalar(LevelSetParameter.lookup("cfl_number")))
{
    reinitialisation();
    calculateH();
}

void Foam::LStwoPhaseMixture::correct()
{
    reinitialisation();
    calculateH();
}

// ************************************************************************* //
