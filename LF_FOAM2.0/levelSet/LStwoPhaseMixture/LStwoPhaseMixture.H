/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Class
    Foam::twoPhaseMixture

Description
    A two-phase mixture model

SourceFiles
    LStwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef LStwoPhaseMixture_H
#define LStwoPhaseMixture_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class LStwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class LStwoPhaseMixture
{
protected:

    // Protected data

        word phase1Name_;
        word phase2Name_;

        volScalarField alpha1_;
        volScalarField alpha2_;
        
        volScalarField d_;
        volScalarField H_;

        dictionary LevelSetParameter;
        scalar EPS;
        scalar CFL;

    // Protected Member Functions
        void calculateH();
        void reinitialisation();

public:

    // Constructors

        //- Construct from components
        LStwoPhaseMixture
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~LStwoPhaseMixture()
    {}


    // Member Functions

        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }
        
        //- Return the phase-fraction of phase 1
        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& alpha1()
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 2
        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha2()
        {
            return alpha2_;
        }
        
        //- Return the levelset function d_
        const volScalarField& d() const
        {
            return d_;
        }

        volScalarField& d()
        {
            return d_;
        }

        //-Return the Heaviside function H_
        const volScalarField& H() const
        {
            return H_;
        }
        
        volScalarField& H()
        {
            return H_;
        }

        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
