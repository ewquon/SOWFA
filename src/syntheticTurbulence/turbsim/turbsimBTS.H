/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::syntheticTurbulence::turbsimBTS

Description
    Handle TurbSim binary full-field time series (.bts) files

\*---------------------------------------------------------------------------*/
#ifndef turbsimBTS_H
#define turbsimBTS_H

#include <fstream>
#include <vector>

#include "vectorList.H"

#include "perturbations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace syntheticTurbulence
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class turbsimBTS
:
    public perturbations
{

private:
    // Private Data

    //- Constants

    // Private Member Functions

    //- Read the entire TurbSim time series
    void read(word fname);
    //void read(word fname, const fvPatch &patch);

public:
    //- Runtime type information
//    TypeName("turbsimBTS");
    
    //- Constructor
//    turbsimBTS
//    (
//        const Time&
//    );
//    turbsimBTS
//    (
//        const Time&,
//        const fvMesh&
//    );
    turbsimBTS
    (
        const fvPatch&
    );
    
    //- Destructor
    virtual ~turbsimBTS()
    {}

    // Public Member Functions

    //- Calculate mean and standard deviation of the perturbation field
    void calcStats();


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam

#endif

// ************************************************************************* //

