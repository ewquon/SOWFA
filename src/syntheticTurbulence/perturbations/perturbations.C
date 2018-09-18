/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "perturbations.H"

namespace Foam
{

namespace syntheticTurbulence
{

// * * * * * * * * * * * * * *  Constructor  * * * * * * * * * * * * * * * * //

perturbations::perturbations
(
    const volVectorField& U
)
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    // Read the dictionary describing the perturbation strategy
    perturbDict_
    (
        IOobject
        (
            "perturbationProperties",
            runTime_.time().constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    perturbType_
    (
        perturbDict_.lookupOrDefault<word>
        (
            "syntheticTurbulenceType",
            "turbsim"
        )
    ),
    perturbPatches_
    (
        perturbDict_.lookup("perturbedPatches")
    ),

    // Initialize field variables
    Ny(0),
    Nz(0),
    dy(0.0),
    dz(0.0),
    Nt(0),
    dt(0.0),
    period(-1.0),
    periodic(false)

{}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam
