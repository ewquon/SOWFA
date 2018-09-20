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
    const Time& runTime
)
:
    // Set the pointer to runTime
    runTime_(runTime),

    // Set the pointer to the mesh
//    mesh_(mesh),

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
    perturbedLayerTop_
    (
        readScalar(perturbDict_.lookup("perturbedLayerTop"))
    ),

    // Initialize field variables (to be read in by derived class)
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

void perturbations::setScaling()
{
    scaling.resize(Ny*Nz);
    Info<< "Setup simple scaling";
    forAll(scaling, I)
    {
        if(points[I].z() < perturbedLayerTop_)
        {
            scaling[I] = vector::one;
        }     
        else  
        {     
            scaling[I] = vector::zero;
        }
    }
}

List<vector> perturbations::getPerturbationsAtTime(scalar t) const
{
    Info<< "Requested perturbations at time " << t << endl;
    List<vector> U;
    if( t < times[0] )
    {
        Info<< "WARNING: t < " << times[0] << endl;
        U = perturb[0];
    }
    else if( !periodic && t > times[times.size()-1] )
    {
        Info<< "WARNING: t > " << times[times.size()-1] << endl;
        U = perturb[times.size()-1];
    }
    else
    {
        if(periodic)
        {
            t -= int(t/period)*period;
            Info<< "  mapped to t = " << t;
        }
        label i;
        for(i=1; i<times.size(); ++i)
        {
            if(times[i] > t) break;
        }
        if(i<times.size())
        {
            Info<< "  between " << times[i-1] << " and " << times[i] << endl;
            U = perturb[i-1]
                + (perturb[i]-perturb[i-1])/(times[i]-times[i-1]) * (t-times[i-1]);
        }
        else
        {
            Info<< "  between " << times[i-1] << " and " << period << endl;
            U = perturb[i-1]
                + (perturb[0]-perturb[i-1])/(period-times[i-1]) * (t-times[i-1]);
        }
    }
    return U;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam
