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
    const fvPatch& p
)
:
    // Set the pointer to runTime
    runTime_(p.boundaryMesh().mesh().time()),

    // Set the pointer to the mesh
    patch_(p),

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

    correctVariance_
    (
        perturbDict_.lookupOrDefault<bool>
        (
            "correctVariance",
            false
        )
    ),
    perturbedLayerCutoff_
    (
        perturbDict_.lookupOrDefault<word>
        (
            "perturbedLayerCutoff",
            "tanh" // simple | tanh
        )
    ),
    transitionThickness_
    (
        perturbDict_.lookupOrDefault<scalar>
        (
            "transitionThickness",
            50.0
        )
    ),
    transitionEdgeScaling_
    (
        perturbDict_.lookupOrDefault<scalar>
        (
            "transitionEdgeScaling",
            0.97
        )
    ),

    perturbedLayerTop_
    (
        readScalar(perturbDict_.lookup("perturbedLayerTop"))
    ),
    perturbationVariance_
    (
        perturbDict_.lookupOrDefault<scalar>
        (
            "perturbationVariance",
            0.0 // not used if correctVariances is false
        )
    ),

    // Initialize field variables (to be set by derived class)
    Ny(0),
    Nz(0),
    dy(0.0),
    dz(0.0),
    Nt(0),
    dt(0.0),
    period(-1.0),
    periodic(false),
    minU(0), minV(0), minW(0),
    maxU(0), maxV(0), maxW(0),
    stdU(0), stdV(0), stdW(0),

    mapperPtr_(NULL)

{
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

inline scalar perturbations::tanhScaling(scalar z) const
{
    return 0.5 * Foam::tanh( 2 * tanhParam * (z-perturbedLayerTop_) / transitionThickness_ ) + 0.5;
}

void perturbations::setScaling()
{
    scaling.resize(Ny*Nz);
    scalar scalingFactor = 1.0;
    if(correctVariance_)
    {
        // TODO: calculate scalingFactor if correctVariances is true
    }
    if(perturbedLayerCutoff_ == "simple")
    {
        Info<< "Setup scaling with simple cutoff" << endl;
        forAll(scaling, I)
        {
            if(points[I].z() < perturbedLayerTop_)
            {
                scaling[I] = scalingFactor; //vector::one;
            }     
            else  
            {     
                scaling[I] = 0.0; //vector::zero;
            }
        }
    }
    else if(perturbedLayerCutoff_ == "tanh")
    {
        Info<< "Setup scaling with tanh cutoff" << endl;
        tanhParam = -Foam::atanh(2*transitionEdgeScaling_ - 1);
        forAll(scaling, I)
        {
            scaling[I] = scalingFactor * tanhScaling(points[I].z());
        }
    }
    else
    {
        FatalError
            << "Cutoff method " << perturbedLayerCutoff_
            << " not recognized" << nl
            << exit(FatalError);
    }
}

void perturbations::setupMapper()
{
    // set up inflow points field (i.e., the "samplePoints")
    pointField inflowPoints(points);

    Info<< inflowPoints << endl;
    Info<< patch_.Cf() << endl;

    // allocate the planar interpolator
    mapperPtr_.reset
    (
        new pointToPointPlanarInterpolation
        (
            inflowPoints, // source
            patch_.Cf(), // destination
            1e-5, // default perturb_ value from TVM* BCs
            false // nearestOnly flag
        )
    );
}

Field<vector> perturbations::getPerturbationsAtTime
(
    scalar t,
    scalar ang
)
{
//    Pout<< "Requested perturbations at time " << t
//        << " with rotation angle "
//        << 180.0/Foam::constant::mathematical::pi * ang << " deg"
//        << endl;

    if(mapperPtr_.empty()) setupMapper();

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
            Info<< "t=" << t;
            t -= int(t/period)*period;
            Info<< " mapped to t=" << t;
        }
        label i;
        for(i=1; i<times.size(); ++i)
        {
            if(times[i] > t) break;
        }
        if(i<times.size())
        {
            Info<< " between " << times[i-1] << " and " << times[i] << endl;
            U = perturb[i-1]
                + (perturb[i]-perturb[i-1])/(times[i]-times[i-1]) * (t-times[i-1]);
        }
        else
        {
            Info<< " between " << times[i-1] << " and " << period << endl;
            U = perturb[i-1]
                + (perturb[0]-perturb[i-1])/(period-times[i-1]) * (t-times[i-1]);
        }
    }

    // enforce cutoff
    forAll(U, faceI)
    {
        U[faceI] *= scaling[faceI];
    }

    // rotate to align streamwise components
    Field<vector> rotatedU(U.size(),vector::zero);
    forAll(U, faceI)
    {
        rotatedU[faceI].x() = U[faceI].x()*Foam::cos(ang) - U[faceI].y()*Foam::sin(ang);
        rotatedU[faceI].y() = U[faceI].x()*Foam::sin(ang) + U[faceI].y()*Foam::cos(ang);
        rotatedU[faceI].z() = U[faceI].z();
    }

    // return field mapped to boundary patch 
    if(rotatedU.size() != mapperPtr_().sourceSize())
    {
        FatalErrorIn
        (
            "perturbations<Type>::"
            "getPerturbationsAtTime()"
        )   << "Number of values (" << rotatedU.size()
            << ") differs from the number of points ("
            <<  mapperPtr_().sourceSize()
            << ")"
            << exit(FatalError);
    }

    return mapperPtr_().interpolate(rotatedU);
}

void perturbations::printScaling()
{
    scalar zval;
    Info<< "Scaling function: [";
    zval = points[0].z();
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    Info<< " ...";
    zval = perturbedLayerTop_ - 2.0*transitionThickness_;
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    zval = perturbedLayerTop_ - 1.0*transitionThickness_;
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    zval = perturbedLayerTop_ - 0.5*transitionThickness_;
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    zval = perturbedLayerTop_;
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    zval = perturbedLayerTop_ + 0.5*transitionThickness_;
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    zval = perturbedLayerTop_ + 1.0*transitionThickness_;
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    zval = perturbedLayerTop_ + 2.0*transitionThickness_;
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    Info<< " ...";
    zval = points[points.size()-1].z();
    Info<< " (" << zval << ", " << tanhScaling(zval) << ")";
    Info<< "]" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam
