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

    perturbedLayerHeight_
    (
        readScalar(perturbDict_.lookup("perturbedLayerHeight"))
    ),
    perturbationVariance_
    (
        perturbDict_.lookupOrDefault<scalar>
        (
            "perturbationVariance",
            1.0 // not used if correctVariances is false
        )
    ),

    perturbationControlTable_
    (
        perturbDict_.lookup("perturbationControlTable")
    ),

    // Perturbation control
    perturbationHeightControl_
    (
        perturbDict_.lookupOrDefault<word>
        (
            "perturbationHeightControl",
            "constant"
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

    // Diagnostics
    haveProbe(false),
    probeID(-1),

    // Tanh scaling parameter
    tanhParam(-Foam::atanh(2*transitionEdgeScaling_ - 1)),

    // Interpolation vars
    mapperPtr_(NULL),
    tlast(-1)

{
    setupControlTable();
    setupProbe();
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void perturbations::setupControlTable()
{
    label N = perturbationControlTable_.size();
    controlTime = List<scalar>(N, 0.0);
    controlHeight = List<scalar>(N, 0.0);
    forAll(perturbationControlTable_, timeI)
    {
        controlTime[timeI] = perturbationControlTable_[timeI][0];
        controlHeight[timeI] = perturbationControlTable_[timeI][1];
    }
}

void perturbations::setupProbe()
{
    // find center of patch (to find probe location), assuming patch is planar
    const boundBox& globalBb = patch_.boundaryMesh().mesh().bounds();
    point domainCenter = 0.5*(globalBb.min() + globalBb.max());
    domainCenter.z() = 100.0; // hard-coded probe height
    point patchCenter(vector::zero);
    if(patch_.size() > 0)
    {
        vector norm = patch_.Sf()[0] / patch_.magSf()[0];
        vector patchDist = ((domainCenter - patch_.Cf()[0]) & norm) * norm;
        patchCenter = domainCenter - patchDist;

        //Pout<< "   patch center: " << patchCenter
        //    << " (" << patch_.name() << ")"
        //    << endl;
    }

    // find probe location
    scalar minDist2(Foam::VGREAT);
    List<scalar> r2(patch_.size());
    if(patch_.size() > 0)
    {
        vector r;
        forAll(patch_.Cf(), faceI)
        {
            r = patch_.Cf()[faceI] - patchCenter;
            r2[faceI] = r & r;
            if(r2[faceI] < minDist2)
            {
                minDist2 = r2[faceI];
            }
        }
        //Pout<< "   minDist2: " << minDist2
        //    << " (" << patch_.name() << ")"
        //    << endl;
    }

    // find proc that contains the center point
    reduce(minDist2, minOp<scalar>());
    //Info<< "  minDist2 = " << minDist2 << endl;
    forAll(patch_.Cf(), faceI)
    {
        if(r2[faceI] == minDist2)
        {
            haveProbe = true;
            probeID = faceI;
            Pout<< "Probe on " << patch_.name()
                << " patch at " << patch_.Cf()[probeID] << endl;
            break;
        }
    }
}


inline scalar perturbations::tanhScaling(scalar z) const
{
    return 0.5*Foam::tanh( 2 * tanhParam * (z-perturbedLayerHeight_) / transitionThickness_ ) + 0.5;
}

inline scalar perturbations::linearScaling(scalar z) const
{
    if(z > perturbedLayerHeight_)
    {
        return 0;
    }
    else
    {
        return (perturbedLayerHeight_ - z) / perturbedLayerHeight_;
    }
}

void perturbations::setScaling()
{
    scaling.resize(Ny*Nz);
    scalar scalingFactor = perturbationVariance_;
    if(correctVariance_)
    {
        // TODO: calculate scalingFactor if correctVariances is true
    }
    if(perturbedLayerCutoff_ == "simple")
    {
        Info<< "Updated scaling with simple cutoff" << endl;
        forAll(scaling, I)
        {
            if(points[I].z() < perturbedLayerHeight_)
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
        Info<< "Updated scaling with tanh cutoff" << endl;
        forAll(scaling, I)
        {
            scaling[I] = scalingFactor * tanhScaling(points[I].z());
        }
    }
    else if(perturbedLayerCutoff_ == "linear")
    {
        Info<< "Updated linear scaling" << endl;
        forAll(scaling, I)
        {
            scaling[I] = scalingFactor * linearScaling(points[I].z());
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

void perturbations::reorientInflowPlane
(
    vector& norm0
)
{
    if(patch_.size() > 0)
    {
        vector norm = patch_.Sf()[0] / patch_.magSf()[0];
        if(norm != norm0)
        {
            vector cross = norm ^ norm0;
            scalar ang = Foam::asin(cross.z());
            Pout<< "Reorienting inflow plane from " << norm0 << " to " << norm
                << " for patch " << patch_.name()
                << " (rotate by "
                << 180.0/Foam::constant::mathematical::pi * ang << " deg)"
                << endl;

            // Rotate...
            List<vector> rotatedPoints(points.size(), vector::zero);
            forAll(points, ptI)
            {
                rotatedPoints[ptI].x() = points[ptI].x()*Foam::cos(ang) - points[ptI].y()*Foam::sin(ang);
                rotatedPoints[ptI].y() = points[ptI].x()*Foam::sin(ang) + points[ptI].y()*Foam::cos(ang);
                rotatedPoints[ptI].z() = points[ptI].z();
            }
            // ...and mirror (poor-man's periodicity) -- acknowledged that fields will be in
            // inconsistent in terms of coherence at the domain corners... but turbsim isn't
            // perfect anyway (not divergence free)
            // TODO: need to come up with a better (consistent) approach for Gabor KS
            points = cmptMag(rotatedPoints); // i.e., absolute value
        }
    }
}

void perturbations::setupMapper()
{
    if(patch_.size() > 0)
    {
        Pout<< "Setting up patch mapper for "
            << patch_.size() << " points on patch "
            << patch_.name()
            << endl;
    }

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

void perturbations::updatePerturbationHeight
(
    scalar t
)
{
    if( t < controlTime[0] )
    {
        Info<< "WARNING: t < " << controlTime[0] << endl;
        perturbedLayerHeight_ = controlHeight[0];
    }
    else if( t > controlTime[controlTime.size()-1] )
    {
        Info<< "WARNING: t > " << controlTime[controlTime.size()-1] << endl;
        perturbedLayerHeight_ = controlHeight[controlTime.size()-1];
    }
    else
    {
        Info<< "Updating perturbation height from control table for t= " << t;
        label i;
        for(i=1; i<controlTime.size(); ++i)
        {
            if(controlTime[i] > t) break;
        }
        if(i<controlTime.size())
        {
            Info<< " between " << controlTime[i-1]
                << " and " << controlTime[i] << endl;
            perturbedLayerHeight_ = controlHeight[i-1]
                + (controlHeight[i]-controlHeight[i-1])/(controlTime[i]-controlTime[i-1])
                    * (t-controlTime[i-1]);
        }
        else
        {
            Info<< " (last time)" << endl;
            perturbedLayerHeight_ = controlHeight[controlTime.size()-1];
        }
    }
    Info<< "  height = " << perturbedLayerHeight_ << " m" << endl;
}

const Field<vector>& perturbations::getPerturbationsAtTime
(
    scalar t,
    scalar ang
)
{
//    Pout<< "Requested perturbations at time " << t
//        << " with rotation angle "
//        << 180.0/Foam::constant::mathematical::pi * ang << " deg"
//        << endl;

    if(t==tlast) return Ulast;

    // first time only
    if(mapperPtr_.empty()) setupMapper();

    // update scaling if height control is "table" or "auto"
    if(perturbationHeightControl_ != "constant")
    {
        updatePerturbationHeight(t);
        setScaling();
    }

    // now interpolate
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
        Info<< "Retrieving inflow for t= " << t;
        scalar t0 = t;
        if(periodic)
        {
            t0 -= int(t0/period)*period;
            Info<< " mapped to t= " << t0;
        }

        tlast = t;
        label i;
        for(i=1; i<times.size(); ++i)
        {
            if(times[i] > t0) break;
        }
        if(i<times.size())
        {
            Info<< " between " << times[i-1] << " and " << times[i] << endl;
            U = perturb[i-1]
                + (perturb[i]-perturb[i-1])/(times[i]-times[i-1]) * (t0-times[i-1]);
        }
        else
        {
            Info<< " between " << times[i-1] << " and " << period << endl;
            U = perturb[i-1]
                + (perturb[0]-perturb[i-1])/(period-times[i-1]) * (t0-times[i-1]);
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

    Ulast = mapperPtr_().interpolate(rotatedU);

    if(haveProbe)
    {
        Pout<< patch_.name() << " boundary probe : " << Ulast[probeID] << endl;
    }

    return Ulast;
}

void perturbations::printScaling()
{
    if(patch_.size() > 0)
    {
        scalar zval;
        Pout<< "tanh scaling function: [";
        zval = points[0].z();
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        Pout<< " ...";
        zval = perturbedLayerHeight_ - 2.0*transitionThickness_;
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        zval = perturbedLayerHeight_ - 1.0*transitionThickness_;
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        zval = perturbedLayerHeight_ - 0.5*transitionThickness_;
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        zval = perturbedLayerHeight_;
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        zval = perturbedLayerHeight_ + 0.5*transitionThickness_;
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        zval = perturbedLayerHeight_ + 1.0*transitionThickness_;
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        zval = perturbedLayerHeight_ + 2.0*transitionThickness_;
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        Pout<< " ...";
        zval = points[points.size()-1].z();
        Pout<< " (" << zval << ", " << tanhScaling(zval) << ")";
        Pout<< "]" << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam
