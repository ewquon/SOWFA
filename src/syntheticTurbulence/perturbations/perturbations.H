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
    Foam::syntheticTurbulence::perturbations

Description
    Class for reading/processing/generating synthetic turbulence fields for
    mitigating the inflow turbulence generation problem through spectral
    enrichment. Perturbations from these synthetic fields are superimposed onto
    the timeVaryingMapped* class of boundary conditions. This base class will
    provide routines to facilitate interpolating in time and mapping the 
    perturbations to the specified boundary patch(es).

\*---------------------------------------------------------------------------*/
#ifndef perturbations_H
#define perturbations_H

#include "fvCFD.H"
#include "pointToPointPlanarInterpolation.H"

//#include "syntheticTurbulence.H" // is this needed?

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace syntheticTurbulence
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: make template<class Type>
class perturbations
{

protected:
    //- Constants

        //- Runtime pointer.
        const Time& runTime_;

        //- Patch pointer.
        const fvPatch& patch_;

        //- Input dictionary
        IOdictionary perturbDict_;

        //- Input properties
        word perturbType_;
        //List<word> perturbPatches_;

        //- Scaling properties
        bool correctVariance_;
        word perturbedLayerCutoff_;
        scalar transitionThickness_;
        scalar transitionEdgeScaling_;

        //- Variable scaling properties
        scalar perturbedLayerHeight_;
        scalar perturbationVariance_;

        List<List<scalar> > perturbationControlTable_;

        //- TODO: Automatic update of scaling parameters
        word perturbationHeightControl_;
    //- Protected field parameters

        //- Inflow plane dimensions
        //  (in plane-aligned coordinates: y=lateral, z=vertical)
        int Ny, Nz;

        //- Inflow plane spacings
        float dy, dz;

        //- Inflow plane times
        int Nt;

        //- Inflow time-step size
        //float dt;
        scalar dt;

        //- Time-series period (< 0 if aperiodic)
        scalar period;
        bool periodic;

        //- Turbulence statistics
        scalar minU, minV, minW;
        scalar maxU, maxV, maxW;
        scalar stdU, stdV, stdW;

    //- Protected field variables

        //- Time vector
        DynamicList<scalar> times;

        //- Coordinate vectors
        DynamicList<vector> points;

        //- List of inflow perturbation fields (one list item per time instance)
        // TODO: make template<class Type>
        DynamicList<List<vector> > perturb;

        //- List of perturbation scaling factors
        DynamicList<scalar> scaling;

    //- Protected functions

        //- Reorient inflow plane to correctly map between synthetic inflow and boundary patches
        //  using pointToPointPlanarInterpolation. This must be called BEFORE the mapper is set up.
        void reorientInflowPlane(vector&);

        //- Set scaling factors
        void setScaling();

private:
    //- Private variables

        //- Probe location (for diagnostics)
        bool haveProbe;
        label probeID;

        //- Scaling parameters
        scalar tanhParam; // "k", where s(z) = 0.5*tanh(k*(z-ztop)/(width/2)) + 0.5

        //- Planar interpolation
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- Saved field (so we don't repeat the same interpolation multiple times)
        scalar tlast;
        Field<vector> Ulast;

        //- Perturbation controls
        List<scalar> controlTime;
        List<scalar> controlHeight;
    
    //- Private functions

        //- Reformat perturbationControlTable for ease of use
        void setupControlTable();

        //- Diagnostics
        void setupProbe();
    
        //- Planar interpolation
        void setupMapper();

        //- Update perturbedLayerHeight from perturbationControlTable
        void updatePerturbationHeight(scalar);

        //- Scaling functions
        inline scalar tanhScaling(scalar) const;

public:
    //- Runtime type information
//    TypeName("perturbations");

// TypeName:
//    ClassName("perturbations");
//    ClassNameNoDebug("perturbations");

// ClassNameNoDebug:
//    static const char* typeName_() { return "perturbations"; }
//    static const ::Foam::word typeName;
//    static int debug;
//    virtual const word& type() const { return typeName; }
    
    //- Constructor
    // initialize with top-level database only (no patch-specific information)
//    perturbations
//    (
//        const Time&
//    );
    // initialize with mesh reference (no patch-specific information)
//    perturbations
//    (
//        const fvMesh&
//    );
    // initialize with patch reference
    perturbations
    (
        const fvPatch&
    );
    
    //- Destructor
    virtual ~perturbations()
    {}

    // Access

        word perturbType() { return perturbType_; }
        scalar perturbedLayerHeight() { return perturbedLayerHeight_; }

    // Public Member Functions

        //- print out the scaling function (for debug)
        void printScaling();

        //- get statistics of the perturbation field
        virtual void calcStats() = 0;

        //- get perturbation field at time, with specified rotation
        // TODO: make template<class Type>
        const Field<vector>& getPerturbationsAtTime
        (
            scalar t,
            scalar ang = 0
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam

#endif

// ************************************************************************* //
