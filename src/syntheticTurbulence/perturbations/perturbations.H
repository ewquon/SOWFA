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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace syntheticTurbulence
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class perturbations
{

protected:
    //- Constants

        //- Runtime pointer.
        const Time& runTime_;

        //- Mesh pointer.
        const fvMesh& mesh_;

        //- Input dictionary
        IOdictionary perturbDict_;

        //- Input properties
        word perturbType_;
        List<word> perturbPatches_;

    //- Field parameters

        //- Inflow plane dimensions
        //  (in plane-aligned coordinates: y=lateral, z=vertical)
        int Ny, Nz;

        //- Inflow plane spacings
        float dy, dz;

        //- Inflow plane times
        int Nt;

        //- Inflow time-step size
        float dt;

        //- Time-series period (< 0 if aperiodic)
        scalar period;
        bool periodic;

    //- Field variables

        //- Time vector
        DynamicList<scalar> t;

        //- Coordinate vectors
        DynamicList<List<vector> > points;

        //- List of inflow perturbation fields (one list item per time instance)
        DynamicList<List<vector> > perturb;


    // Private Member Functions

public:
    
    //- Constructor
    perturbations
    (
        const volVectorField& U
    );
    
    //- Destructor
    virtual ~perturbations()
    {}

    // Public Member Functions

        //- get statistics of the inflow field
        virtual void calcStats();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam

#endif

// ************************************************************************* //
