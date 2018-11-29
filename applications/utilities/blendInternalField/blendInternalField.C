/*---------------------------------------------------------------------------*\
This file was modified or created at the National Renewable Energy
Laboratory (NREL) on January 6, 2012 in creating the SOWFA (Simulator for
Offshore Wind Farm Applications) package of wind plant modeling tools that
are based on the OpenFOAM software. Access to and use of SOWFA imposes
obligations on the user, as set forth in the NWTC Design Codes DATA USE
DISCLAIMER AGREEMENT that can be found at
<http://wind.nrel.gov/designcodes/disclaimer.html>.
\*---------------------------------------------------------------------------*/

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

Application
    blendInternalField

Description
    Updates the internal field based on specified (temperature, for now) field
    information. The specified field F(z) is blended with the existing field 
    F(x,y,z). This is experimental and design to be used in tandem with
    `boundaryDataConvert -enforceLapseRate`; it depends on constant/lapseDict
    that is used by that code.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interpolateXY.H"
#include "interpolateSplineXY.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Read in the existing solution files.   

Info << "Reading field T" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


// Update the interior fields.

Info << "Updating internal T field..." << endl;
forAll(T,cellI)
{
    //T[cellI] = Tbottom;
}


// Update the boundary field.
//Info << "Updating boundaries..." << endl;
//T.correctBoundaryConditions();


// Write out the updated fields.
Info<< "Writing field T" << endl;
T.write(); 


return 0;
}


// ************************************************************************* //

