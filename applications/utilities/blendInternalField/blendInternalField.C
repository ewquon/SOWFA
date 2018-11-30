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

// Read dictionary used by boundaryDataConvert
IOdictionary lapseDict
(
    IOobject
    (
        "lapseDict",
        runTime.time().constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
scalar lapseRate = lapseDict.lookupOrDefault<scalar>("lapseDict",1.0);
lapseRate /= 1000.0;  // convert to K/m
scalar blendStart = readScalar(lapseDict.lookup("blendStart"));
scalar blendLayerHeight = lapseDict.lookupOrDefault<scalar>("blendLayerHeight",100.0);
scalar blendEnd = blendStart + blendLayerHeight;
Info<< "\nA lapse rate of " << lapseRate << " K/m"
    << " will be enforced blending from " << blendStart
    << " to " << blendEnd
    << endl << endl;

// for selecting a layer of cells to determine Tbottom
scalar avgtol = lapseDict.lookupOrDefault<scalar>("averagingTolerance",1.0);

// Find bottom of the inversion layer (assuming grid is structured Cartesian,
// i.e., if we find the nearest level to the inversion, we'll have the entire
// layer of cells)
scalar nearest(VGREAT);
forAll(T,cellI)
{
    scalar dz = mesh.C()[cellI].z() - blendStart;
    if((dz > 0) && (dz < nearest))
    {
        nearest = dz;
    }
}
scalar zStart = nearest + blendStart;

// Find Tbottom as the layer mean
scalar Tsum(0);
label layercount(0);
forAll(T,cellI)
{
    scalar dz = mesh.C()[cellI].z() - zStart;
    if((dz > 0) && (abs(dz) < avgtol))
    {
        Tsum += T[cellI];
        layercount++;
    }
}
dimensionedScalar Tbottom("Tbottom", dimTemperature, Tsum/layercount);
Info<< "Calculated " << Tbottom
    << " at z = " << zStart
    << " averaged over " << layercount << " cells"
    << endl;

volScalarField blending
(
    IOobject
    (
        "blending",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "blending",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        0.0
    )
);

volScalarField Tlapse
(
    IOobject
    (
        "Tlapse",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Tlapse",dimTemperature,0.0)
);

forAll(T, cellI)
{
    scalar z = mesh.C()[cellI].z();
    // specified lapse rate region
    if(z > blendEnd)
    {
        //Tlapse[cellI] = z - blendStart;
        Tlapse[cellI] = z - zStart;
        blending[cellI] = 1.0;
    }
    // blended region
    //else if(z >= blendStart)
    else if(z >= zStart)
    {
        //Tlapse[cellI] = z - blendStart;
        Tlapse[cellI] = z - zStart;
        //blending[cellI] = (z - blendStart) / (blendEnd - blendStart);
        //blending[cellI] = Tlapse[cellI] / blendLayerHeight;
        blending[cellI] = Tlapse[cellI] / (blendEnd - zStart);
    }
    // DEBUG:
    //Info<< z << " " << blending[cellI] << " " << Tlapse[cellI] << endl;
}
Tlapse = lapseRate*Tlapse + Tbottom;


// Update the interior fields.

Info << "Updating internal T field..." << endl;
T = blending*Tlapse + (1-blending)*T;


// Update the boundary field.
//Info << "Updating boundaries..." << endl;
//T.correctBoundaryConditions();


// Write out the updated fields.
Info<< "Writing field T" << endl;
T.write(); 
Tlapse.write();
blending.write();


Info<< "\nEnd." << endl;
return 0;
}


// ************************************************************************* //

