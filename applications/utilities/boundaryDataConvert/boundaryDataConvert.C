/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    boundaryDataConvert

Description
    Convert boundaryDataPre data (sampled from boundaries) to binary
    boundaryData files for the timeVaryingMapped* class of BCs.

Notes
    boundaryData/<patchName>/points are face centers

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "IFstream.H"
#include "OFstream.H"

#include "vector.H"
#include "pointIOField.H"
#include "AverageIOField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void writeBoundaryDataField
(
    const Time& runTime,
    const fileName prePath,
    const fileName outDir,
    const word fieldName,
    const word timeName,
    const word patchName,
    const List<label>& pointOrder
)
{
    //- read source data
    fileName sourcePath
    (
        prePath
       /timeName
       /patchName
       /pTraits<Type>::typeName+"Field"
       /fieldName
    );
    List<Type> sampledField;
    IFstream(sourcePath)() >> sampledField;

    //- write target data
//
// This writes all the data out correctly but the header incorrectly labels the
// object class as <Type>Field instead of <Type>AverageField; ultimately
// resulting in a runtime error.
//
//    IOField<Type> outputField
//    (
//        IOobject
//        (
//            fieldName,
//            runTime.constant(),
//            "boundaryData"/patchName/timeName,
//            runTime,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE,
//            false
//        )
//    );
//    mkDir(outputField.path());
//    OFstream os
//    (
//        outputField.objectPath(),
//        runTime.writeFormat()
//    );
//    outputField.writeHeader(os);
//    os << "// Average" << endl;
//    os << pTraits<Type>::zero << endl;
//    os << sampledField;
//    //outputField.writeFooter(os);

    AverageIOField<Type> outputField
    (
        IOobject
        (
            fieldName,
            runTime.constant(),
            "boundaryData"/patchName/timeName,
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        sampledField.size()
    );
    mkDir(outputField.path());
    //outputField = sampledField;  // not sure why you can't assign a list to pointField which should be a derived list
    forAll(outputField, faceI)
    {
        //outputField[faceI] = sampledField[faceI];
        outputField[faceI] = sampledField[pointOrder[faceI]];
    }
    outputField.write();

    Info<< "    wrote " << outputField.objectPath() << endl;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "pointsOnly",
        "only output boundaryData points"
    );
    argList::addOption
    (
        "movedPoints",
        "path",
        "directory containing points (in foamFile format) displaced by moveDynamicMesh; these will be the boundaryData points"
    );
    argList::addOption
    (
        "ref",
        "path",
        "directory containing points (in foamFile format) to use as a reference for mapping between old and new points; if specified, boundary points will be reordered"
    );
    argList::addOption
    (
        "patches",
        "fileNameList",
        "list of boundary patches for which to convert to proper boundaryData"
    );
    argList::addBoolOption
    (
        "enforceLapseRate",
        "modify the potential temperature field above a specified height to enforce a specified lapse rate; depends on dictionary constant/lapseDict"
    );

    #include "setRootCase.H"

    // TODO: add boundaryDataPre data location option
    fileName prePath("postProcessing/boundaryDataPre");
    fileName outDir("postProcessing/boundaryData");

    // TODO: option to specify fields to process
    List<word> scalarFields(1, "T");
    List<word> vectorFields(1, "U");

    const bool pointsOnly = args.optionFound("pointsOnly");
    fileName dispPath;
    const bool replacePoints = args.optionReadIfPresent
    (
        "movedPoints",
        dispPath
    );
    fileName refPath;
    const bool haveRef = args.optionReadIfPresent
    (
        "ref",
        refPath
    );
    const bool enforceLapseRate = args.optionFound("enforceLapseRate");
    scalar lapseRate;
    scalar blendStart(0); // suppress compiler warning
    scalar blendEnd;
    scalar blendLayerHeight;

    #include "createTime.H"
//    instantList timeDirs = timeSelector::select0(runTime, args);
//    #include "createNamedMesh.H"
    instantList outputTimes = Time::findTimes(prePath);
    Info<< endl
        << outputTimes.size() << " times sampled: "
        << outputTimes[0].value() << " .. "
        << outputTimes[outputTimes.size()-1].value()
        << endl;

    // read additional parameters for enforcing the lapse rate
    if(enforceLapseRate)
    {
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
        lapseRate = lapseDict.lookupOrDefault<scalar>("lapseDict",1.0);
        lapseRate /= 1000.0;  // convert to K/m
        blendStart = readScalar(lapseDict.lookup("blendStart"));
        blendLayerHeight = lapseDict.lookupOrDefault<scalar>("blendLayerHeight",100.0);
        blendEnd = blendStart + blendLayerHeight;
        Info<< endl
            << "A lapse rate of " << lapseRate
            << " will be enforced blending from " << blendStart
            << " to " << blendEnd
            << endl << endl;
    }

    //
    // Read sampled boundary patches from the first time
    //
    instant time0 = outputTimes[0];
    fileNameList dlist
    (
        readDir(prePath/time0.name(), fileName::DIRECTORY)
    );
    List<word> patches;
    if (args.optionFound("patches"))
    {
        args.optionLookup("patches")() >> patches;
        Info<< "Specified boundary patches: " << patches << endl;
    }
    else
    {
        Info<< "Found patches: " << endl;
        forAll(dlist,patchI)
        {
            Info<< dlist[patchI] << endl;
            patches.append(dlist[patchI]);
        }
    }

    //
    // Write out boundary points at face centers, assuming they're invariant
    //
    List<List<label> > order(patches.size());
    scalar zStart(-1); // for enforceLapseRate
    forAll(patches, patchI)
    {
        word patchName = patches[patchI];

        //- read in face centres
        //  These should be foamFile format, output from boundaryDataPre
        //  sampling output during the simulation; one unfortunate consequence
        //  of the coprocessing is that the ordering of points in the sampling
        //  output is not necessarily standard.
        fileName faceCentersPath
        (
            prePath / time0.name() / patchName / "faceCentres"
        );
        List<vector> faceCenters;
        IFstream(faceCentersPath)() >> faceCenters;

        //- check ranges
        scalar xMin(VGREAT);
        scalar yMin(VGREAT);
        scalar zMin(VGREAT);
        scalar xMax(-VGREAT);
        scalar yMax(-VGREAT);
        scalar zMax(-VGREAT);
        scalar nearest(VGREAT);
        forAll(faceCenters, faceI)
        {
            vector face = faceCenters[faceI];
            if (face.x() < xMin)
            {
                xMin = face.x();
            }
            else if (face.x() > xMax)
            {
                xMax = face.x();
            }
            if (face.y() < yMin)
            {
                yMin = face.y();
            }
            else if (face.y() > yMax)
            {
                yMax = face.y();
            }
            if (face.z() < zMin)
            {
                zMin = face.z();
            }
            else if (face.z() > zMax)
            {
                zMax = face.z();
            }
            if(enforceLapseRate)
            {
                // find layer of cells to calculate an averaged value
                // and start blending from
                scalar dz = face.z() - blendStart;
                if((dz > 0) && (dz < nearest))
                {
                    nearest = dz;
                }
            }
        }
        Info<< patchName << " face centers in"
            << " (" << xMin << ", " << xMax << ")"
            << " (" << yMin << ", " << yMax << ")"
            << " (" << zMin << ", " << zMax << ")"
            << endl;

        //-
        if(enforceLapseRate)
        {
            nearest += blendStart;
            if(zStart < 0)
            {
                zStart = nearest;
                Info<< "averaging layer at z = " << zStart << endl;
            }
            else if(nearest != zStart)
            {
                Info<< "WARNING: mismatch in nearest z value: "
                    << nearest << endl;
            }
        }

        //- find order from reference sampling patch
        //  These should be foamFile format, sampled from a reconstructed mesh
        order[patchI] = List<label>(faceCenters.size());
        if (haveRef)
        {
            fileName refPointsFile
            (
                refPath / patchName / "faceCentres"
            );
            List<vector> refPoints;
            IFstream(refPointsFile)() >> refPoints;
            forAll(refPoints, refI)
            {
                label index(-1);
                for( int i=0; i < faceCenters.size(); i++ )
                {
                    if(Foam::mag(refPoints[refI] - faceCenters[i]) < 1e-5)
                    {
                        index = i;
                        break;
                    }
                }
                if (index < 0)
                {
                    Info<< "Warning: ref point not matched with actual face center"
                        << endl;
                }
                order[patchI][refI] = index;
            }
            Info<< "Note: Reordering assuming the moved mesh and the reference"
                << " mesh have the same standard OpenFOAM ordering." << endl;
        }
        else
        {
            // don't reorder
            forAll(faceCenters, faceI)
            {
                order[patchI][faceI] = faceI;
            }
            Info<< "Note: no reordering will be performed." << endl;
        }

        //- read in moved points, if provided
        //  These should be foamFile format, sampled from a reconstructed mesh
        if (replacePoints)
        {
            fileName pointsFile
            (
                dispPath / patchName / "faceCentres"
            );
            IFstream(pointsFile)() >> faceCenters;
            Info<< "Read moved points from " << pointsFile << endl;
        }

        //- now write out an openfoam IO object
        pointIOField pts
        (
            IOobject
            (
                "points",
                runTime.constant(),
                "boundaryData"/patchName,
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            faceCenters.size()
        );
        //pts = faceCenters;  // not sure why you can't assign a list to pointField which should be a derived list
        forAll(pts, ptI)
        {
            pts[ptI] = faceCenters[ptI];
        }
        mkDir(pts.path());
        pts.write();
        Info<< "Wrote " << pts.objectPath() << endl << endl;

    }

    //
    // Write out requested fields over all times
    //
    if (!pointsOnly)
    {
        forAll(outputTimes, timeI)
        {
            word timeName(outputTimes[timeI].name());
            Info<< "\nt = " << outputTimes[timeI].value() << endl;

            forAll(patches, patchI)
            {
                word patchName(patches[patchI]);
                Info<< "Processing boundary " << patchName << endl;

                //- process scalar fields
                forAll(scalarFields, fieldI)
                {
                    writeBoundaryDataField<scalar>
                    (
                        runTime,
                        prePath,
                        outDir,
                        scalarFields[fieldI],
                        timeName,
                        patchName,
                        order[patchI]
                    );
                }

                //- process vector fields
                forAll(vectorFields, fieldI)
                {
                    writeBoundaryDataField<vector>
                    (
                        runTime,
                        prePath,
                        outDir,
                        vectorFields[fieldI],
                        timeName,
                        patchName,
                        order[patchI]
                    );
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
