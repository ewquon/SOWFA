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
    const word patchName    
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
        outputField[faceI] = sampledField[faceI];
    }
    outputField.write();

    Info<< "        wrote " << outputField.objectPath() << endl;
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
        "cellDisplacement",
        "path",
        "directory containing sampled boundary vectorField 'cellDisplacement' in foamFile format"
    );
    argList::addOption
    (
        "patches",
        "fileNameList",
        "list of boundary patches for which to convert to proper boundaryData"
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
    const bool addDisplacement = args.optionReadIfPresent
    (
        "cellDisplacement",
        dispPath
    );

    #include "createTime.H"
//    instantList timeDirs = timeSelector::select0(runTime, args);
//    #include "createNamedMesh.H"
    instantList outputTimes = Time::findTimes(prePath);
    Info<< outputTimes.size() << " times sampled: "
        << outputTimes[0].value() << " .. "
        << outputTimes[outputTimes.size()-1].value()
        << endl;

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
    Info<< "boundary patches: " << patches << endl;

    //
    // Write out boundary points at face centers, assuming they're invariant
    // TODO: translate plane here, if needed
    //
    forAll(patches, patchI)
    {
        word patchName = patches[patchI];

        //- read in face centres
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
        }
        Info << patchName << " face centers in"
            << " (" << xMin << ", " << xMax << ")"
            << " (" << yMin << ", " << yMax << ")"
            << " (" << zMin << ", " << zMax << ")"
            << endl;

        //- read in cell displacements if provided
        List<vector> displacement(faceCenters.size(), vector::zero);
        if (addDisplacement)
        {
            fileName cellDisplacementPath
            (
                dispPath / patchName / "vectorField" / "cellDisplacement"
            );
            IFstream(cellDisplacementPath)() >> displacement;
            Info<< "Read " << cellDisplacementPath << endl;
            //Info<< displacement[0] << " " << displacement[displacement.size()-1] << endl;
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
            pts[ptI] = faceCenters[ptI] + displacement[ptI];
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
        forAll(patches, patchI)
        {
            word patchName(patches[patchI]);
            Info<< "\nProcessing boundary " << patchName << endl;

            forAll(outputTimes, timeI)
            {
                word timeName(outputTimes[timeI].name());
                Info<< "  t = " << outputTimes[timeI].value() << endl;

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
                        patchName
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
                        patchName
                    );
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
