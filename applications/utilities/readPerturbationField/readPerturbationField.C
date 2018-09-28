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

Application
    readPerturbationField

Description
    For synthetic turbulence data read

\*---------------------------------------------------------------------------*/

#include <fstream>
#include <vector>

#include "fvCFD.H"
#include "vectorList.H"

#include "turbsimBTS.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Utility for testing the syntheticTurbulence library"
    );

#   include "addRegionOption.H"
    argList::addOption
    (
        "fields",
        "wordList",
        "perturbed patches corresponding to timeVaryingMapped BCs"
    );

#   include "setRootCase.H"
    List<word> pertPatches = args.optionRead<wordList>("fields");
    Info<< "Specified patches: " << pertPatches << endl;

#   include "createTime.H"
    IOdictionary pertDict
    (
        IOobject
        (
            "perturbationProperties",
            runTime.time().constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word pertType
    (
        pertDict.lookupOrDefault<word>("syntheticTurbulenceType", "turbsim")
    );

#   include "createNamedMesh.H"

    forAll(pertPatches, patchI)
    {
        label patchID(mesh.boundary().findPatchID(pertPatches[patchI]));
        Info<< "Patch " << pertPatches[patchI] << " id =" << patchID << endl;

        if(pertType == "turbsim")
        {
            // TODO: runtime selectable
            syntheticTurbulence::turbsimBTS ts(mesh.boundary()[patchID]);

            ts.printScaling();
            ts.calcStats();

            // test interpolation between inflow planes
            Field<vector> U;
            U = ts.getPerturbationsAtTime(-1); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(0); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(1); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(1.025); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(1.05); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(599.95); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(599.975); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(600); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(601.025); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;
            U = ts.getPerturbationsAtTime(1201.025); Info<< "U = [ " << U[0] << " " << U[1] << " ... ]" << endl;

        }
        else
        {
            Info<< "Synthetic turbulence type not recognized: " << pertType << endl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
