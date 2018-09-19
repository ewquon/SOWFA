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

/*
// Read a TurbSim binary full-field time series (.bts) file and map the
// fluctuating velocity field to a specified boundary patch
void read_bts(word fname)
//void read_bts(word fname, const fvPatch &patch)
{
    std::ifstream bts(fname.c_str(), std::ios::binary);
    if(bts)
    {
        bts.seekg(0, std::ios_base::end);
        Info<< "Reading " << fname
            << " (" << bts.tellg() << " bytes)"
            << endl;
        bts.seekg(0, std::ios_base::beg);

        short id;

        int NY, NZ, Ntower;
        int N; // number of times
        float dt,dy,dz;

        float uhub;
        float zhub, zbot;

        float Vslope[3];
        float Vintercept[3];

        bts.read(reinterpret_cast<char*>(&id), sizeof(id));
        //assert( id==7 || id==8 )
        if(id==7) Info << "non-periodic" << endl;
        else if (id==8) Info << "periodic" << endl;
        else Info << "Warning: Unknown BTS file ID: " << id << endl;

        // read resolution settings
        bts.read(reinterpret_cast<char*>(&NZ), sizeof(NZ));
        bts.read(reinterpret_cast<char*>(&NY), sizeof(NY));
        Info<< "NY, NZ = " << NY << " " << NZ << endl;

        bts.read(reinterpret_cast<char*>(&Ntower), sizeof(Ntower));
        Info<< "Ntower = " << Ntower << endl;

        bts.read(reinterpret_cast<char*>(&N), sizeof(N));
        Info<< "N = " << N << endl;

        bts.read(reinterpret_cast<char*>(&dy), sizeof(dy));
        bts.read(reinterpret_cast<char*>(&dz), sizeof(dz));
        bts.read(reinterpret_cast<char*>(&dt), sizeof(dt));
        Info<< "dy, dz = " << dy << " " << dz << endl;
        Info<< "dt = " << dt << endl;
        Info<< "period = " << N * dt << endl;

        // read reference values
        bts.read(reinterpret_cast<char*>(&uhub), sizeof(uhub));
        bts.read(reinterpret_cast<char*>(&zhub), sizeof(zhub));
        bts.read(reinterpret_cast<char*>(&zbot), sizeof(zbot));
        Info<< "uhub = " << uhub << endl;
        Info<< "zhub = " << zhub << endl;
        Info<< "zbot (NOT USED) = " << zbot << endl;

        // read scaling factors
        Info<< "Vslope, intercept : " << endl;
        for( int i=0; i < 3; ++i )
        {
            bts.read(reinterpret_cast<char*>(&Vslope[i]), sizeof(float));
            bts.read(reinterpret_cast<char*>(&Vintercept[i]), sizeof(float));
            Info<< "  " << Vslope[i] << " " << Vintercept[i] << endl;
        }

        // read turbsim info string
        int nchar;
        bts.read(reinterpret_cast<char*>(&nchar), sizeof(nchar));
        char infostr[nchar];
        bts.read(infostr, nchar);
        Info<< infostr << endl;

        // read normalized data and dimensionalize it
        std::vector<short> ivals(3*NY*NZ*N);
//        std::vector<short> ivals_tow(3*Ntower*N);
//        std::vector<float> U(3*NY*NZ*N);
//        std::vector<float> Utower(3*Ntower*N);
        std::vector<float> U(NY*NZ*N);
        std::vector<float> V(NY*NZ*N);
        std::vector<float> W(NY*NZ*N);

        Info<< "Reading normalized velocity data...";
        Foam::flush(Info);
        bts.read(reinterpret_cast<char*>(ivals.data()), ivals.size()*sizeof(short));
        Info<< " done!" << endl;

        int n=0;

        Info<< "Dimensionalizing velocity data...";
        int m = 0;
        for(int itime=0; itime<N; ++itime)
        {
            for(int k=0; k<NZ; ++k)
            {
                for(int j=0; j<NY; ++j)
                {
                    U[n] = (ivals[m] - Vintercept[0]) / Vslope[0] - uhub;
                    V[n] = (ivals[m+1] - Vintercept[1]) / Vslope[1];
                    W[n] = (ivals[m+2] - Vintercept[2]) / Vslope[2];
                    m += 3;
                    ++n;
                }
            }
        }
        Info<< " done!" << endl;

//        if(Ntower > 0)
//        {
//            Info<< "Reading normalized tower data...";
//            bts.read(reinterpret_cast<char*>(ivals_tow.data()), ivals_tow.size()*sizeof(short));
//            Info<< " done!" << endl;
//
//            Info<< "Dimensionalizing grid data...";
//            i = 0;
//            for(int itime=0; itime<N; ++itime)
//            {
//                for(int j=0; j<Ntower; ++j)
//                {
//                    for(int idim=0; idim<3; ++idim)
//                    {
//                        Utower[i] = (ivals_tow[i] - Vintercept[idim]) / Vslope[idim];
//                        if(U[i] < umin[idim]) umin[idim] = U[i];
//                        if(U[i] > umax[idim]) umax[idim] = U[i];
//                        ++i;
//                    }
//                }
//            }
//            Info<< " done!" << endl;
//        }

        // calculate coordinates
        float y[NY], z[NZ];
        vectorList points(NY*NZ);
        float t[N];

        Info<< "y: " << NY << " [";
        for(int j=0; j<NY; ++j)
        {
            //y = -0.5*(NY-1)*dy + j*dy; // original TurbSim plane
            y[j] = j*dy;

            if((j < 3) || (j >= NY-2))
            {
                Info<< " " << y[j];
            }
            else if(j==3)
            {
                Info << " ...";
            }
        }
        Info<< " ]" << endl;

        Info<< "z: " << NZ << " [";
        for(int k=0; k<NZ; ++k)
        {
            //z[k] = k*dz + zbot; // original TurbSim plane
            z[k] = k*dz;

            if((k < 3) || (k >= NZ-2))
            {
                Info<< " " << z[k];
            }
            else if(k==3)
            {
                Info << " ...";
            }
        }
        Info<< " ]" << endl;

        Info<< "t: " << N << " [";
        for(int i=0; i<N; ++i)
        {
            t[i] = i*dt;

            if((i < 3) || (i >= N-2))
            {
                Info<< " " << t[i];
            }
            else if(i==3)
            {
                Info << " ...";
            }
        }
        Info<< " ]" << endl;

        n = 0;
        for(int j=0; j<NY; ++j)
        {
            for(int k=0; k<NZ; ++k)
            {
                points[n] = vector(0.0, y[j]+dy/2, z[k]+dz/2);
                ++n;
            }
        }

        // DEBUG
        //forAll(patch.Cf(), faceI)
        //{
        //    Info<< patch.Cf()[faceI] << endl;
        //}

        // create list of vector fields (one list item per time)
        // reorder the list to be consistent with column-major ordering
        List<List<vector> > perturb;
        perturb = List<List<vector> >(N, List<vector>(NY*NZ, vector::zero));
        Info<< "Reordering data...";
        Foam::flush(Info);
        label lin_idx;
        label tidx_offset(0);
        for(int itime=0; itime<N; ++itime)
        {
            n = 0;
            for(int j=0; j<NY; ++j)
            {
                for(int k=0; k<NZ; ++k)
                {
                    // do this the naive way to make sure we're getting the indexing right...
                    lin_idx = j + k*NY + tidx_offset;
                    perturb[itime][n] = vector(U[lin_idx], V[lin_idx], W[lin_idx]);
                    ++n;
                }
            }
            tidx_offset += NY*NZ;
        }
        Info<< " done." << endl;

        // calculate statistics
        Info<< "Calculating statistics...";
        Foam::flush(Info);
        float ustd =    0, vstd =    0, wstd =    0;
        float umin =  9e9, vmin =  9e9, wmin =  9e9;
        float umax = -9e9, vmax = -9e9, wmax = -9e9;
        vector Uvec;
        forAll(perturb[0],faceI)
        {
            float uu_sum=0, vv_sum=0, ww_sum=0;
            for(int itime=0; itime<N; ++itime)
            {
                Uvec = perturb[itime][faceI];
                // note: U,V,W are fluctuations
                if(Uvec.x() < umin) umin = Uvec.x();
                if(Uvec.x() > umax) umax = Uvec.x();
                if(Uvec.y() < vmin) vmin = Uvec.y();
                if(Uvec.y() > vmax) vmax = Uvec.y();
                if(Uvec.z() < wmin) wmin = Uvec.z();
                if(Uvec.z() > wmax) wmax = Uvec.z();
                uu_sum += Uvec.x() * Uvec.x();
                vv_sum += Uvec.y() * Uvec.y();
                ww_sum += Uvec.z() * Uvec.z();
            }
            ustd += Foam::sqrt(uu_sum/N);
            vstd += Foam::sqrt(vv_sum/N);
            wstd += Foam::sqrt(ww_sum/N);
        }
        ustd = ustd / (NY*NZ);
        vstd = vstd / (NY*NZ);
        wstd = wstd / (NY*NZ);
        Info<< " done!" << endl;

        Info<< "U min,max,stdev = "
            << umin << " " << umax << " " << ustd << endl;
        Info<< "V min,max,stdev = "
            << vmin << " " << vmax << " " << vstd << endl;
        Info<< "W min,max,stdev = "
            << wmin << " " << wmax << " " << wstd << endl;

    }
    else
    {
        FatalError
            << fname << " not found" << nl
            << exit(FatalError);
    }
}
*/

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
//#   include "createNamedMesh.H"

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

    word pertType(pertDict.lookupOrDefault<word>("syntheticTurbulenceType", "turbsim"));
    List<word> pertPatches(pertDict.lookup("perturbedPatches"));
    Info<< "Patches to perturb: " << pertPatches << endl;

    if(pertType == "turbsim")
    {
/*
        forAll(pertPatches, patchI)
        {
            //label patchID(mesh.boundary().findPatchID(pertPatches[patchI]));
            //Info<< "Patch " << pertPatches[patchI] << " id =" << patchID << endl;

            //Info<< mesh.boundary()[pertPatches[patchI]] << endl;

            word pertFieldName(pertDict.lookup("turbsimField"));
            fileName fpath = runTime.time().constant() / "boundaryData" / pertFieldName+".bts";
            read_bts(fpath);
            //read_bts(fpath, mesh.boundary()[pertPatches[patchI]]);

            break; // DEBUG
        }
*/

        syntheticTurbulence::turbsimBTS ts(runTime);
        //turbsimBTS ts(runTime);
        ts.calcStats();

    }
    else
    {
        Info<< "Synthetic turbulence type not recognized: " << pertType << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
