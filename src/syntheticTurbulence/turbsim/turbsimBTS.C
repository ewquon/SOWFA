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

#include "turbsimBTS.H"

namespace Foam
{

namespace syntheticTurbulence
{

// * * * * * * * * * * * * * *  Constructor  * * * * * * * * * * * * * * * * //

turbsimBTS::turbsimBTS
(
    const Time& runTime
)
:
    perturbations(runTime)
{
    word fieldName(perturbDict_.lookup("turbsimField"));
    fileName fpath = runTime_.time().constant() / "boundaryData" / fieldName+".bts";
    read(fpath);
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void turbsimBTS::read(word fname)
//void read(word fname, const fvPatch &patch)
{
    std::ifstream bts(fname.c_str(), std::ios::binary);
    if(bts)
    {
        bts.seekg(0, std::ios_base::end);
        Info<< "Reading " << fname
            << " (" << bts.tellg() << " bytes)"
            << endl;
        bts.seekg(0, std::ios_base::beg);

        // TurbSim IO vars
        short id;
        float uhub;
        float Vslope[3];
        float Vintercept[3];

        // unused variables from TurbSim
        float zhub, zbot;
        int Ntower;

        // working vars
        int idx; // flattened index

        bts.read(reinterpret_cast<char*>(&id), sizeof(id));
        //assert( id==7 || id==8 )
        if(id==7)
        {
            Info << "non-periodic" << endl;
        }
        else if (id==8)
        {
            Info << "periodic" << endl;
            periodic = true;
        }
        else
        {
            Info << "Warning: Unknown BTS file ID: " << id << endl;
        }

        // read resolution settings
        bts.read(reinterpret_cast<char*>(&Nz), sizeof(Nz));
        bts.read(reinterpret_cast<char*>(&Ny), sizeof(Ny));
        Info<< "NY, NZ = " << Ny << " " << Nz << endl;

        bts.read(reinterpret_cast<char*>(&Ntower), sizeof(Ntower));
        Info<< "Ntower = " << Ntower << endl;

        bts.read(reinterpret_cast<char*>(&Nt), sizeof(Nt));
        Info<< "N = " << Nt << endl;

        bts.read(reinterpret_cast<char*>(&dy), sizeof(dy));
        bts.read(reinterpret_cast<char*>(&dz), sizeof(dz));
        bts.read(reinterpret_cast<char*>(&dt), sizeof(dt));
        Info<< "dy, dz = " << dy << " " << dz << endl;
        Info<< "dt = " << dt << endl;

        period = Nt * dt;
        Info<< "period = " << period << endl;
        if(!periodic) period *= -1;

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
        std::vector<short> ivals(3*Ny*Nz*Nt);
        std::vector<float> U(Ny*Nz*Nt);
        std::vector<float> V(Ny*Nz*Nt);
        std::vector<float> W(Ny*Nz*Nt);

        Info<< "Reading normalized velocity data...";
        Foam::flush(Info);
        bts.read(reinterpret_cast<char*>(ivals.data()), ivals.size()*sizeof(short));
        Info<< " done!" << endl;

        Info<< "Dimensionalizing velocity data...";
        int icomp = 0; // component-wise flattened index
        idx = 0; // full flattened index
        for(int itime=0; itime<Nt; ++itime)
        {
            for(int k=0; k<Nz; ++k)
            {
                for(int j=0; j<Ny; ++j)
                {
                    U[icomp] = (ivals[idx] - Vintercept[0]) / Vslope[0] - uhub;
                    V[icomp] = (ivals[idx+1] - Vintercept[1]) / Vslope[1];
                    W[icomp] = (ivals[idx+2] - Vintercept[2]) / Vslope[2];
                    idx += 3;
                    ++icomp;
                }
            }
        }
        Info<< " done!" << endl;

        // calculate inflow plane y-coordinates
        float y[Ny];
        Info<< "y: " << Ny << " [";
        for(int j=0; j<Ny; ++j)
        {
            //y = -0.5*(NY-1)*dy + j*dy; // original TurbSim plane
            y[j] = j*dy - dy/2;

            if((j < 3) || (j >= Ny-2))
            {
                Info<< " " << y[j];
            }
            else if(j==3)
            {
                Info << " ...";
            }
        }
        Info<< " ]" << endl;

        // calculate inflow plane z-coordinates
        float z[Nz];
        Info<< "z: " << Nz << " [";
        for(int k=0; k<Nz; ++k)
        {
            //z[k] = k*dz + zbot; // original TurbSim plane
            z[k] = k*dz - dz/2;

            if((k < 3) || (k >= Nz-2))
            {
                Info<< " " << z[k];
            }
            else if(k==3)
            {
                Info << " ...";
            }
        }
        Info<< " ]" << endl;

        // calculate inflow time vector
        t.setSize(Nt);
        Info<< "t: " << Nt << " [";
        for(int i=0; i<Nt; ++i)
        {
            t[i] = i*dt;

            if((i < 3) || (i >= Nt-2))
            {
                Info<< " " << t[i];
            }
            else if(i==3)
            {
                Info << " ...";
            }
        }
        Info<< " ]" << endl;

        // setup list of inflow perturbation points
        points.setSize(Ny*Nz);
        idx = 0;
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                points[idx] = vector(0.0, y[j], z[k]);
                ++idx;
            }
        }

        // Create list of vector fields and reorder the list to have column-
        // major ordering consistent with OpenFOAM
        //perturb = List<List<vector> >(N, List<vector>(NY*NZ, vector::zero));
        Info<< "Reordering data...";
        Foam::flush(Info);
        label tidx_offset(0);
        for(int itime=0; itime<Nt; ++itime)
        {
            List<vector> tempPlane(Ny*Nz, vector::zero);
            int n = 0;
            for(int j=0; j<Ny; ++j)
            {
                for(int k=0; k<Nz; ++k)
                {
                    // do this the naive way to make sure we get the indexing right...
                    idx = j + k*Ny + tidx_offset;
                    tempPlane[n] = vector(U[idx], V[idx], W[idx]);
                    ++n;
                }
            }
            tidx_offset += Ny*Nz;
            perturb.append(tempPlane);
        }
        Info<< " done." << endl;

        // skip reading turbsim tower points

    }
    else
    {
        FatalError
            << fname << " not found" << nl
            << exit(FatalError);
    }
}

void turbsimBTS::calcStats()
{
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
        for(int itime=0; itime<Nt; ++itime)
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
        ustd += Foam::sqrt(uu_sum/Nt);
        vstd += Foam::sqrt(vv_sum/Nt);
        wstd += Foam::sqrt(ww_sum/Nt);
    }
    ustd = ustd / (Ny*Nz);
    vstd = vstd / (Ny*Nz);
    wstd = wstd / (Ny*Nz);
    Info<< " done!" << endl;

    Info<< "U min,max,stdev = "
        << umin << " " << umax << " " << ustd << endl;
    Info<< "V min,max,stdev = "
        << vmin << " " << vmax << " " << vstd << endl;
    Info<< "W min,max,stdev = "
        << wmin << " " << wmax << " " << wstd << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam

