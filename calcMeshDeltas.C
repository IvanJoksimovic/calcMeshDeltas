/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    calcMeshDeltas

Description

Calculates different mesh metric quantities, i.e. characteristic length scale of the mesh:

- deltaCbrt, length scale based on the cubic root of mech volume
- deltaMax, length scale based on the maximum distance between faces (if mesh is isotropic, this is the same as deltaCbrt)
- deltaMin, length scale based on the minimum distance between faces (if mesh is isotropic, this is the same as deltaCbrt)

call:
calcMeshDelta -latestTime 

Optionally, it cacluates the Kolmogorov length scale, if epsilon is provided

call:
calcMeshDelta -latestTime -epsilonName epsilon




\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculate mesh deltas + optionally calculate Kolmogorov length scale for incompressible flows"
    );

    argList::addOption
    (
        "epsilonName",
        "word",
        "Specify the name of the epsilon field"
    );

    //argList::noParallel();
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField deltaCbrt
    (
        IOobject
        (
            "deltaCbrt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,1,0,0,0,0,0), Zero),
        calculatedFvPatchField<scalar>::typeName
    );

    volScalarField deltaMax
    (
        IOobject
        (
            "deltaMax",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,1,0,0,0,0,0), Zero),
        calculatedFvPatchField<scalar>::typeName
    );

    volScalarField deltaMin
    (
        IOobject
        (
            "deltaMin",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(0,1,0,0,0,0,0), Zero),
        calculatedFvPatchField<scalar>::typeName
    );

    deltaCbrt.ref() = cbrt(mesh.V());

    const cellList& cells = mesh.cells();
    scalarField hmax(cells.size());
    scalarField hmin(cells.size());


    forAll(cells,celli)
    {
        scalar deltaMaxTmp = 0.0;
        scalar deltaMinTmp = 1e18;

        const labelList& cFaces = mesh.cells()[celli];
        const point& centrevector = mesh.cellCentres()[celli];

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const point& facevector = mesh.faceCentres()[facei];
            scalar tmpMax = mag(facevector - centrevector);
            scalar tmpMin = mag(facevector - centrevector);
            if (tmpMax > deltaMaxTmp)
            {
                deltaMaxTmp = tmpMax;
            }
            if (tmpMin < deltaMinTmp)
            {
                deltaMinTmp = tmpMin;
            }

        }

        hmax[celli] = 2*deltaMaxTmp;
        hmin[celli] = 2*deltaMinTmp;
    }

    deltaMax.primitiveFieldRef() = hmax;
    deltaMin.primitiveFieldRef() = hmin;

    Info << "Writting " << deltaCbrt.name() <<" to " << runTime.timeName() << endl;
    deltaCbrt.write();

    Info << "Writting " << deltaMax.name() <<" to " << runTime.timeName() << endl;
    deltaMax.write();

    Info << "Writting " << deltaMin.name() <<" to " << runTime.timeName() << endl;
    deltaMin.write();

    if (!args.found("epsilonName"))
    {
    }
    else
    {
        word epsilonName(args.lookup("epsilonName").toString());

        IOobject header
        (
            epsilonName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        bool foundEpsilon = header.good();// mesh.foundObject<volScalarField>(epsilonName);

        if(!foundEpsilon)
        {
          Info << "Not calculating deltaKolmogorov, no field " << epsilonName <<" found in " << runTime.timeName() << endl;
          Info << mesh.thisDb().sortedNames() << endl;
        }
        else
        {

          Info << "Reading " << epsilonName << endl;

          volScalarField epsilon
          (
              IOobject
              (
                  epsilonName,
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::NO_WRITE
              ),
              mesh
          );

          Info << "Reading viscosity" << endl;

          IOdictionary transportProperties
          (
              IOobject
              (
                  "transportProperties",
                  runTime.constant(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::NO_WRITE
              )
          );

          dimensionedScalar nu
          (
              "nu",
              dimViscosity,
              transportProperties
          );

          if(nu.value() == 0)
          {
            Info << "ERROR: Read constant zero viscosity, viscosity must be defined in constant/transportProperties" << endl;
            Info << "Terminating program!" << endl;
            return 0;
          }

          Info << "Calculating deltaKolmogorov" << endl;
          volScalarField deltaKolmogorov
          (
              IOobject
              (
                  "deltaKolmogorov",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              pow(pow(nu,3)/epsilon,0.25)
          );

          //deltaKolmogorov = pow(pow(nu,3)/epsilon,0.25);
          Info << "Writting " << deltaKolmogorov.name() <<" to " << runTime.timeName() << endl;
          deltaKolmogorov.write();
        }


    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
