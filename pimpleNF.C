/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    pimplemfnondimen

Description
    Large time-step transient solver for incompressible, turbulent flow, using
    the PIMPLE (merged PISO-SIMPLE) algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    //dimensionSet::debug = 0;
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
	const labelUList& owner6 = mesh.owner();
	const labelUList& neighbour6 = mesh.neighbour();
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
	    
	    
	    #include "sourcemem.H"
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

//wait for fully developed
	if (mesh.time().value() > 20) 
	{
	   abc=abc+1;
          if (abc==1)
	   {
////initial condition for interfield volume fraction
		forAll (mesh.C(),celli)
        	{
		   if (mesh.C()[celli].y()>0.0032)
		   {
			vf[celli]=4e-5;
			
		   }
		   if (mesh.C()[celli].y()<=0.0032)
		   {
			vf[celli]=0;
			
		   }
		}
	   }
	   
//update diffusion coefficient
 	   #include "dVfcal.H"
	    double vf6=0;
	    surfaceScalarField dVff = fvc::interpolate(dVf);

//all surfaces inside membrane including top, no flux and diffusion
	    forAll(neighbour6, facen)
   	    {
         	if  (srcm[neighbour6[facen]]!=0 || srcm[owner6[facen]]!=0 || srcu[neighbour6[facen]]!=0 || srcu[owner6[facen]]!=0)
		{
     			phi[facen]=0;
			dVff[facen]=0;
			
		}
	    }

//convection diffusion equation
            fvScalarMatrix vfEqn
            (
                fvm::ddt(vf)
              + fvm::div(phi, vf)
              - fvm::laplacian(dVff, vf)
              ==
                fvOptions(vf) 
            );

            vfEqn.relax();
            fvOptions.constrain(vfEqn);
            vfEqn.solve();
            fvOptions.correct(vf);
	    #include "sourceU.H"
	    forAll(mesh.C(),celli)
	   {
		if (vf[celli]>=vf6)
		{
			vf6=vf[celli];
		}
	   }
	   Info<<vf6<<endl;
}
	
	
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
