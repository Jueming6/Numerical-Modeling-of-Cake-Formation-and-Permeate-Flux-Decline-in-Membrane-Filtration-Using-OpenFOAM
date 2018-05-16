// define the source term for U in microfiltration
int fg=0;
const volScalarField nu1 = turbulence->nu();

forAll (mesh.C(),celli)
{

	   if (vf[celli]>=0.6)
	   {
		fg=1;
		srcu[celli]=fg*nu1[celli]/nu1[celli]*5.0*pow((1-0.4),2)*pow((3/a.value()),2)/pow(0.4,3);
	   }
	   
	   fg=0;

}


forAll (mesh.boundaryMesh(), patchI)
{
    	   forAll(srcu.boundaryField()[patchI], faceI)
           {
		
	        if (vf.boundaryField()[patchI][faceI]>=0.6)
	        {
		    fg=1;
		    srcu.boundaryFieldRef()[patchI][faceI]=fg*nu1.boundaryField()[patchI][faceI]/nu1.boundaryField()[patchI][faceI]*5.0*pow((1-0.4),2)*pow(3/a.value(),2)/pow(0.4,3);
	        }
	   
	   fg=0;
	   }
        
}


