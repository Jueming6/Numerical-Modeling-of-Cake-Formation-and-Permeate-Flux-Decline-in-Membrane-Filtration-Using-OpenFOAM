// define the source term for membrane in microfiltration
int fg2=0;
double Rm=1.62e+2;
double thi=1.2;
const volScalarField nu2 = turbulence->nu();

forAll (mesh.C(),celli)
{
//find membrane region
	   if (mesh.C()[celli].x()>=(15)  && mesh.C()[celli].x()<=(45) && mesh.C()[celli].y()>(10) && mesh.C()[celli].y()<(11.2)) 
	   {
		fg2=1;

	   }
           srcm[celli]=fg2*nu2[celli]*Rm/thi;
	   fg2=0;
}


