// define the source term for membrane in nanofiltration
int fg2=0;
double Rm=8.11e+11;
double thi=0.0012;
int numb=0;
const volScalarField nu2 = turbulence->nu();
//find cake thickness and Rc

forAll(mesh.C(),cellh)
{
	if(vf[cellh]>0.63)
	{
		ymax[cellh]=mesh.C()[cellh].y();
		forAll(mesh.C(),cellj)
		{
			if ( abs(mesh.C()[cellj].x()-mesh.C()[cellh].x())<1e-8 && vf[cellj]>0.63 && mesh.C()[cellj].y()>ymax[cellh])
			{
				ymax[cellh]=mesh.C()[cellj].y();
			}
		}
		rc[cellh]=((ymax[cellh]-0.0032)/(0.75e-6)*2+1)/2*(0.75e-6)*(1.6e+17);
	}
}
forAll (mesh.C(),celli)
{
	   if (mesh.C()[celli].x()>=(0.02)  && mesh.C()[celli].x()<=(0.07) && mesh.C()[celli].y()>(0.002) && mesh.C()[celli].y()<(0.0032)) 
	   {
		fg2=1;
		srcm[celli]=fg2*nu2[celli]*Rm/thi;
//if with same x position there is a cake, rm'=rm+rc'
		forAll(mesh.C(),cellm)
		{
			if ( abs(mesh.C()[cellm].x()-mesh.C()[celli].x())<1e-8 && vf[cellm]>0.63)
			{
				srcm[celli]=fg2*nu2[celli]*Rm/thi+fg2*nu2[celli]*rc[cellm]/(thi+ymax[cellm]-0.0032);
			}
		}
	   }
		
}
