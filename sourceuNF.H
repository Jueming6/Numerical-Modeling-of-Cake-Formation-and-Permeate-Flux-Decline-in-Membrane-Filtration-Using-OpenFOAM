// define the source term for U in nanofiltration
const volScalarField nu1 = turbulence->nu();

//find cake thickness,Rc,rc'
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
		srcu[cellh]=nu1[cellh]*rc[cellh]/(0.0012+ymax[cellh]-0.0032);
	}
}

