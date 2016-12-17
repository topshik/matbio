/****************************************************************************************************************
* This is a collection of procedures for obtaining correlation functions from spatial patterns			*
* 														*
* Last update 4.4.98.												*
****************************************************************************************************************/


/****************************************************************************************************************
* Procedure for calculating mean path of total population sizes of stochastic realizations.			*
* This procedure was developed in consultation with Ulf.  The idea is to achieve a path without having to	*
* do any averaging of the realizations (and without having to hold large arrays).				*
****************************************************************************************************************/
void	mean_path_tpss()
{
	FILE	**in;			/* the file pointers are to be held in a dynamic array			*/

	float	time[niterate];		/* current time in each realization					*/
	int	nnow[niterate][nspp];	/* current   # individuals of each sp in each realization		*/
	int	nold[niterate][nspp];	/* previous  # individuals of each sp in each realization		*/
	float	ndelta[nspp];		/* change in # individuals of each sp in realization being updated 	*/
	float	npath[nspp];		/* current   # individuals of each sp in mean path			*/

	char	filename[40];
	int	sp, iterate, choose;
	float	tnow;


	in   = (FILE**) malloc (niterate * sizeof(FILE*));
	out6 = fopen("output.tpssmean.dat", "w");

    /***Open the tpss output files***/
	for (iterate=0; iterate<niterate; ++iterate)
	{	sprintf(filename, "output.%d.tpss.dat", iterate);
		in[iterate] = fopen(filename, "r");
	}

    /***Read the first row of each file (time 0)***/
	for (iterate=0; iterate<niterate; ++iterate)
	{	fscanf(in[iterate], "%g", &time[iterate]);
		for (sp=1; sp<nspp; ++sp)	fscanf(in[iterate], "%4d", &nold[iterate][sp]);
	}

    /***Set the initial values of mean path (Can do this from any iterate as all start the same)***/
	for (sp=1; sp<nspp; ++sp)	npath[sp] = (float) nold[0][sp];

    /***Record initial state of mean path***/
	fprintf(out6, "  %8.4f  ", time[0]);
	for (sp=1; sp<nspp; ++sp)	fprintf(out6, "%10.4f", npath[sp]);
	fprintf(out6, "\n");
	fflush(out6);

    /***Read each file to find the first event***/
	for (iterate=0; iterate<niterate; ++iterate)
	{	fscanf(in[iterate], "%g", &time[iterate]);
		for (sp=1; sp<nspp; ++sp)	fscanf(in[iterate], "%4d", &nnow[iterate][sp]);
	}

    /***Step through all the realizations picking the next event***/
	while(tnow<tmax)
	{
	    /***Find the time and realization in which the next event happens***/
		tnow = 1000.0;
		for (iterate=0; iterate<niterate; ++iterate)
		if (time[iterate] < tnow)
		{	tnow   = time[iterate];
			choose = iterate;
		}

	    /***Obtain the new number of individuals in mean path.  We take the change in***/
	    /***number of this path and divide by the number of realizations***/
	    /***This is the only averaging involved***/
		for (sp=1; sp<nspp; ++sp)	ndelta[sp]         = ((float)nnow[choose][sp] - (float)nold[choose][sp]) / (float)niterate;
		for (sp=1; sp<nspp; ++sp)	npath[sp]          =  npath[sp] + ndelta[sp];
		for (sp=1; sp<nspp; ++sp)	nold[choose][sp]   =  nnow[choose][sp];
	
	    /***Read the next line of the file of the chosen realization***/
		fscanf(in[choose], "%g", &time[choose]);
		for (sp=1; sp<nspp; ++sp)	fscanf(in[choose], "%4d", &nnow[choose][sp]);

	    /***Record the new value of the mean path***/
/*		fprintf(out6, "%d ", choose); 
*/		fprintf(out6, "%8.4f  ", tnow);
		for (sp=1; sp<nspp; ++sp)	fprintf(out6, "%10.4f", npath[sp]);
		fprintf(out6, "\n");
	}

    /***Close the tpss output files***/
	free(in);
	fflush(out6);
	fclose(out6);
}


/****************************************************************************************************************
* Procedure for calculating mean path of radial correlation functions of stochastic realizations.		*
****************************************************************************************************************/
void	mean_path_rcfs()
{
	FILE	**in;			/* the file pointers are to be held in a dynamic array			*/

	float	rcfsstore[niterate][binmax][1+(nspp-1)*(nspp-1)];	/*Array holding rcfs (all iterations) at current time	*/
	float	rcfsmean[binmax][1+(nspp-1)*(nspp-1)];			/*Array holding rcfs means (across iterations		*/

	int	iterate, bin, sp;
	char	filename[40];

	in   = (FILE**) malloc (niterate * sizeof(FILE*));
	out7 = fopen("output.rcfsmean.dat", "w");

    /***Open the rcfs output files for reading***/
	for (iterate=0; iterate<niterate; ++iterate)
	{	sprintf(filename, "output.%d.rcfs.dat", iterate);
		in[iterate] = fopen(filename, "r");
	}

    /***Read the first rcfs of each file (time 0)***/
	for (iterate=0; iterate<niterate; ++iterate)
	for (bin=0; bin<binmax; ++bin)
	{	fscanf(in[iterate], "%g", &t);
		for (sp=0; sp<1+(nspp-1)*(nspp-1); ++sp)	fscanf(in[iterate], "%g", &rcfsstore[iterate][bin][sp]);
	}

    /***Set rcfsmean to zero***/
	for (sp=0; sp<1+(nspp-1)*(nspp-1); ++sp)
	for (bin=0; bin<binmax; ++bin)
	for (iterate=0; iterate<niterate; ++iterate)		rcfsmean[bin][sp] = 0.0;

    /***Calculate the mean values of rcfs***/
	for (sp=0; sp<1+(nspp-1)*(nspp-1); ++sp)
	for (bin=0; bin<binmax; ++bin)
	for (iterate=0; iterate<niterate; ++iterate)		rcfsmean[bin][sp] = rcfs[bin][sp] + rcfsstore[iterate][bin][sp];
	for (sp=0; sp<1+(nspp-1)*(nspp-1); ++sp)
	for (bin=0; bin<binmax; ++bin)				rcfsmean[bin][sp] = rcfs[bin][sp] / niterate;

    /***Output the radial correlation functions***/
	for (bin=0; bin<binmax; ++bin)
	{	for (sp=0; sp<1+(nspp-1)*(nspp-1); ++sp)	fprintf(out7, "%8.4f", rcfsmean[bin][sp]);
		fprintf(out7, "\n");
	}
	fprintf(out7, "\n");

    /***Close the rcfs output files***/
	free(in);
	fflush(out7);
	fclose(out7);
}


/****************************************************************************************************************
* Procedure for calculating spatial correlation functions							*
****************************************************************************************************************/
void	spatial_correlation_functions()
{
	repP	targP;		/* pointer to the target individual			*/
	repP	nbrP;		/* pointer to a neighbour of target individual (targP)	*/
	float	x, y;		/* x and y coordinates of neighbour relative to target	*/
	float	distance;	/* distance from neighbour to target			*/
	int	xbin, ybin, i, j, k;
	int	targsp, nbrsp;
	float	scale;


    /***Initialize the spatial correlation functions***/
	for (xbin=0;  xbin < 2*binmax;         ++xbin)
	for (ybin=0;  ybin < 2*binmax;         ++ybin)
	for (i=0;     i    < 1+(nspp-1)*(nspp-1); ++i)		scfs[xbin][ybin][i] = 0.0;

    /***Loop for target individuals***/
	for (targsp=1; targsp<nspp; ++targsp)
	if (param[targsp].n > 0)
	{
		targP = param[targsp].firstP;
		while (targP != NULL)
		{
		    /***Loop for neighbours around a target individual***/
			for (nbrsp=1; nbrsp<nspp; ++nbrsp)
			if (param[nbrsp].n > 0)
			{	
				nbrP = param[nbrsp].firstP;
				while (nbrP != NULL)
				{
				    /***Target individual itself is not counted as a neighbour!***/
					if (nbrP != targP)
					{
					    /***Calculate the location of neighbour relative to i***/
						neighbour_distance(targP, nbrP, rcfsmax, &x, &y, &distance);
						
					    /***Column of spatial correlation function array to which targP and nbrP refer***/
						i = (targP->sp - 1) * (nspp - 1) + nbrP->sp;

					    /***Bin in which to place the neighbour***/
						xbin = binmax + x/bininc;
						ybin = binmax + y/bininc;
						if (xbin>=0 && xbin < 2*binmax && ybin>=0 && ybin < 2*binmax)
						{
							scfs[xbin][ybin][i] = scfs[xbin][ybin][i] + 1.0;

						    /***Detailed record of xcfs calculation***/
							if (scfs_flag == 1)
							{	fprintf(out4, "targsp:%d  nbrsp:%d  i:%d  targP:%d  nbrP:%d  x:%7.4f  y:%7.4f  distance:%6.4f", 
								targsp, nbrsp, i, targP, nbrP, x, y, distance);
								fprintf(out4, "   scfs[%2d][%2d][%d]:%6.3f\n", xbin-binmax, ybin-binmax, i, scfs[xbin][ybin][i]);
								fflush(out4);
							}
						}
					}
					nbrP = nbrP->nextP;
				}
			}
			if (scfs_flag ==1)	fprintf(out4, "\n");
			targP = targP->nextP;
		}
	}

    /***Record the spatial correlation functions before normalisation***/
	if (scfs_flag ==1)
	for (i=1; i<1+(nspp-1)*(nspp-1); ++i)
	{	fprintf(out4, "Spatial correlation function %d before normalization\n", i);
		for (xbin=0; xbin<2*binmax; ++xbin)
		{	for (ybin=0; ybin<2*binmax; ++ybin)
			fprintf(out4, "%10.4f", scfs[xbin][ybin][i]);
			fprintf(out4, "\n");
		}
	}
	fflush(out4);

    /***Normalise the spatial correlation functions***/
	if (scfs_flag ==1)	fprintf(out4, "\nScaling factors for spatial correlation functions\n");
	for (i=1; i<nspp; ++i)
	for (j=1; j<nspp; ++j)
	{
	    /***Scaling to allow for area of cell encompassed by the present bin***/
		scale = 1.0 / (bininc*bininc);

	    /***Scaling to allow for density  of species i and j (scfs relative to mean field expectation)***/
		if (i>0 && j>0)		scale = scale * arena_area * arena_area / (float) (param[i].n * param[j].n);
		else			scale = 0.0;

	    /***Record of scaling factor***/
		if (scfs_flag ==1)
		fprintf(out4, "xbin:%d  ybin:%d  k:%d  param[%d].n:%d  param[%d].n:%d  scale:%6.4f\n", 
			xbin, ybin, k, i, param[i].n, j, param[j].n, scale, k);

		for (xbin=0; xbin<2*binmax; ++xbin)
		for (ybin=0; ybin<2*binmax; ++ybin)
		{
		    /***Multiply scfs coefficient by scaling factor***/
			k = (i - 1) * (nspp - 1) + j;
			scfs[xbin][ybin][k] = scfs[xbin][ybin][k] * scale;
		}
	}

    /***Record the spatial correlation functions after normalization***/
	if (scfs_flag ==1)
	for (i=1; i<1+(nspp-1)*(nspp-1); ++i)
	{	fprintf(out4, "Spatial correlation function %d after normalization\n", i);
		for (xbin=0; xbin<2*binmax; ++xbin)
		{	for (ybin=0; ybin<2*binmax; ++ybin)
			fprintf(out4, "%10.4f", scfs[xbin][ybin][i]);
			fprintf(out4, "\n");
		}
	}
	fflush(out4);

    /***Record the spatial correlation functions (keeps only the current scfs array)***/
	rewind(out5);
	if (2.0*rcfsmax > xmax-xmin)
	fprintf(out5, "\nWARNING:  maximum rcfs exceeds arena length / 2. This will generate errors\n");
	fprintf(out5, "#Spatial correlation functions at time:%8.4f\n", t);
	for (xbin=0; xbin<2*binmax; ++xbin)
	{	for (ybin=0; ybin<2*binmax; ++ybin)
		{	fprintf(out5, "%10.4f", bininc * ((float)xbin - (float)binmax));
			fprintf(out5, "%10.4f", bininc * ((float)ybin - (float)binmax));
			for (i=1; i<1+(nspp-1)*(nspp-1); ++i)	fprintf(out5, "%10.4f", scfs[xbin][ybin][i]);
			fprintf(out5, "\n");
		}
		fprintf(out5, "\n");
	}
	fflush(out5);
}


/****************************************************************************************************************
* Procedure for calculating radial correlation functions							*
****************************************************************************************************************/
void	radial_correlation_functions()
{
	repP	targP;		/* pointer to the target individual			*/
	repP	nbrP;		/* pointer to a neighbour of target individual (targP)	*/
	float	x, y;		/* x and y coordinates of neighbour relative to target	*/
	float	distance;	/* distance from neighbour to target			*/
	int	bin, i, j, k;
	int	targsp, nbrsp;
	float	scale, r1, r2;


    /***Initialize the radial correlation functions***/
	for (bin=0;  bin < binmax;          ++bin)
	for (i=0;    i   < 1+(nspp-1)*(nspp-1); ++i)		rcfs[bin][i] = 0.0;
	for (bin=1;  bin < binmax; ++bin)			rcfs[bin][0] = rcfs[bin-1][0] + bininc;

    /***Loop for target individuals***/
	for (targsp=1; targsp<nspp; ++targsp)
	{	
		targP = param[targsp].firstP;
		while (targP != NULL)
		{	
		    /***Loop for neighbours around a target individual***/
			for (nbrsp=1; nbrsp<nspp; ++nbrsp)
			{	
				nbrP = param[nbrsp].firstP;
				while (nbrP != NULL)
				{
				    /***Target individual itself is not counted as a neighbour!***/
					if (nbrP != targP)
					{
					    /***Calculate the distance from target individual to neighbour***/
						neighbour_distance(targP, nbrP, rcfsmax, &x, &y, &distance);
						
					    /***Column of radial correlation function array to which targP and nbrP refer***/
						i = (targP->sp - 1) * (nspp - 1) + nbrP->sp;

					    /***Bin in which to place the neighbour***/
						bin = distance/bininc;
						if (bin < binmax)
						{
							rcfs[bin][i] = rcfs[bin][i] + 1.0;

						    /***Detailed record of rcfs calculation***/
							if (rcfs_flag == 1)
							{	fprintf(out4, "targsp:%d  nbrsp:%d  i:%d  targP:%d  nbrP:%d  distance:%6.4f", 
								targsp, nbrsp, i, targP, nbrP, distance);
								fprintf(out4, "   rcfs[%2d][%d]:%6.3f\n", bin, i, rcfs[bin][i]);
								fflush(out4);
							}
						}
					}
					nbrP = nbrP->nextP;
				}
			}
			if (rcfs_flag ==1)	fprintf(out4, "\n");
			targP = targP->nextP;
		}
	}

    /***Record the radial correlation functions before normalisation***/
	if (rcfs_flag ==1)
	{	fprintf(out4, "Radial correlation function before normalisation\n");
		for (bin=0; bin<binmax; ++bin)
		{	for (i=0; i<1+(nspp-1)*(nspp-1); ++i)	fprintf(out4, "%10.4f", rcfs[bin][i]);
			fprintf(out4, "\n");
			if (bin == binmax - 1)	fprintf(out4, "\n");
		}
	}
	fflush(out4);

    /***Normalise the radial correlation functions***/
	for (i=1; i<nspp; ++i)
	for (j=1; j<nspp; ++j)
	{
		for (bin=0; bin < binmax; ++bin)
		{
		    /***Scaling to all for area of ring encompassed by the present bin***/
			r1    = bininc * (float) bin;
			r2    = bininc * (float) (bin + 1);
			scale = 1.0/(pi*(r2*r2 - r1*r1));

		    /***Scaling to allow for density of species i and j (rcfs relative to mean field expectation)***/
			if (i>0 && j>0)		scale = scale * arena_area * arena_area /(float) (param[i].n * param[j].n);
			else			scale = 0.0;

		    /***Multiply rcfs coefficient by scaling factor***/
			k = (i - 1) * (nspp - 1) + j;
			rcfs[bin][k] = rcfs[bin][k] * scale;

		    /***Record of scaling factor***/
			if (rcfs_flag ==1)
			fprintf(out4, "bin:%d  i:%d  j:%d  k:%d  r1:%6.3f  r2:%6.3f param[%d].n:%d  param[%d].n:%d  scale:%6.3f\n", 
			bin, i, j, k, r1, r2, i, param[i].n, j, param[j].n, scale, k);
		}
		if (rcfs_flag == 1)	fprintf(out4, "\n");
	}

    /***Record the radial correlation functions***/
	if (2.0*rcfsmax > xmax-xmin)
	fprintf(out3, "\nWARNING:  maximum rcfs exceeds arena length / 2. This will generate errors\n");
/*	fprintf(out3, "#Radial correlation functions at time:%8.4f\n", t);
*/	for (bin=0; bin<binmax; ++bin)
	{	fprintf(out3, "%10.4f", t);
		fprintf(out3, "%10.4f", rcfs[bin][0]+0.5*bininc);
		for (i=1; i<1+(nspp-1)*(nspp-1); ++i)	fprintf(out3, "%10.4f", rcfs[bin][i]);
		fprintf(out3, "\n");
	}
	fprintf(out3, "\n");
	fflush(out3);
}

/****************************************************************************************************************
* Procedure for calculating third moment functions								*
* These are given as functions of two radii: the distances i->j and i->k					*
****************************************************************************************************************/
void	third_moment_functions()
{
	repP	targP;		/* pointer to the target individual i			*/
	repP	jnbrP;		/* pointer to neighbour j of target individual i	*/
	repP	knbrP;		/* pointer to neighbour k of target individual i	*/
	float	x, y;		/* x and y coordinates of neighbour relative to target	*/
	float	jdistance;	/* distance from a j neighbour to target		*/
	float	kdistance;	/* distance from a k neighbour to target		*/
	int	jbin, kbin, ii, i, j, k;
	int	targsp, jnbrsp, knbrsp;
	float	scale, jscale, kscale, jr1, jr2, kr1, kr2;


    /***Initialize the third moments***/
	for (jbin=0;  jbin < binmax;            ++jbin)
	for (kbin=0;  kbin < binmax;            ++kbin)
	for (i=0;        i < 1+(nspp-1)*(nspp-1);  ++i)		mom3[jbin][kbin][i] = 0.0;
	for (jbin=1;  jbin < binmax;            ++jbin)		mom3[jbin][0][0] = mom3[jbin-1][0][0] + bininc;
	for (kbin=1;  kbin < binmax;            ++kbin)		mom3[0][kbin][0] = mom3[0][kbin-1][0] + bininc;

    /***Loop for target i individuals***/
	for (targsp=1; targsp<nspp; ++targsp)
	{	
		targP = param[targsp].firstP;
		while (targP != NULL)
		{	
		    /***Loop for neighbours j around i***/
			for (jnbrsp=1; jnbrsp<nspp; ++jnbrsp)
			{	
				jnbrP = param[jnbrsp].firstP;
				while (jnbrP != NULL)
				{
				    /***Target individual itself is not counted as a neighbour!***/
					if (jnbrP != targP)
					{
					    /***Loop for neighbours k around i***/
						for (knbrsp=1; knbrsp<nspp; ++knbrsp)
						{
							knbrP = param[knbrsp].firstP;
							while (knbrP != NULL)
							{
							    /***Target individual itself and neighbour j are not counted as a neighbours!***/
								if (knbrP != targP && knbrP != jnbrP)
								{
								    /***Calculate the distance from target individual to neighbour j and k***/
									neighbour_distance(targP, jnbrP, rcfsmax, &x, &y, &jdistance);
									neighbour_distance(targP, knbrP, rcfsmax, &x, &y, &kdistance);

								    /***Column of third moment array to which targP, jnbrP and knbrP refer***/
									ii = (targP->sp - 1) * (nspp - 1) + (jnbrP->sp - 1) * (nspp - 1) + knbrP->sp;

								    /***bins for neighbour j and k***/
									jbin = jdistance/bininc;
									kbin = kdistance/bininc;
									if (jbin < binmax && kbin < binmax)
									{
										mom3[jbin][kbin][ii] = mom3[jbin][kbin][ii] + 1.0;

				    /***Detailed record of third moment calculation***/
					if (mom3_flag == 1)
					{	fprintf(out4, "targsp:%d  jnbrsp:%d  knbrsp:%d  ii:%d  targP:%d  jnbrP:%d  knbrP:%d  jdistance:%6.4f  kdistance:%6.4f", 
								targsp, jnbrsp, knbrsp, ii, targP, jnbrP, knbrP, jdistance, kdistance);
						fprintf(out4, "   mom3[%d][%d][%d]:%6.3f\n", jbin, kbin, ii, mom3[jbin][kbin][ii]);
						fflush(out4);
					}
									}
								}
								knbrP = knbrP->nextP;
							}
						}
					}
					jnbrP = jnbrP->nextP;
				}
			}
			if (mom3_flag ==1)	fprintf(out4, "\n");
			targP = targP->nextP;
		}
	}

    /***Normalise the third moments***/
	for (i=1; i<nspp; ++i)
	for (j=1; j<nspp; ++j)
	for (k=1; k<nspp; ++k)
	{
		for (jbin=0; jbin < binmax; ++jbin)
		for (kbin=0; kbin < binmax; ++kbin)
		{
		    /***Scaling to allow for area of rings encompassed by the present jbin and kbin***/
			jr1    = bininc * (float)  jbin;
			jr2    = bininc * (float) (jbin + 1);
			jscale = 1.0 / (pi*(jr2*jr2 - jr1*jr1));
			kr1    = bininc * (float)  kbin;
			kr2    = bininc * (float) (kbin + 1);
			kscale = 1.0 / (pi*(kr2*kr2 - kr1*kr1));
			scale  = jscale * kscale;

		    /***Scaling to allow for density of species i, j and k (third moment relative to mean field expectation)***/
			scale = scale * arena_area * arena_area * arena_area /(float) (param[i].n * param[j].n * param[k].n);

		    /***Multiply  third moment by scaling factor***/
			ii = (i - 1) * (nspp - 1) + (j - 1) * (nspp - 1) + k;
			mom3[jbin][kbin][ii] = mom3[jbin][kbin][ii] * scale;

		    /***Record of scaling factor***/
			if (mom3_flag ==1)
			fprintf(out4, "jbin:%d  kbin:%d  ii:%d  ",   jbin, kbin, ii);
			fprintf(out4, "jr1:%6.3f  jr2:%6.3f jscale:%6.3f  kr1:%6.3f  kr2:%6.3f  kscale:%6.3f  ", jr1, jr2, jscale, kr1, kr2, kscale);
			fprintf(out4, "n%d:%d  n%d:%d  n%d:%d  scale:%6.3f\n", i, param[i].n, j, param[j].n, k, param[k].n, scale);
		}
		if (mom3_flag == 1)	fprintf(out4, "\n");
	}

    /***Record the third moments as a table jbin by kbin***/
/*	fprintf(out8, "#Third moment functions at time:%8.4f\n", t);
	for (i=1; i<nspp; ++i)
	for (j=1; j<nspp; ++j)
	for (k=1; k<nspp; ++k)
	{
		fprintf(out8, "#Species combination i:%d j:%d k:%d \n", i, j, k);
		ii = (i - 1) * (nspp - 1) + (j - 1) * (nspp - 1) + k;

		fprintf(out8, "          ");
		for (kbin=0; kbin<binmax; ++kbin)	fprintf(out8, "%10.4f", mom3[0][kbin][0]);
		fprintf(out8, "\n");

		for (jbin=0; jbin<binmax; ++jbin)
		{	fprintf(out8, "%10.4f", mom3[jbin][0][0]);
			for (kbin=0; kbin<binmax; ++kbin)	fprintf(out8, "%10.4f", mom3[jbin][kbin][i]);
			fprintf(out8, "\n");
		}
		fprintf(out8, "\n");
		fflush(out8);
	}
*/
    /***Record the third moments as columns for gnuplot***/
/*	fprintf(out8, "#Third moment functions at time:%8.4f\n", t);
*/	for (jbin=0; jbin<binmax; ++jbin)
	{	for (kbin=0; kbin<binmax; ++kbin)
		{	fprintf(out8, "%10.4f", mom3[jbin][0][0]+0.5*bininc);
			fprintf(out8, "%10.4f", mom3[0][kbin][0]+0.5*bininc);
			for (i=1;  i<1+(nspp - 1)*(nspp - 1)*(nspp - 1); ++i)	fprintf(out8, "%10.4f", mom3[jbin][kbin][i]);
			fprintf(out8, "\n");
			fflush(out8);
		}
		fprintf(out8, "\n");
		
	}
}
