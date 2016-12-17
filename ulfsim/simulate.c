/****************************************************************************************************************
* This program carries out a realisation of a multspecies community in continuous space and time.		*
* 														*
* Last update 7.4.98.												*
****************************************************************************************************************/

#include	<stdlib.h>
#include	<math.h>
#include	<stdio.h>
#include	<string.h>
#include	"vogl.h"
#include	"vodevice.h"


#define		NULL	0
#define		error	1e-9		/* rounding error allowed for type float			*/
#define		pi	3.1415927


#define		xmin		0.0	/* lower x-coordinate for arena					*/
#define		xmax		1.0	/* upper x-coordinate for arena					*/
#define		ymin		0.0	/* lower y-coordinate for arena					*/
#define		ymax		1.0	/* upper y-coordinate for arena					*/


#define		nspp		3	/* Number of species + 1					*/
#define		niterate	1	/* Number of realizations from a given initial state		*/
#define		grid		0	/* 0: realizations from 1 starting point	1: from grid	*/
					/* only to be used when nspp = 3				*/
#define		n1start		100	/* Parameters for setting the initial population sizes for grid */
#define		n1stop		200	/* Only apply if grid = 1					*/
#define		n1inc		50
#define		n2start		100
#define		n2stop		200
#define		n2inc		50


#define		tmax		0.0	/* maximum time for simulation					*/
#define		maxevent	1	/* # events between calculating individual event probabilities	
					   accelerates the simulation -- with loss of accuracy		*/
#define		tmaxstat	0.5	/* Time between calculating correlation functions		*/


#define		birth_move	1	/* 0: no movement at birth			1: move		*/
#define		bmove_kernel	2	/* 0: birth movement kernel flat
					   1: birth movement kernel tent		
					   2: birth movement kernel Gaussian				*/
#define		bnbr_kernel	2	/* 0: kernel for neighbour effects on birth flat
					   2: kernel for neighbour effects on birth Gaussian		*/
#define		dnbr_kernel	2	/* 0: kernel for neighbour effects on death flat
					   2: kernel for neighbour effects on death Gaussian		*/
#define		move_kernel	1	/* 0: adult movement kernel flat
					   1: adult movement kernel tent		
					   2: adult movement kernel Gaussian				*/


#define		rcfsmax		0.25	/* Maximum radius for computing correlation functions		*/
#define		binmax		10	/* Number of bins for rcfs and mom3				*/
					/* Number for scfs is 2*binmax for x,y				*/


FILE		*in, 
		*out1,			/* file for spatial pattern					*/
		*out2,			/* file for total population sizes				*/
		*out3,			/* file for radial correlation functions			*/
		*out4,			/* file for details in program					*/
		*out5,			/* file for spatial correlation functions			*/
		*out6,			/* file for mean path of population sizes			*/
		*out7,			/* file for mean path of rcfs					*/
		*out8;			/* file for third moment functions				*/


int		state_flag  = 0;	/* 1: detailed output for system state on		0: off	*/
int		nbr_flag    = 0;	/* 1: detailed output for neighbourhood calculations on	0: off	*/
int		choose_flag = 0;	/* 1: detailed output for choose event  calculations on	0: off	*/
int		birth_flag  = 0;	/* 1: detailed output for birth    calculations on	0: off	*/
int		death_flag  = 0;	/* 1: detailed output for death    calculations on	0: off	*/
int		move_flag   = 0;	/* 1: detailed output for movement calculations on	0: off	*/
int		scfs_flag   = 0;	/* 1: detailed output for scfs     calculations on	0: off	*/
int		rcfs_flag   = 0;	/* 1: detailed output for rcfs     calculations on	0: off	*/
int		mom3_flag   = 0;	/* 1: detailed output for moment 3 calculations on	0: off	*/
int		vogl_flag   = 0;	/* 1:          output for vogl graphics on		0: off	*/


/****************************************************************************************************************
* Global types													*
****************************************************************************************************************/

/***This structure holds the basic information on an individual***/
/***It is designed to create a bidirectional list for a species***/
struct	rep
	{	int	sp;				/* species 	*/
		float	x;				/* x-coordinate	*/
		float	y;				/* y-coordinate */
		float	b;				/* probability of birth    per unit time */
		float	d;				/* probability of death    per unit time */
		float	m;				/* probability of movement per unit time */
		int	bnbr;				/* Number of neighbours in birth nbrhood */
		int	dnbr;				/* Number of neighbours in death nbrhood */
		float	bwt;				/* Weight of neighbours in death nbrhood */
		float	dwt;				/* Weight of neighbours in death nbrhood */
		struct	rep	*nextP, *prevP;		/* makes structure a bidirectional list  */
	};

typedef	struct	rep	REP;	/*  REP  is a new data type for declaring structure of type rep */
typedef		REP	*repP;	/* *repP is pointer to the data type REP.  (C book p613)	*/


/***This structure holds information on properties of a species***/
struct	species_properties
	{	int	n;		/* Number of individuals					  */
		float	C1;		/* Density of individuals (info passed across from moment dynam)  */
		float	B0;		/* Intrinsic probability of birth per unit time			  */
		float	D0;		/* Intrinsic probability of death per unit time			  */
		float	B1[nspp];	/* Parameters for effect of neighbours on births		  */
		float	D1[nspp];	/* Parameters for effect of neighbours on births		  */
		float	Wbrmax[nspp];	/* Maximum radius over which neighbours affect births		  */
		float	Wbsd[nspp];	/* Standard deviation for neighbour eff on birth (Gaussian kernel)*/
		float	Wdrmax[nspp];	/* Maximum radius over which neighbours affect deaths		  */
		float	Wdsd[nspp];	/* Standard deviation for neighbour eff on death (Gaussian kernel)*/
		float	Mbrmax;		/* Maximum radius over which newborn individuals are displaced	  */
		float	Mbsd;		/* Standard deviation for birth movements (Gaussian kernel)	  */
		float	Mrmax;		/* Maximum radius of movement					  */
		float	Msd;		/* Standard deviation for adult movements (Gaussian kernel)	  */
		float	Mintegral;	/* Total probability of movement per unit time			  */
		float	BD;		/* Difference between total probability of birth and death	  */
		repP	firstP;		/* Pointer to first individual					  */
		repP	lastP;		/* Pointer to last  individual					  */
	};


/****************************************************************************************************************
* Global variables												*
****************************************************************************************************************/

/***This variable holds the current state of the system for all species***/
struct	rep	state[nspp];


/***This variable holds the properties of the species***/
struct	species_properties	param[nspp];


/***Maximum neighbourhood size -- used for dealing with interactions across the boundary***/
float	nbrmax;


/***Area of arena***/
float	arena_area;


/***Size of bin increments***/
float	bininc = rcfsmax/binmax;


/***Counter for stochastic realization number***/
int	iterate;


/***Time***/
float	t;		/* time	*/


/***Spatial correlation functions***/
/*The origin is between bin binmax-1 nd binmax*/
/* Column 0:	x coordinate
   Column 1:	y coordinate
   Column 2..:	spatial correlation	*/

float	scfs[2*binmax][2*binmax][1 + (nspp-1)*(nspp-1)];


/***Radial correlation functions***/
/* Column 0:	radial distance
   Column 1..:	radial correlation	*/

float	rcfs[binmax][1 + (nspp-1)*(nspp-1)];


/***Third momentfunctions***/

float	mom3[binmax][binmax][1 + (nspp-1)*(nspp-1)*(nspp-1)];


/***Variables for random number generator***/
double	drand48();	/* Function for random number generator */
int	seed;		/* Seed for random number generator	*/


/****************************************************************************************************************
* Procedures for vogl declared externally									*
****************************************************************************************************************/

#include	"simulate.vogl.c"


/****************************************************************************************************************
*  Function for generating random numbers on a normal distribution on [0,1]					*
*  From a method suggested by David.										*
****************************************************************************************************************/
double random_normal()
{
	int i;
	double theta, y, value1, value2;

	y     = drand48();
	y     = -2.0 * log(y);
	theta = drand48();
	theta = theta * 2.0 * pi;

	value1 = sqrt(y) * cos(theta);	/*only this one being used*/
	value2 = sqrt(y) * sin(theta);

	return(value1);
}


/****************************************************************************************************************
* Ulf's procedure for reporting errors when reading in parameters						*
****************************************************************************************************************/

void	report(format)
char	format[];
{
	printf("\nError while reading parameter file: %s not found.\n\n", format);
	exit(1);
}


/****************************************************************************************************************
* Ulf's procedure for inputing parameters (so that stochastic simulation and moment dynamics use the same	*
* input file)													*
****************************************************************************************************************/

void initialize_parameters()
{
	int	sp, i, bin, nt, iterate;
	char	name [40], format[40];

    /***Calculate area of arena***/
	arena_area = (xmax - xmin) * ( ymax - ymin);

    /***Loop for inputing parameters one species at a time***/
	for (sp=1; sp<nspp; ++sp)
	{
	    /***Read C1 (this is the starting density used by the moment dynamics program)***/
		sprintf(name, "C1[%d]", sp);			sprintf(format, "%s %%g\n", name);
		if (fscanf(in, format, &param[sp].C1)!=1)	report(name);

	    /***Convert starting density to a number (this program works on numbers)***/
		param[sp].n = (int) param[sp].C1 * arena_area;

	    /***Read B0's***/
		sprintf(name, "B0[%d]", sp);			sprintf(format, "%s %%g\n", name);
		if (fscanf(in, format, &param[sp].B0)!=1)	report(name);

	    /***Read D0's***/
		sprintf(name, "D0[%d]", sp);			sprintf(format, "%s %%g\n", name);    
		if (fscanf(in, format, &param[sp].D0)!=1)	report(name);    

	    /***Read B1's***/
		for (i=1;  i<nspp;  ++i)
		{	sprintf(name, "B1[%d][%d]", sp, i);		sprintf(format, "%s %%g\n", name);      
			if (fscanf(in, format, &param[sp].B1[i])!=1)	report(name);    
		}

	    /***Read D1's***/
		for (i=1;  i<nspp;  ++i)
		{	sprintf(name, "D1[%d][%d]", sp, i);		sprintf(format, "%s %%g\n", name);          
			if(fscanf(in, format, &param[sp].D1[i])!=1)	report(name);
		}

	    /***Read Wbrmax's***/
		for (i=1;  i<nspp;  ++i)
		{	sprintf(name,"Wbrmax[%d][%d]", sp, i);		sprintf(format, "%s %%g\n", name);
			if (fscanf(in, format,&param[sp].Wbrmax[i])!=1)	report(name);    
		}

	    /***Read Wbsd's***/
		for (i=1;  i<nspp;  ++i)
		{	sprintf(name,"Wbsd[%d][%d]", sp, i);		sprintf(format, "%s %%g\n", name);
			if (fscanf(in, format,&param[sp].Wbsd[i])!=1)	report(name);    
		}

	    /***Read Wdrmax's***/
		for (i=1;  i<nspp;  ++i)
		{	sprintf(name, "Wdrmax[%d][%d]", sp, i);		sprintf(format, "%s %%g\n", name);    
			if (fscanf(in, format,&param[sp].Wdrmax[i])!=1)	report(name);    
		}

	    /***Read Wdsd's***/
		for (i=1;  i<nspp;  ++i)
		{	sprintf(name,"Wdsd[%d][%d]", sp, i);		sprintf(format, "%s %%g\n", name);
			if (fscanf(in, format,&param[sp].Wdsd[i])!=1)	report(name);    
		}

	    /***Read Mbrmax's***/
		sprintf(name, "Mbrmax[%d]", sp);		sprintf(format, "%s %%g\n", name);              
		if (fscanf(in, format, &param[sp].Mbrmax)!=1)	report(name);

	    /***Read Mbsd's***/
		sprintf(name, "Mbsd[%d]", sp);			sprintf(format, "%s %%g\n", name);    
		if (fscanf(in, format, &param[sp].Mbsd)!=1)	report(name);    

	    /***Read Mrmax's***/
		sprintf(name, "Mrmax[%d]", sp);			sprintf(format, "%s %%g\n", name);              
		if (fscanf(in, format, &param[sp].Mrmax)!=1)	report(name);

	    /***Read Msd's***/
		sprintf(name, "Msd[%d]", sp);			sprintf(format, "%s %%g\n", name);    
		if (fscanf(in, format, &param[sp].Msd)!=1)	report(name);    

	    /***Read Mintegral's***/
		sprintf(name, "Mintegral[%d]", sp);		sprintf(format, "%s %%g\n", name);              
		if (fscanf(in, format,&param[sp].Mintegral)!=1)	report(name);
	}

    /***Calculate the size of the largest interaction neighbourhood (of Wbrmax and Wdrmax) across all species***/
    /***This is used for dealing with interactions across the boundary***/
	nbrmax = 0.0;
	for (sp=1; sp<nspp; ++sp)
	for (i=1;  i<nspp;  ++i)
	{	if (param[sp].Wbrmax[i] > nbrmax)	nbrmax = param[sp].Wbrmax[i];
		if (param[sp].Wdrmax[i] > nbrmax)	nbrmax = param[sp].Wdrmax[i];
	}
}


/****************************************************************************************************************
* This procedure outputs the parameters of the species								*
****************************************************************************************************************/

void	output_parameters()
{
	int	sp, i, bin;

    /***Record the parameters of species***/
	fprintf(out4, "Parameters of species\n");
	for (sp=1; sp<nspp; ++sp)
	{	fprintf(out4, "\nSPECIES:%d  Number of indviduals:%4d\n", sp, param[sp].n);
		fprintf(out4, "B0:%7.4f  Mbrmax:%7.4f  Mbsd:%7.4f\n", param[sp].B0, param[sp].Mbrmax, param[sp].Mbsd);
		fprintf(out4, "D0:%7.4f\n", param[sp].D0);
		for (i=1; i<nspp; ++i)		fprintf(out4, "B1[%d][%d]    :%7.4f    ", sp, i, param[sp].B1[i]);
		fprintf(out4, "\n");
		for (i=1; i<nspp; ++i)		fprintf(out4, "Wbrmax[%d][%d]:%7.4f    ", sp, i, param[sp].Wbrmax[i]);
		fprintf(out4, "\n");
		for (i=1; i<nspp; ++i)		fprintf(out4, "Wbsd[%d][%d]  :%7.4f    ", sp, i, param[sp].Wbsd[i]);
		fprintf(out4, "\n");
		for (i=1; i<nspp; ++i)		fprintf(out4, "D1[%d][%d]    :%7.4f    ", sp, i, param[sp].D1[i]);
		fprintf(out4, "\n");
		for (i=1; i<nspp; ++i)		fprintf(out4, "Wdrmax[%d][%d]:%7.4f    ", sp, i, param[sp].Wdrmax[i]);
		fprintf(out4, "\n");
		for (i=1; i<nspp; ++i)		fprintf(out4, "Wdsd[%d][%d]  :%7.4f    ", sp, i, param[sp].Wdsd[i]);
		fprintf(out4, "\n");
		fprintf(out4, "Mintegral:   %7.4f    Mrmax:%7.4f         Msd:%7.4f\n", param[sp].Mintegral, param[sp].Mrmax, param[sp].Msd);
	}

	fprintf(out4, "\nArena area                :%8.4f\n", arena_area);
	fprintf(out4, "\nMaximum interaction radius:%8.4f\n", nbrmax);
	fprintf(out4, "\nMaximum rcfs radius       :%8.4f\n", rcfsmax);
	fprintf(out4, "\n");
	if (bmove_kernel == 0)	fprintf(out4, "Birth movement kernel     :  flat\n");
	if (bmove_kernel == 1)	fprintf(out4, "Birth movement kernel     :  tent\n");
	if (bmove_kernel == 2)	fprintf(out4, "Birth movement kernel     :  Gaussian\n");
	if (bnbr_kernel == 0)	fprintf(out4, "Birth interation kernel   :  flat\n");
	if (bnbr_kernel == 2)	fprintf(out4, "Birth interation kernel   :  Gaussian\n");
	if (dnbr_kernel == 0)	fprintf(out4, "Death interation kernel   :  flat\n");
	if (dnbr_kernel == 2)	fprintf(out4, "Death interation kernel   :  Gaussian\n");

    /***Warn the user if the neighbourhood diameter is larger than the x-length of arena***/
	if (2.0*nbrmax > xmax-xmin)
	fprintf(out4, "\nWARNING:  radius of interaction neighbourhood exceeds arena length / 2\n");
}


/****************************************************************************************************************
* This procedure outputs the current state of the system							*
****************************************************************************************************************/

void	output_system_state()
{
	int	sp, i;
	repP	nowP;

    /***Output bidirectional list species by species***/
	fprintf(out4, "\n");
	for (sp=1; sp<nspp; ++sp)
	{	nowP = param[sp].firstP;
		fprintf(out4, "Species:%d  firstP:%d  lastP:%d\n", sp, param[sp].firstP, param[sp].lastP);
		fprintf(out4, "Species Individual   prevP     nowP    nextP  x-coord  y-coord    birth    death movement");
		fprintf(out4, " nbrb nbrd   nbrbwt   nbrdwt\n");
		i = 0;
		while (nowP != NULL)
		{	fprintf(out4, "%8d %8d %8d %8d %8d %8.4f %8.4f %8.4f %8.4f %8.4f %4d %4d %8.4f %8.4f\n", 
				nowP->sp, i, nowP->prevP, nowP, nowP->nextP, nowP->x, nowP->y, nowP->b, nowP->d, nowP->m, 
				nowP->bnbr, nowP->dnbr, nowP->bwt, nowP->dwt);
			i = i + 1;
			nowP = nowP->nextP;
		}
	}
/*	fflush(out4);
*/
}


/****************************************************************************************************************
* Procedure for placing spatial pattern in output file	for gnuplot						*
****************************************************************************************************************/
void	output_arena()
{
	char	filename[40];
	repP	nowP;
	int	sp, i;


    /***Output bidirectional list species by species***/
	fprintf(out4, "\n");
	for (sp=1; sp<nspp; ++sp)
	{	nowP = param[sp].firstP;
		i = 0;

	    /***Open a file for output of current species***/
		sprintf(filename, "output.%d.area.dat", sp);
		out1 = fopen(filename, "w");

		while (nowP != NULL)
		{	fprintf(out1, "%8d %8d %8.4f %8.4f %8.4f %8.4f %8.4f %4d %4d %8.4f %8.4f\n", 
			nowP->sp, i, nowP->x, nowP->y, nowP->b, nowP->d, nowP->m, nowP->bnbr, nowP->dnbr, nowP->bwt, nowP->dwt);
			i = i + 1;
			nowP = nowP->nextP;
		}
		fflush(out1);
		fclose(out1);
	}
}


/****************************************************************************************************************
* This procedure initializes the location of individuals in the arena						*
****************************************************************************************************************/

void	initialize_arena()
{
	int	sp, i;
	repP	firstP, lastP, nowP, prevP;

    /***Initialize individuals for the arena***/
	for (sp=1; sp<nspp; ++sp)
	{
	    /***Set the first and last pointers to null***/
		param[sp].firstP = NULL;
		param[sp].lastP  = NULL;

	    /***Set the start pointer as long as individuals actually exist***/
		if (param[sp].n != 0)
		{	param[sp].firstP = &state[sp];
			prevP = NULL;
			nowP  = param[sp].firstP;
			
		    /***Set up bidirectional list for all individuals by means of pointers***/
			for (i=1;   i<=param[sp].n;  ++i)
			{
			    /***Allcate the species of the individual***/
				nowP->sp = sp;
			
			    /***Allocate spatial locations of individuals uniformly at random (Poisson distribution)***/
				nowP->x = xmin + (xmax - xmin) * drand48();
				nowP->y = ymin + (ymax - ymin) * drand48();

			    /***Create the pointer links to previous and next individual***/
				nowP->prevP = prevP;
				nowP->nextP = NULL;
				if (i != param[sp].n)
				{	nowP->nextP = (repP) malloc(sizeof(REP));
					if (!nowP->nextP)
					{	fprintf(out4, "Unable to allocate memory to pointer: exiting from program\n");
						exit(1);
					}
				}

			    /***Update the pointers***/
				prevP = nowP;
				nowP  = nowP->nextP;
			}
			param[sp].lastP = prevP;
		}
	}

    /***Output initial state of system***/
	if (state_flag == 1)	output_system_state();

}


/****************************************************************************************************************
* Procedure for calculating distance between target and neighbour						*
* The complication here is that one has to allow for the periodic boundaries					*
****************************************************************************************************************/

void	neighbour_distance(targP, nbrP, radius, x, y, distance)
repP	targP;		/* pointer to the target individual		       */
repP	nbrP;		/* pointer to a neighbour of target individual (targP) */
float	radius;		/* radius of neighbourhood			       */
float	*x;		/* x coordinate of neighbour relative to target	       */
float	*y;		/* y coordinate of neighbour relative to target	       */
float	*distance;	/* distance from neighbour to target		       */
{
	int	i, j;
	int	N, E, S, W;
	float	X[3], Y[3], sqdist[5];	/* vector containing 4 squared distances */

    /***Start by setting distance and elements of sqdist to some large value***/
	for (i=1; i<3; ++i)	X[i]      = (xmax-xmin) * (ymax-ymin) * 1000.0;
	for (i=1; i<3; ++i)	Y[i]      = (xmax-xmin) * (ymax-ymin) * 1000.0;
	for (i=1; i<5; ++i)	sqdist[i] = (xmax-xmin) * (ymax-ymin) * 1000.0;
				*distance = (xmax-xmin) * (ymax-ymin) * 1000.0;

    /***Distance measurement (1): measured within the arena***/
    /***Location of neighbour relative to target***/
	X[1] = nbrP->x - targP->x;
	Y[1] = nbrP->y - targP->y;
	sqdist[1] = X[1]*X[1] + Y[1]*Y[1];
	*distance = sqdist[1];

    /***Establish which edges the target individual lies close to, and whether the neighbour lies clase to the opposite edge***/
	N = 0;	E = 0;	S = 0;	W = 0;
	if (targP->y + radius > ymax && nbrP->y - radius < ymin)	N = 1;	/* Target on north, neighbour on south */
	if (targP->x + radius > xmax && nbrP->x - radius < xmin)	E = 1;	/* Target on east,  neighbour on west  */
	if (targP->y - radius < ymin && nbrP->y + radius > ymax)	S = 1;	/* Target on south, neighbour on north */
	if (targP->x - radius < xmin && nbrP->x + radius > xmax)	W = 1;	/* Target on west,  neighbour on east  */

    /***If target and neighbour lie close to opposite boundaries 'move' the neighbour across***/
	if (N)		Y[2] =   (nbrP->y + (ymax-ymin))  -  targP->y;
	if (E)		X[2] =   (nbrP->x + (xmax-xmin))  -  targP->x;
	if (S)		Y[2] =   (nbrP->y - (ymax-ymin))  -  targP->y;
	if (W)		X[2] =   (nbrP->x - (xmax-xmin))  -  targP->x;
	sqdist[2] = X[2]*X[2] + Y[1]*Y[1];
	sqdist[3] = X[1]*X[1] + Y[2]*Y[2];
	sqdist[4] = X[2]*X[2] + Y[2]*Y[2];

    /***Take the square distance as the shortest of the four alternative placements of the neighbour***/
    /***as long as at least one of the boundary conditions is satisfied***/
	j = 1;
	if (N || E || S || W)	for (i=1; i<5; ++i)	if (sqdist[i] < *distance)
							{	*distance = sqdist[i];
								j = i;
							}

    /***Record the 'moved' coordinates of the neighbour relative to the target individual***/
	if (j==1)	{	*x = X[1];	*y = Y[1];	}
	if (j==2)	{	*x = X[2];	*y = Y[1];	}
	if (j==3)	{	*x = X[1];	*y = Y[2];	}
	if (j==4)	{	*x = X[2];	*y = Y[2];	}

    /***Record the distance***/
	if (*distance > -error)	*distance = sqrt(*distance);
	else
	{	fprintf(out4, "Distance (%7.4f) between target:%d and neighbour:%d is negative\n", *distance, targP, nbrP);
		fprintf(out4, "Exiting from program\n");
		exit(1);
	}
}


/****************************************************************************************************************
* Procedure for calculating weights on a target individual for births and deaths  due to neighbours		*
****************************************************************************************************************/

void	neighbourhood_weights(targP)
repP	targP;		/* Pointer to the target individual */
{
	int	sp, i;
	repP	nbrP;		/* pointer to a neighbour of target individual (targP)	*/
	float	x, y;		/* x and y coordinates of neighbour relative to target	*/
	float	distance;	/* distance from neighbour to target			*/
	float	rmax;		/* distance over which interaction kernel extends	*/
	float	var;		/* variance of Gaussian interaction kernel		*/
	float	norm;		/* Normalization of Gaussian kernel (truncated at rmax)	*/
	float	weight;		/* weight given to neighbour at distance r		*/

	targP->bnbr = 0;	targP->dnbr = 0;
	targP->bwt  = 0.0;	targP->dwt  = 0.0;

    /***Loop for calculating total birth and death weightings on a target individual by counting over all of its neighbours***/
	for (sp=1; sp<nspp; ++sp)
	{	nbrP = param[sp].firstP;
		while (nbrP != NULL)
		{
		    /***Target individual itself is not counted as a neighbour!***/
			if (nbrP != targP)
			{
			    /***Calculate the distance from target individual to neighbour***/
				neighbour_distance(targP, nbrP, nbrmax, &x, &y, &distance);
			
			    /***If the neighbour lies in the birth neighbourhood of the target individual, add it on to the birth weight***/
			    /***This is done only if rmax > 0 and if there is density-dependent birth***/
				rmax = param[targP->sp].Wbrmax[nbrP->sp];
				if (distance < rmax && rmax > error && param[targP->sp].B1[nbrP->sp] > error)
				{	targP->bnbr = targP->bnbr + 1;
					if (bnbr_kernel == 0)	weight = 1.0 / (pi * rmax * rmax);					/*Flat kernel*/
					if (bnbr_kernel == 2)
					{			var    = param[targP->sp].Wbsd[nbrP->sp] * param[targP->sp].Wbsd[nbrP->sp];
								norm   = 2.0 * pi * var * (1.0 - exp(-rmax*rmax / (2.0*var)));		/*Gaussian kernel*/
								weight = exp(-distance*distance / (2.0*var)) / norm;
					}
					targP->bwt  = targP->bwt  + param[targP->sp].B1[nbrP->sp] * weight;
				}
				
			    /***If the neighbour lies in the death neighbourhood of the target individual, add it on to the death weight***/
			    /***This is done only if rmax > 0 and if there is density-dependent death***/
				rmax = param[targP->sp].Wdrmax[nbrP->sp];
				if (distance < rmax && rmax > error && param[targP->sp].D1[nbrP->sp] > error)
				{	targP->dnbr = targP->dnbr + 1;
					if (dnbr_kernel == 0)	weight = 1.0 / (pi * rmax * rmax);					/*Flat kernel*/
					if (dnbr_kernel == 2)
					{			var    = param[targP->sp].Wdsd[nbrP->sp] * param[targP->sp].Wdsd[nbrP->sp];
								norm   = 2.0 * pi * var * (1.0 - exp(-rmax*rmax / (2.0*var)));		/*Gaussian kernel*/
								weight = exp(-distance*distance / (2.0*var)) / norm;
					}
					targP->dwt  = targP->dwt  + param[targP->sp].D1[nbrP->sp] * weight;
				}

			    /***Record details of calculations on neighbours***/
				if (nbr_flag == 1)
				{	fprintf(out4, "\ttargP:%d  nbrP:%d   species:%d  x:%7.4f  y:%7.4f  distance:%6.4f", 
						targP, nbrP, nbrP->sp, x, y, distance);
					if (distance < nbrmax)	fprintf(out4, "*");	else	fprintf(out4, " ");
					fprintf(out4, "  bnbr:%d  dnbr:%d  wtb:%7.4f  wtd:%7.4f \n", nbrP->bnbr, nbrP->dnbr, nbrP->bwt, nbrP->dwt);
					fflush(out4);
				}
			}
			nbrP = nbrP->nextP;
		}
	}
}


/****************************************************************************************************************
* Procedure for choosing probabilities of birth death and movement per unit time of each target individual	*
****************************************************************************************************************/

void	individual_event_probabilities()
{
	int	sp, i;
	repP	targP;		/* pointer to the target individual whose neighbourhood is to be evaluated */

    /***Births and deaths are allowed to be dependent on neighbourhood***/
    /***But the routines for doing this use a lot of time -- only enter them if B0 or D0 > 0***/
	for (sp=1; sp<nspp; ++sp)
	if (param[sp].B0 > error || param[sp].D0 > error)
	{	targP = param[sp].firstP;
		while (targP != NULL)
		{	
		    /***Record details of calculations on neighbours***/
			if (nbr_flag == 1)
			fprintf(out4, "\nTarget individual: species:%d  targP:%d\n", targP->sp, targP);
			fflush(out4);

		    /***Go to the neighbourhood of each individual and add up the weights***/
			neighbourhood_weights(targP);

		    /***Calculate the probability per unit time of the target individual's birth and death***/
			targP->b = param[sp].B0 + targP->bwt;
			targP->d = param[sp].D0 + targP->dwt;

		    /***Jump out of prgram if a negative birth or death rate has been encountered***/
			if (targP->b < 0.0)
			{	printf("Negative birth rate encountered:  program halted\n");
				exit(1);
			}
			if (targP->d < 0.0)
			{	printf("Negative death rate encountered:  program halted\n");
				exit(1);
			}

			targP = targP->nextP;
		}
	}

    /***Movements are assumed to be independent of neighbourhood***/
	for (sp=1; sp<nspp; ++sp)
	{	targP = param[sp].firstP;
		while (targP != NULL)
		{	targP->m = param[sp].Mintegral;
			targP    = targP->nextP;
		}
	}
}


/****************************************************************************************************************
* Procedure for choosing births, deaths and movements at random							*
****************************************************************************************************************/

void	choose_event(dt, event, pickP)
float	*dt;		/* time increment to the event				*/
char	*event;		/* event chosen in each iteration: b, d, m		*/
repP	*pickP;		/* pointer to the individual chosen for event		*/
{
	int	sp, i;
	repP	nowP;
	float	Bsum[nspp], Dsum[nspp];
	float	B, D, M, sum, rn, add;

    /***Probability per unit time of birth event***/
	for (sp=1; sp<nspp; ++sp)
	{	Bsum[sp] = 0.0;
		nowP  = param[sp].firstP;
		while (nowP != NULL)
		{	Bsum[sp] = Bsum[sp] + nowP->b;
			nowP  = nowP->nextP;
		}
	}

    /***Probability per unit time of death event***/
	for (sp=1; sp<nspp; ++sp)
	{	Dsum[sp] = 0.0;
		nowP  = param[sp].firstP;
		while (nowP != NULL)
		{	Dsum[sp] = Dsum[sp] + nowP->d;
			nowP  = nowP->nextP;
		}
	}

    /***Probability per unit time of movement event***/
	M = 0.0;
	for (sp=1; sp<nspp; ++sp)
	{	nowP = param[sp].firstP;
		while (nowP != NULL)
		{	M = M + nowP->m;
			nowP = nowP->nextP;
		}
	}

    /***Keep a record of the difference between total birth and death probability p.u.t***/
	for (sp=1; sp<nspp; ++sp)
	param[sp].BD = Bsum[sp] - Dsum[sp];

    /***Sum of probabilities of all events***/
	B = 0.0;
	D = 0.0;
	for (sp=1; sp<nspp; ++sp)
	{	B = B + Bsum[sp];
		D = D + Dsum[sp];
	}
	sum = B + D + M;

    /***Random time to next event (exponential distribution of waiting times)***/
	if (sum > error)	*dt = -log(1.0 - drand48())/sum;
/*	else
	{	*dt    = tmax - t;
		*event = 'o';
		*pickP = NULL;
printf("Test dt:%6.3f event:%c\n", *dt, *event);
		goto     label_end;
	}
*/
    /***Jump out of program if sum of event probabilities = 0***/
	else	
	{	printf("Sum of event probabilities = 0:  program halted\n");
		exit(1);
	}

    /***Choose at random a birth, death or encounter event***/
	rn = drand48();
	if (rn <  B/sum)			*event = 'b';
	if (rn >= B/sum && rn < (B+D)/sum)	*event = 'd';
	if (rn >= (B+D)/sum)			*event = 'm';

    /***Initialize variables for working out the individual to which the event happens***/
	rn   = drand48();
	add  = 0.0;

    /***Choose at random the individual to which event happens***/
	switch (*event)
	{
		case 'b': /*birth event*/
			for (sp=1; sp<nspp; ++sp)
			{	nowP = param[sp].firstP;
				while (nowP != NULL)
				{	add = add + nowP->b/B;
					if (rn < add)
					{	add = add - nowP->b/B;
						goto label_b;
					}
					nowP  = nowP->nextP;
				}
			}
		label_b:break;

		case 'd': /*death event*/
			for (sp=1; sp<nspp; ++sp)
			{	nowP = param[sp].firstP;
				while (nowP != NULL)
				{	add = add + nowP->d/D;
					if (rn < add)
					{	add = add - nowP->d/D;
						goto label_d;
					}
					nowP  = nowP->nextP;
				}
			}
		label_d:break;

		case 'm': /*movement event*/
			for (sp=1; sp<nspp; ++sp)
			{	nowP = param[sp].firstP;
				while (nowP != NULL)
				{	add = add + nowP->m/M;
					if (rn < add)
					{	add = add - nowP->m/M;
						goto label_m;
					}
					nowP  = nowP->nextP;
				}
			}
		label_m:break;
	}
	*pickP   = nowP;

    /***Record of calculations***/
	label_end:
	if (choose_flag == 1)
	{	fprintf(out4, "\nTime:%5.3f  State at end of choose_event procedure\n", t+*dt);
		fprintf(out4, "B:%6.4f  D:%6.4f  M:%6.4f  ", B, D, M);
		fprintf(out4, "param[%d].firstP:%d rn:%6.4f  add:%6.4f  ", sp, param[nowP->sp].firstP, rn, add);
		fprintf(out4, "species:%d  event:%c  pickP:%d\n", nowP->sp, *event, *pickP);
		fflush(out4);
	}
}


/****************************************************************************************************************
* Birth of individual												*
****************************************************************************************************************/

void	birth(pickP)
repP	pickP;		/* pointer to the individual chosen for event */
{
	double	x, y, r, a;
	repP	birthP;	/* pointer to the newborn individual	      */

    /***Allocate memory for newborn individual***/
	birthP = NULL;
	birthP = (repP) malloc(sizeof(REP));
	if (!birthP)
	{	fprintf(out4, "Unable to allocate memory to pointer: exiting from program\n");
		exit(1);
	}

    /***WARNING!  It is essential to alter the pointers in the correct sequence below***/

    /***If the new individual is at the end of the list, the pointer to the next individual is set to NULL***/
    /*** and a record of the tail pointer is made***/
	if (pickP->nextP == NULL)
	{	birthP->nextP = NULL;
		param[pickP->sp].lastP = birthP;
	}

    /***If the new individual is not at the end of the list, it is linked to the next individual ***/
	else
	{	birthP->nextP = pickP->nextP;
		pickP->nextP->prevP = birthP;
	}

   /***Lastly the links are made between the parent and the new individual***/
	birthP->prevP = pickP;
	pickP->nextP  = birthP;

    /***The new individual is obviously of the same species as the parent***/
	birthP->sp = pickP->sp;

    /***Model newborn individuals born at the location of parent***/
	if (!birth_move)
	{
	    /***Here the new individual is assumed to be born at the location of the parent***/
		birthP->x  = pickP->x;
		birthP->y  = pickP->y;
	}

    /***Model newborn individuals born displaced from the location of parent***/
	if (birth_move)
	{
	    /***Flat movement kernel: Distances x and y taken from uniform distibutions***/
	    /***		      Upper limit set by radius Mbrmax***/
		if (bmove_kernel == 0)
		{	do
			{	x = (2.0 * drand48() - 1.0) * param[pickP->sp].Mbrmax;
				y = (2.0 * drand48() - 1.0) * param[pickP->sp].Mbrmax;
			}	while	(x*x + y*y > param[pickP->sp].Mbrmax * param[pickP->sp].Mbrmax);
			birthP->x = pickP->x + x;
			birthP->y = pickP->y + y;
		}

	    /***Tent movement kernel: Distance moved taken from uniform distribution up to a radius Mbrmax***/
	    /***                      Angle of movement taken from a uniform distribution***/
		if (bmove_kernel == 1)
		{	r = drand48() * param[pickP->sp].Mbrmax;
			a = drand48() * 2.0 * pi;
			birthP->x = pickP->x + r * cos(a);
			birthP->y = pickP->y + r * sin(a);
		}

	    /***Gaussian movement kernel: Distance moved taken from Gaussian distribution up to a radius Mbrmax***/
	    /***                          Angle of movement taken from a uniform distribution***/
		if (bmove_kernel == 2)
		{	do
			{	r = random_normal() * param[pickP->sp].Mbsd;
			}	while (error > r || r > param[pickP->sp].Mbrmax);
			a = drand48() * 2.0 * pi;
			birthP->x = pickP->x + r * cos(a);
			birthP->y = pickP->y + r * sin(a);
		}

		if (bmove_kernel!=0 && bmove_kernel!=1 && bmove_kernel!=2)
		{	printf("Error in choice of movement kernel: Program terminating\n");	exit(1);
		}

	    /***Wrap-around edges***/
		if (birthP->x < xmin)	birthP->x = birthP->x + (xmax - xmin);
		if (birthP->x > xmax)	birthP->x = birthP->x - (xmax - xmin);
		if (birthP->y < ymin)	birthP->y = birthP->y + (ymax - ymin);
		if (birthP->y > ymax)	birthP->y = birthP->y - (ymax - ymin);
	}

    /***The other parameters are initialised at zero***/
	birthP->b    = 0.0;
	birthP->d    = 0.0;
	birthP->m    = 0.0;
	birthP->bnbr = 0;
	birthP->dnbr = 0;
	birthP->bwt  = 0.0;
	birthP->dwt  = 0.0;

    /***The number of individuals of this species is increased by 1***/
	param[pickP->sp].n = param[pickP->sp].n + 1;

    /***Record of calculations***/
	if (birth_flag == 1)
	{	fprintf(out4, "\nTime:%5.3f  State at end of birth procedure\n", t);
		fprintf(out4, " pickP->sp:%d   pickP->x:%6.4f   pickP->y:%6.4f   pickP->prevP:%d   pickP:%d   pickP->nextP:%d\n", 
		pickP->sp, pickP->x, pickP->y, pickP->prevP, pickP, pickP->nextP);
		fprintf(out4, "birthP->sp:%d  birthP->x:%6.4f  birthP->y:%6.4f  birthP->prevP:%d  birthP:%d  birthP->nextP:%d\n",
		birthP->sp, birthP->x, birthP->y, birthP->prevP, birthP, birthP->nextP);
	}

    /***Terminate program if newborn individual is outside the arena***/
	if (birthP->x < xmin || birthP->x > xmax || birthP->y < ymin || birthP->y > ymax)
	{	printf("birthP:%d  x:%8.4f  y:%8.4f\n", birthP, birthP->x, birthP->y);
		printf("Program terminating\n");
		exit(1);
	}
}


/****************************************************************************************************************
* Death of individual												*
****************************************************************************************************************/

void	death(pickP)
repP	pickP;		/* pointer to the individual chosen for event */
{

    /***If the dead individual is in the interior of the list, the the two adjoining individauls are linked***/
	if (pickP->nextP !=NULL && pickP->prevP !=NULL)
	{	pickP->prevP->nextP = pickP->nextP;
		pickP->nextP->prevP = pickP->prevP;
	}

    /***If the dead individual is at the head of the list, the next individual becomes the head***/
	if (pickP->prevP == NULL && pickP->nextP != NULL)
	{	pickP->nextP->prevP = NULL;
		param[pickP->sp].firstP = pickP->nextP;
	}

    /***If the dead individual is at the tail of the list, the previous individual becomes the tail***/
	if (pickP->nextP == NULL && pickP->prevP != NULL)
	{	pickP->prevP->nextP = NULL;
		param[pickP->sp].lastP = pickP->prevP;
	}

    /***If the dead individual is at the head AND the tail of the list, the species is extinct***/
	if (pickP->nextP == NULL && pickP->prevP == NULL)
	{	param[pickP->sp].firstP = NULL;
		param[pickP->sp].lastP  = NULL;
	}

    /***The number of individuals of this species is decreased by 1***/
	param[pickP->sp].n = param[pickP->sp].n - 1;

    /***Record of calculations***/
	if (death_flag == 1)
	{	fprintf(out4, "\nTime:%5.3f  State at end of death procedure\n", t);
		fprintf(out4, " pickP->sp:%d   pickP->x:%6.4f   pickP->y:%6.4f   pickP:%d\n",
		pickP->sp, pickP->x, pickP->y, pickP);
	}
}


/****************************************************************************************************************
* Move individual												*
****************************************************************************************************************/

void	movement(pickP)
repP	pickP;		/* pointer to the individual chosen for event */
{
	double	x, y, r, a;

    /***Flat movement kernel: Distances x and y taken from uniform distibutions***/
    /***		      Upper limit set by radius Mrmax***/
	if (move_kernel == 0)
	{	do
		{	x = (2.0 * drand48() - 1.0) * param[pickP->sp].Mrmax;
			y = (2.0 * drand48() - 1.0) * param[pickP->sp].Mrmax;
		}	while	(x*x + y*y > param[pickP->sp].Mrmax * param[pickP->sp].Mrmax);
		pickP->x = pickP->x + x;
		pickP->y = pickP->y + y;
	}

    /***Tent movement kernel: Distance moved taken from uniform distribution up to a radius Mrmax***/
    /***                      Angle of movement taken from a uniform distribution***/
	if (move_kernel == 1)
	{	r = drand48() * param[pickP->sp].Mrmax;
		a = drand48() * 2.0 * pi;
		pickP->x = pickP->x + r * cos(a);
		pickP->y = pickP->y + r * sin(a);
	}

    /***Gaussian movement kernel: Distance moved taken from Gaussian distribution up to a radius Mbrmax***/
    /***                          Angle of movement taken from a uniform distribution***/
	if (move_kernel == 2)
	{	do
		{	r = random_normal() * param[pickP->sp].Msd;
		}	while (error > r || r > param[pickP->sp].Mrmax);
		a = drand48() * 2.0 * pi;
		pickP->x = pickP->x + r * cos(a);
		pickP->y = pickP->y + r * sin(a);
	}

	if (move_kernel!=0 && move_kernel!=1 && move_kernel!=2)
	{	printf("Error in choice of movement kernel: Program terminating\n");	exit(1);
	}

    /***Wrap-around edges***/
	if (pickP->x < xmin)	pickP->x = pickP->x + (xmax - xmin);
	if (pickP->x > xmax)	pickP->x = pickP->x - (xmax - xmin);
	if (pickP->y < ymin)	pickP->y = pickP->y + (ymax - ymin);
	if (pickP->y > ymax)	pickP->y = pickP->y - (ymax - ymin);

    /***Record of calculations***/
	if (move_flag == 1)
	{	fprintf(out4, "\nTime:%5.3f  State at end of movement procedure\n", t);
		fprintf(out4, " pickP->sp:%d   pickP->x:%6.4f   pickP->y:%6.4f   pickP:%d  \n", 
		pickP->sp, pickP->x, pickP->y, pickP);
	}

    /***Terminate program if an individual has moved outside the arena***/
	if (pickP->x < xmin || pickP->x > xmax || pickP->y < ymin || pickP->y > ymax)
	{	printf("pickP:%d  x:%8.4f  y:%8.4f\n", pickP, pickP->x, pickP->y);
		printf("Program terminating\n");
		exit(1);
	}
}


/****************************************************************************************************************
* Procedures for calculating spatial statistics declared externally						*
****************************************************************************************************************/

#include	"simulate.correlations.c"


/****************************************************************************************************************
* Stochastic simulation												*
****************************************************************************************************************/

void	simulate()
{
	int	i, j, sp;

	float	dt;		/* time increment to the event				*/
	int	nt;		/* index for time in running average			*/
	float	tlower;		/* lower bound on time for calculating mean numbers	*/
	float	sum;		/* sum for calculating mean numbers at nt		*/
	char	event;		/* event chosen in each iteration: b, d, m		*/
	int	nevent;		/* counter for number of events	for updating event prbs	*/
	float	tstat;		/* counter for number of events	for correl functions	*/
	repP	pickP;		/* pointer to the individual chosen for event		*/

    /***Set up starting conditions for the iterations***/
	t      = 0.0;
	tlower = 0.0;
	sum    = 0.0;
	nt     = 0;
	nevent = 0;
	tstat  = 0.0;

    /***Initialize the arena***/
	initialize_arena();
	output_arena();

    /***Initialize the individual event probabilities***/
	individual_event_probabilities();

    /***Record the initial numbers in output file***/
	fprintf(out2, "%6.4f  ", t);
	for (sp=1; sp<nspp; ++sp)	fprintf(out2, "%3d  ", param[sp].n);
	fprintf(out2, "\n");

    /***Calculate initial spatial and radial correlation functions***/
	spatial_correlation_functions();
	radial_correlation_functions();

	t = 0.0;

    /***Begin the iterations***/
	while (t<=tmax)
	{
	    /***Output state of system***/
		if (state_flag == 1)	output_system_state();

	    /***Calculate the birth, death and movement probabilities per unit time of each individual***/
	    /***maxevent gives the number of iterations before individual birth, death, movement probabilities are updated***/
		nevent = nevent + 1;
		if (nevent == maxevent)
		{	individual_event_probabilities();

		    /***Graphical display of state of system***/
			if (vogl_flag == 1)	display_vogl();

			nevent = 0;
		}

	    /***Choose a birth, death or movement event at random***/
		choose_event(&dt, &event, &pickP);

	    /***Update the time***/
		t     = t     + dt;
		tstat = tstat + dt;

	    /***Update the state of the system***/
		switch	(event)
		{	case 'b':	birth(pickP);		break;
			case 'd':	death(pickP);		break;
			case 'm':	movement(pickP);	break;
			case 'o':	break;
		}

	    /***Check that numbers have not gone negative (this can happen if event probabilities not updated every iteration***/
		for (sp=1; sp<nspp; ++sp)	if (param[sp].n < 0)
		{	param[sp].n = 0;
			printf("Warning: species %d numbers have gone negative (returned to zero)\n", sp);
		}

	    /***Record summary results on screen***/
		printf("Realization:%2d  ", iterate);
		printf("t:%6.4f  sp:%d  x:%6.4f  y:%6.4f  %c  ", t, pickP->sp, pickP->x, pickP->y, event);
		for (sp=1; sp<nspp; ++sp)	printf("B-D:%6.3f  ", param[sp].BD);
		for (sp=1; sp<nspp; ++sp)	printf("%3d  ", param[sp].n);
/*		printf("  address:%d  ", pickP);
*/		printf("\n");

	    /***Record total population sizes in output file***/
		fprintf(out2, "%6.4f  ", t);
		for (sp=1; sp<nspp; ++sp)	fprintf(out2, "%3d  ", param[sp].n);
		fprintf(out2, "\n");
		fflush(out2);

	    /***Record of statistics made after every nmaxevent***/
		if (tstat >= tmaxstat)
		{
		    /***Put copy of spatial pattern into area file***/
			output_arena();

		    /***Calculate spatial and radial correlation functions***/
			spatial_correlation_functions();
			radial_correlation_functions();

			tstat = 0.0 +tstat - tmaxstat;
		}
	}

    /***Output of 2nd and third moments at end of simulation***/
/*	output_arena();
	radial_correlation_functions();
	third_moment_functions();
*/
    /***Pause at end of simulation***/
	if (vogl_flag == 1)	system("sleep 5");
}


/****************************************************************************************************************
* Main body of program												*
****************************************************************************************************************/

void	main()
{
	char	filename[20];
	int	sp, n1, n2, nt;

    /***Open files (channels out1, out2, out3 out6 and out7 are reserved)***/
	in   = fopen("input.dat","r");
	out4 = fopen("output.detail.dat", "w");
	out5 = fopen("output.scfs.dat", "w");

    /***Seed the random number generator***/
	printf("\nSeed for random number generator (4 digits):    ");
	scanf("%4d", &seed);
	srand48(seed);
	fprintf(out4, "Seed for random number generator:%d\n\n", seed);

    /***Initialize graphical display***/
	if (vogl_flag == 1)	initialize_vogl();

    /***Initialize the parameters of species and the arena***/
	initialize_parameters();
	output_parameters();

    /***Run a specified number of stochastic realizations from a single initial value***/
	if (nspp!=3 || grid!=1)	
	{
		for (iterate=0; iterate<niterate; ++iterate)
		{
		    /***Set starting numbers from starting density used in moment dynamics (this program works on numbers)***/
			for (sp=1; sp<nspp; ++sp)	param[sp].n = (int) param[sp].C1 * arena_area;

		    /***Open files for tpss, rcfs and mom3 output***/
			sprintf(filename, "output.%d.tpss.dat", iterate);
			out2 = fopen(filename, "w");
			sprintf(filename, "output.%d.rcfs.dat", iterate);
			out3 = fopen(filename, "w");
			sprintf(filename, "output.%d.mom3.dat", iterate);
			out8 = fopen(filename, "w");

		    /***Stochastic simulation***/
			simulate();
			
		    /***Close the output files***/
			fflush(out2);
			fflush(out3);
			fflush(out8);
			fclose(out2);
			fclose(out3);
			fclose(out8);
		}

	    /***Obtain mean path of realizations***/
		mean_path_tpss();
	}

    /***Run a grid of stochastic realizations for constructing a 2-species phase plane***/
	n1 = n1start;
	while (nspp==3 && grid==1 && n1 <= n1stop)
	{	n2 = n2start;
		while (n2 <= n2stop)
		{
			for (iterate=0; iterate<niterate; ++iterate)
			{
			    /***Set starting numbers***/
				param[1].n = n1;
				param[2].n = n2;

			    /***Open files for tpss and rcfs output***/
				sprintf(filename, "output.%d.%d.%d.tpss.dat", n1, n2, iterate);
				out2 = fopen(filename, "w");
				sprintf(filename, "output.%d.%d.%d.rcfs.dat", n1, n2, iterate);
				out3 = fopen(filename, "w");
				sprintf(filename, "output.%d.mom3.dat", iterate);
				out8 = fopen(filename, "w");

			    /***Stochastic simulation***/
				simulate();

			    /***Close the output files***/
				fflush(out2);
				fflush(out3);
				fflush(out8);
				fclose(out2);
				fclose(out3);
				fclose(out8);
			}

		    /***Obtain mean path of realizations***/
			mean_path_tpss();

			n2 = n2 + n2inc;
		}
		n1 = n1 + n1inc;
	}

	fflush(out4);
	fflush(out5);
	fflush(out6);
	fflush(out7);
	fclose(in);
	fclose(out4);
	fclose(out5);
	fclose(out6);
	fclose(out7);
}
