/****************************************************************************************************************
* This is a collection of procedures using vogl graphics to provide a display for stochastic simulation of	*
* for stochastic simulation of a multispecies community in continuous space and	time.				*
*														*
* Last update 11.2.98.												*
****************************************************************************************************************/

/***Global constants for vogl graphics***/

#define		color_code 1	/*0: white to black	1: black to red		*/
#define		window     120.0/*Size of window				*/
#define		txt_margin 0.3	/*Gives a margin around edge for text		*/
#define		sep_margin 0.1	/*Gives a margin between grid displays		*/
#define		col_margin 0.01	/*Gives a margin for drawing around the plot	*/
#define		nmax       25.0 /*Scaling of numbers to maximum colour number	*/



/****************************************************************************************************************
* This procedure carries out colour mixing.  Colouring depends on color_code:					*
*	0:	white to black	(for postscript)								*
*	1:	black to red	(for vogl)									*
* Procedure supplied by Ulf.											*
****************************************************************************************************************/

void fraction2rgb(fraction, rP, gP, bP)
double	fraction, *rP, *gP, *bP;
{
	int	huesector;
	double	hue, huetune, mix_up, mix_do, r, g, b;

	if(color_code==1)
	{
		hue=1.0-fraction;
		if(hue<0.0)	hue=0.0;
		if(hue>1.0)	hue=1.0;
		huesector=(int)floor(hue*5.0);
		huetune=hue*5.0-huesector;
		mix_up=huetune;
		mix_do=1.0-huetune;
		mix_up=pow(mix_up,1.0/2.5);
		mix_do=pow(mix_do,1.0/2.5);
		switch(huesector)
		{
			case 0 : r=1.0   ;g=mix_up;b=0.0   ;break; /* red    to yellow */
			case 1 : r=mix_do;g=1.0   ;b=0.0   ;break; /* yellow to green  */
			case 2 : r=0.0   ;g=1.0   ;b=mix_up;break; /* green  to cyan   */
			case 3 : r=0.0   ;g=mix_do;b=1.0   ;break; /* cyan   to blue   */
			case 4 : r=0.0   ;g=0.0   ;b=mix_do;break; /* blue   to black  */
			default: r=0.0   ;g=0.0   ;b=0.0   ;break;
		}
	}

	else
	{
		r=1.0-fraction;
		g=1.0-fraction;
		b=1.0-fraction;    
	}

	*rP=r;*gP=g;*bP=b;
}


/****************************************************************************************************************
* This procedure generates a color scale from black to red for coding numbers.					*
* Procedure supplied by Ulf.											*
****************************************************************************************************************/

void colors_init()
{
	int	index;
	double	r, g, b;

    /***Colours 0 to 7 are predefined***/

    /***colors 8 to 15 are less bright versions of predefined colours***/
	mapcolor( 8,  0,  0,  0);
	mapcolor( 9,100,  0,  0);
	mapcolor(10,  0,100,  0);
	mapcolor(11, 50, 50,  0);
	mapcolor(12,  0,  0,100);
	mapcolor(13, 50,  0, 50);
	mapcolor(14,  0, 50, 50);
	mapcolor(15,153,153,153);  

    /***colors 16 to 79 are mixed as below (64 colors)***/
	for(index=16; index<=79; index++)
	{	fraction2rgb((double)(index-16)/(double)(80-16),&r,&g,&b);
		mapcolor(index,(int)(r*255.0),(int)(g*255.0),(int)(b*255.0));
	}
}


/****************************************************************************************************************
* This procedure initializes a window for vogl graphics								*
****************************************************************************************************************/

void	initialize_vogl()
{
	double	width, height;

    /***Initialize a window***/
	width  = txt_margin + 4.0 * (1.0 + sep_margin);
	height = txt_margin + 4.0 * (1.0 + sep_margin) + txt_margin;
	prefposition(0, (int)(width*window), 0, (int)(height*window));
	vinit("X11");
	winopen("Stochastic simulation");	/*Window must have a name before display appears*/

    /***Set background colour***/
	color(BLACK);
	clear();

    /***Initialize colour scheme***/
	colors_init();

    /***Set up a coordinate frame***/
	ortho2(0.0, width, 0.0, height);
}


/****************************************************************************************************************
* This procedure displays the status of the simulation.								*
****************************************************************************************************************/

void	display_status()
{
	int	sp;
	double	xshift, yshift;
	double	x_lower, x_upper, y_lower, y_upper;
	double	x1, y1;
	char	labelling[50];

    /***Set the origin for output of status information***/
	xshift = txt_margin;
	yshift = txt_margin + 4.0*(1.0 + sep_margin) + txt_margin/2.0;

    /***Draw boundary of status information***/
	x_lower = xshift - col_margin;
	y_lower = yshift - col_margin;
	x_upper = xshift + col_margin + (nspp-1.0)/2.0;
	y_upper = yshift + col_margin + txt_margin;
	color(BLACK);
	polymode(PYM_LINE);
	rectf(x_lower, y_lower, x_upper, y_upper);

    /***Insert the text***/
	color(3);
	hfont("times.rb");
	hcentertext(0);
	htextsize(txt_margin/3.0, txt_margin/2.2);
	htextang(0.0);
	sprintf(labelling, "Community with %d species", nspp-1);
	x1  = xshift;
	y1  = yshift;
	move2(x1, y1);
	hcharstr(labelling);

	sprintf(labelling, "Numbers:");
	yshift = txt_margin + 4.0*(1.0 + sep_margin);
	x1  = xshift;
	y1  = yshift;
	move2(x1, y1);
	hcharstr(labelling);
	x1  = xshift + 0.3;

	for (sp=1; sp<nspp; ++sp)
	{	color(sp);
		sprintf(labelling, " %3d", param[sp].n);
		x1  = x1 + 0.2;
		move2(x1, y1);
		hcharstr(labelling);
	}

	color(3);
	sprintf(labelling, "Time:");
	yshift = txt_margin + 4.0*(1.0 + sep_margin) - txt_margin/2.0;
	x1  = xshift;
	y1  = yshift;
	move2(x1, y1);
	hcharstr(labelling);

	sprintf(labelling, "%8.2f", t);
	yshift = txt_margin + 4.0*(1.0 + sep_margin) - txt_margin/2.0;
	x1  = xshift + 0.3;
	y1  = yshift;
	move2(x1, y1);
	hcharstr(labelling);
}


/****************************************************************************************************************
* This procedure displays the simulation.								*
****************************************************************************************************************/

void	display_simulation()
{
	int	sp;
	repP	nowP;
	float	x, y, xshift, yshift;
	double	cell_color;
	char	labelling[50];

    /***Set the origin for labelling***/
	xshift = txt_margin;
	yshift = sep_margin;

    /***Insert the labelling***/
	color(3);
	hfont("times.rb");
	hcentertext(0);
	htextsize(txt_margin/3.0, txt_margin/2.2);
	htextang(0.0);
	sprintf(labelling, "Simulation");
	x  = xshift;
	y  = yshift;
	move2(x, y);
	hcharstr(labelling);

    /***Set the origin for the plot***/
	xshift = txt_margin + 0.0;
	yshift = txt_margin + 0.0;

    /***Draw the plot boundary***/
	color(7);
	polymode(PYM_LINE);
	rectf(xshift-col_margin, yshift-col_margin, xshift+2.0+col_margin, yshift+2.0+col_margin);


    /***Loop for each species***/
	for (sp=1; sp<nspp; ++sp)
	{
	    /***Plot the location of each individual of current species***/
		color(sp);
		polymode(PYM_FILL);
		nowP = param[sp].firstP;
		while (nowP != NULL)
		{	x = nowP->x;
			y = nowP->y;
			x = xshift + 2*x/(xmax-xmin);
			y = yshift + 2*y/(ymax-ymin);
			move2(x, y);
			circf(x, y, 0.015);
			nowP = nowP->nextP;
		}
	}
}


/****************************************************************************************************************
* This procedure displays the current rcfs of the simulation							*
****************************************************************************************************************/

void	display_rcfs()
{
	int	i, j, bin;
	double	x1, y1, xshift, yshift;
	double	x_lower, x_upper, y_lower, y_upper;
	double	rcf_block;
	double	yscale = 0.2;		/*Vertical scale for rcfs		*/
	char	labelling[50];

    /***Set the origin for labelling of block of radial correlations***/
	xshift    = txt_margin + 2.0*(1.0 + sep_margin);
	yshift    = txt_margin + 4.0*(1.0 + sep_margin);

    /***Insert the labelling***/
	color(3);
	hfont("times.rb");
	hcentertext(0);
	htextsize(txt_margin/3.0, txt_margin/2.2);
	htextang(0.0);
	sprintf(labelling, "Radial correlation functions");
	x1  = xshift;
	y1  = yshift;
	move2(x1, y1);
	hcharstr(labelling);

    /***Loop for each pair of species***/
	for (i=1; i<nspp; ++i)
	for (j=1; j<nspp; ++j)
	{
	    /***Set the origin for the block of rcf graphs***/
		rcf_block = 2.0;	/*Starts the block half way across window*/
		xshift    = txt_margin + rcf_block * (1.0 + sep_margin);
		yshift    = txt_margin +       2.0 * (1.0 + sep_margin);

	    /***Set the origin for the current rcf graph***/
		xshift = xshift +       (double)(i-1) * (rcf_block/(double)(nspp-1) + 0.5*sep_margin);
		yshift = yshift + (double)(nspp-1-j)  * (      2.0/(double)(nspp-1) + 0.5*sep_margin);

	    /***Draw the plot boundary***/
		x_lower = xshift - col_margin;
		y_lower = yshift - col_margin;
		x_upper = xshift + col_margin + rcf_block/(nspp-1.0);
		y_upper = yshift + col_margin +       2.0/(nspp-1.0);
		color(7);
		polymode(PYM_LINE);
		rectf(x_lower, y_lower, x_upper, y_upper);

	    /***Draw the line rcf = 1***/
		color(15);
		x1  = xshift;
		y1  = yshift + 1.0 * yscale;
		move2(x1, y1);
		x1  = xshift + (rcf_block/(nspp-1.0));
		y1  = yshift + 1.0 * yscale;
		draw2(x1, y1);

	    /***Draw the simulated radial correlations for species under optimization***/
		color(3);
		x1  = xshift;
		y1  = yshift + rcfs[0][1+(i-1)*(nspp-1)+(j-1)] * yscale;
		if (y1 > y_upper)	y1 = y_upper;
		move2(x1, y1);
		for (bin=0; bin<binmax; ++bin)
		{	x1  = xshift + (rcf_block/(nspp-1.0)) * rcfs[bin][0]/rcfs[binmax-1][0];
			y1  = yshift + rcfs[bin][1+(i-1)*(nspp-1)+(j-1)] * yscale;
			if (y1 > y_upper)	y1 = y_upper;
			draw2(x1, y1);
		}
	}
}


/****************************************************************************************************************
* This procedure displays the results of the stochastic simulation						*
****************************************************************************************************************/

void	display_vogl()
{
    /***Puts modifications to graphs into a backbuffer for synchronous updating***/
	backbuffer();

    /***Clear the vogl window***/
    /***(backbuffer only overwrites -- it doesn't start from blank screen)***/
	color(BLACK);
	clear();

    /***Display the following information in the vogl window***/
	display_status();
	display_simulation();
/*	display_tpss();
*/	display_rcfs();

    /***Move backbuffer to become the front buffer***/
	swapbuffers();

}
