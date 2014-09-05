#include "crlprop.h"


/* parse programm line arguments, set-up parameter, initialise and run the simulation */
int main(int argc, const char* argv[])
{
    int ret = 0;
int i;
    double energy_keV, distance1_m, f_m, distance2_m;
    struct parameters parameters;
    struct s2c s2c;
    struct insidecrl insidecrl;
    struct c2f c2f;
    int tag = 0;  //No read, create field
    
    if (argc != 5)
    {
	fprintf(stderr, "usage: %s <energy/keV> <distance1/m> <f/m> <distance2/m>\n", argv[0]);
	ret = -1;
	goto cleanup;
    }

    energy_keV  = atof(argv[1]); /* energy in keV */
    distance1_m = atof(argv[2]); /* distance from point-source to crl device in m */
    f_m         = atof(argv[3]); /* focal distance of individual lens in m */
    distance2_m = atof(argv[4]); /* distance from crl device to detector in m */

    parameters.xray.energy        = energy_keV;
    parameters.xray.wavelength    =   12.398/energy_keV*1e-10; /* conversion factor from keV to angstrom to metre */
    parameters.xray.wavenumber    =    2*M_PI / parameters.xray.wavelength;

    parameters.source.distance    = distance1_m;

    parameters.crl.f              = f_m;
    parameters.crl.aperture       =    1e-3; /* 1 mm */
    parameters.crl.separation     =   10e-6; /* 10 µm */
    parameters.crl.number         =    1;
    parameters.crl.deltafactor    =   3.53e-6;
    parameters.crl.R		  =  f_m*2.*M_PI*parameters.crl.deltafactor;

    parameters.detector.distance  = distance2_m;
    parameters.detector.number    = 1000;
    parameters.detector.width     =  100e-6;
    parameters.detector.intensity = malloc(parameters.detector.number * sizeof(double)); 
    print_parameters(&parameters, stdout);


    /* now propagate from point-source to CRL device */
    copy_xray(&parameters.xray,   &s2c.xray);
    copy_source(&parameters.source, &s2c.source);
    source_to_crl(&s2c, parameters.crl.aperture);			

i=0; if (i==0) tag=tag;
show_field(s2c.field, "s2c.field (main) after running source_to_crl");
printf("values[1]=%f \n", cimag(s2c.field->values[1]));

	
    /* propagate through the CRL device */   //(just phase-change according to lens width profile)
    /* TODO */
    copy_xray(&s2c.xray, &insidecrl.xray);
    copy_crl(&parameters.crl, &insidecrl.crl);

//IF NO READ
if (tag==0){
    insidecrl.field = malloc(sizeof(struct field));
    insidecrl.field->size = malloc(s2c.field->dimensions*sizeof(int));
    insidecrl.field->values = malloc(s2c.field->components*sizeof(complex double));
    copy_field(s2c.field, insidecrl.field);
}   
    crl_inside(&insidecrl,tag);

    /* finally, propagate from CRL device to focus */
    /* TODO */
    copy_xray(&insidecrl.xray, &c2f.xray);
    copy_crl(&parameters.crl, &c2f.crl);
    copy_detector(&parameters.detector,&c2f.detector);
    crl_to_focus(&c2f);


    /* write result to a file */
    /* TODO */

    /* call gnuplot script to generate plots */
    system ("gnuplot gpscript.gp");

cleanup:
    return ret;
}//main



/* print parameters to FILE* f; if f==NULL, print to stdout */
int print_parameters(struct parameters* para, FILE* f)
{
    if (f==NULL)
	f = stdout;

    fprintf(f, "Parameter Settings\n");
    fprintf(f, "===========================\n");
    fprintf(f, "\n");

    fprintf(f, "X-Ray Settings\n");
    fprintf(f, "---------------------------\n");
    fprintf(f, "energy:        %8.2f keV\n", para->xray.energy);
    fprintf(f, "wavelength:    %8.2f Å\n",   para->xray.wavelength*1e10);
    fprintf(f, "wavenumber:    %8.2f Å⁻¹\n", para->xray.wavenumber*1e-10);
    fprintf(f, "\n");

    fprintf(f, "Source Settings\n");
    fprintf(f, "---------------------------\n");
    fprintf(f, "distance:      %8.2f m\n",   para->source.distance);
    fprintf(f, "\n");

    fprintf(f, "CRL Device Settings\n");
    fprintf(f, "---------------------------\n");
    fprintf(f, "focal length:  %8.2f m\n",   para->crl.f);
    fprintf(f, "aperture:      %8.2f mm\n",  para->crl.aperture*1e3);
    fprintf(f, "separation:    %8.2f µm\n",  para->crl.separation*1e6);
    fprintf(f, "offset:        %8.2f µm\n",  para->crl.offset*1e6);
    fprintf(f, "deltafactor:   %8.2f µm\n",  para->crl.deltafactor*1e6);
    fprintf(f, "radius R:      %8.2f µm\n",  para->crl.R*1e6);
    fprintf(f, "number/lenses: %5d\n",       para->crl.number);
    fprintf(f, "\n");

    fprintf(f, "Detector Settings\n");
    fprintf(f, "---------------------------\n");
    fprintf(f, "distance:      %8.2f m\n",   para->detector.distance);
    fprintf(f, "number/pixels: %5d\n",       para->detector.number);
    fprintf(f, "width:         %8.2f µm\n",  para->detector.width*1e6);
    fprintf(f, "pixel size:    %8.2f nm\n",  para->detector.width/para->detector.number*1e9);
    fprintf(f, "\n");

    return 0;
}


/* copy struct xray elements */
int copy_xray(struct xray* in, struct xray* out)
{
    int ret = 0;

    if (in == NULL)
    {
	fprintf(stderr, "error: in points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    if (out == NULL)
    {
	fprintf(stderr, "error: out points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    out->energy     = in->energy;
    out->wavelength = in->wavelength;
    out->wavenumber = in->wavenumber;

cleanup:
    return ret;
}

/* copy struct crl elements */
int copy_crl(struct crl* in, struct crl* out)
{
    int ret = 0;

    if (in == NULL)
    {
	fprintf(stderr, "error: in points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    if (out == NULL)
    {
	fprintf(stderr, "error: out points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    out->f = in->f;
    out->aperture = in->aperture;
    out->separation = in->separation;
    out->number = in->number;
    out->deltafactor = in->deltafactor;
    out->R = in->R;

cleanup:
    return ret;
}

/* copy struct source elements */
int copy_source(struct source* in, struct source* out)
{
    int ret = 0;

    if (in == NULL)
    {
	fprintf(stderr, "error: in points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    if (out == NULL)
    {
	fprintf(stderr, "error: out points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    out->distance = in->distance;

cleanup:
    return ret;
}

/* copy struct detector elements */
int copy_detector(struct detector* in, struct detector* out)
{
    int ret = 0;

    if (in == NULL)
    {
	fprintf(stderr, "error: in points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    if (out == NULL)
    {
	fprintf(stderr, "error: out points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    out->distance = in->distance;
    out->number = in->number;
    out->width = in->width;
    out->intensity = in->intensity;

cleanup:
    return ret;
}

/* copy field */
int copy_field(struct field* in, struct field* out)
{
    int ret = 0;
//    int dim=in->dimensions;
//    int i;
    if (in == NULL)
    {
	fprintf(stderr, "error: in points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    if (out == NULL)
    {
	fprintf(stderr, "error: out points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    out->dimensions = in->dimensions;
    out->components = in->components;
    out->size = malloc(in->dimensions * sizeof(int));
    memcpy(out->size, in->size, in->dimensions * sizeof(int));
    out->values = malloc(in->components * sizeof(complex double));
    memcpy(out->values , in->values, in->components * sizeof(complex double));

cleanup:
    return ret;
}

/* propagate point-source to CRL device */
int source_to_crl(struct s2c* arg, double L)
{
    int ret = 0;

    complex double* u = NULL;		//array to create field
    int N = Nmin;			//Number of points 
    int n = 1;				//dimensions (components)
    int i;				//counter
    double delta=L/N;		        //discretisation grid space (sampling)
    double posy=0;			//lens plane vertical position(y) from axis, to give an offset
    double dz=arg->source.distance;     //distance z
    double dy=0;			//distance y 
    double A;				//Amplitude
    double ph;           		//phase
    double Re, Im;			//Real and Imaginary components of field
    double wvl=(arg->xray.wavelength);	//wavelength
    char* fname="entrance_plane.txt";   //file name

    struct field field;
    field.size=NULL;
    field.values=NULL;

    if (arg == NULL)
    {
	fprintf(stderr, "error: arg points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }


  /*Delta-N loop, max phase difference - Sets N and delta*/
  optimizeDelta(&N, &delta, L, wvl, dz, posy);
 
    u = malloc(N*sizeof(complex double));
    if (u == NULL)
    {
	fprintf(stderr, "error: malloc failed (%s:%d) [%s]\n", __FILE__, __LINE__, strerror(errno));
	ret = -1;
	goto cleanup;
    }

    /**/
    A=1;
    for (i=0;i<N;i++)
    {
        dy=-0.5*L+i*delta+posy;
	ph=getPhase(wvl,dy,dz);
    Re=A*cos(ph);
    Im=A*sin(ph);
    u[i]=Re+I*Im;    
    }

    /* copy u to field; first get some memory */
    field.dimensions = 1;
    field.size = malloc(field.dimensions*sizeof(int));

printf("Start allocation in arg s2c.field src_to_crl...\n");//EXAMPLE
    arg->field = malloc(sizeof(struct field));//EXAMPLE
    arg->field->size = malloc(field.dimensions*sizeof(int));//EXAMPLE
printf("%d int allocated to size \n", field.dimensions);


    /* calculate total number of points in all dimensions */
    for (i=0; i<field.dimensions; i++)
    {
	field.size[i] = N;	//as many components as dimensions
	n *= field.size[i];
    }

    /*copy n to field.components*/ 
    field.components = n;
    field.values = malloc(n*sizeof(complex double));
    arg->field->values = malloc(n*sizeof(complex double));
printf("allocation done\n\n");
    /* copy u to field.values */
    for (i=0; i<n;i++){
    field.values[i]=u[i];		
    }

    copy_field(&field, arg->field);
show_field(&field, "source_to_crl");
show_field(arg->field, "s2c.field inside source_to_crl function");
    /* now save the entrance field to a file */
    write_field_to_file(&field, fname);

cleanup:
    if (field.size)   free(field.size);   field.size=NULL;
    if (field.values) free(field.values); field.values=NULL;
    if (u) free(u); u=NULL;
    return ret;
}


/* propagate field through CRL device (apply phase-shift in real space)*/
//11111111111111111111111111111111111111111111111111111111111111111111111111111
int crl_inside(struct insidecrl* arg, int tag){
    int ret = 0;
      void *foo=arg;   int foo2=tag;  arg=foo;   foo2=foo2*0;

/*
    double y;				           //postition on lens
    double wvl=arg->xray.wavelength;		   //wavelength
    double L = arg->crl.aperture;      		   //aperture
    //double posy = arg->crl.offset;      	   //lens off-axis
    double R = arg->crl.R;		      	   //radius  
    double deltafactor = arg->crl.deltafactor;     //delta factor (phase shift)
    int i; 					   //counter
    int N=Nmin+0*tag;					   //Number of points 
    int n = 1;					   //dimensions (components)
    double* w;    				   //array holdind lens profile values
    double delta;
    double phshift;
    char* fnamewrite="crl_plane.txt";   	   //file name
    char* fnameread;
    tag=-1;
    bool noread=(tag<0);

    struct field field;
    field.size=NULL;
    field.values=NULL;

    if (arg == NULL)
    {
	fprintf(stderr, "error: arg points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }
    
    field.dimensions=1;
    field.size = malloc(field.dimensions*sizeof(int));
    
    fprintf(stderr, "warning: no scaling nor fft used in %s.\n", __FUNCTION__);
    if (noread) {
    printf("arg.components= %d \n", arg->field->components);
    field.values = malloc((arg->field->components)*sizeof(complex double));
    copy_field(arg->field, &field);
    } 
    
    else{ 
    fprintf(stderr, "warning: field is read from file in %s.\n \
    \t Problems (dim > 1) could appear \n\n", __FUNCTION__);

    fnameread="entrance_plane.txt";
    //Read values , get N and allocate memory field.values
    read_field_from_file(&field, fnameread);
    }
    n=field.components;

   //  retrieve size from total number of points in all dimensions 
    for (i=0; i<field.dimensions; i++)
    {
	//Problems with inverse process dim!=1. Usage of pow nth root is an option. 
        //field.size[i] = N;	//as many components as dimensions
	//n *= field.size[i];
      N=n;
    field.size[i]=N;
    }
    delta=L/N;

    // calculate crl profile and apply phase shift 
    w = (double*) malloc(N*sizeof(double));
    for (i=0;i<N;i++){
    y=-0.5*L+delta*i;
    w[i]=0.5*y*y/R;
    phshift=fmod(2.*M_PI*deltafactor*2.*w[i]/wvl,2*M_PI);

    field.values[i] *= cexp(I*phshift);
    }

    // now save the crl field to a file 
    write_field_to_file(&field, fnamewrite);

cleanup:
    if (field.size)   free(field.size);   field.size=NULL;
    if (field.values) free(field.values); field.values=NULL;
    if (w) free(w); w=NULL; */
    return ret;
}
/* propagate field from CRL to focus(detector) */
//22222222222222222222222222222222222222222222222222222222222222222222222222222
int crl_to_focus(struct c2f* arg){
    int ret=0;
      void *foo=arg;   arg=foo;   
    /*
    double y;				           //position on lens
    double wvl=arg->xray.wavelength;		   //wavelength
    double L = arg->crl.aperture;      		   //aperture
    // double posy = arg->crl.offset;      	   //lens off-axis
    int i; 					   //counter
    int N=Nmin;					   //Number of points 
    int n = 1;					   //dimensions (components)
    // double f = arg->crl.f;
    double distance = arg->detector.distance;
    double delta;
    double deltaf;
    double complex quadratic;
    double complex constphase;			   //constant phase
    double complex scaling;
    char* fname="det_plane.txt";   //file name

    fftw_complex* in, *out;
    fftw_plan p;
    fftw_complex val;

    struct field field;
    field.size=NULL;
    field.values=NULL;

    if (arg == NULL)
    {
	fprintf(stderr, "error: arg points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }
    field.dimensions=1;
    field.size = malloc(field.dimensions*sizeof(int));
   
    fprintf(stderr, "warning: scaling and correct use of fft to be reviewed in %s.\n", __FUNCTION__);
    fprintf(stderr, "warning: field is read from file in %s.\n \
    \t Problems (dim > 1) could appear \n\n", __FUNCTION__);

    //Read values, get N and allocate memory (field.size) inside function
    read_field_from_file(&field, "crl_plane.txt");
    n=field.components;

    // retrieve size from total number of points in all dimensions 
    for (i=0; i<field.dimensions; i++)
    {
	//Problems with inverse process dim!=1. Usage of pow nth root is an option. 
        //field.size[i] = N;	//as many components as dimensions
	//n *= field.size[i];
      N=n;
    field.size[i]=N;
    }
    delta=L/N;
    deltaf=1./(delta*N);

    in  = fftw_malloc(sizeof(fftw_complex) * n);
    out = fftw_malloc(sizeof(fftw_complex) * n);
    p   = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // fftw  and quadratic phase
   for (i=0; i<n; i++)
    {
        y=-0.5*L+delta*i;
        quadratic=cexp(I*M_PI*y*y/(wvl*distance));
	in[i]=quadratic*field.values[i];
    }

    //fftexecute
    fftw_execute(p);
    //fftshift
    if (fftshift)
    for (i=0; i<N/2; i++){
    val=out[i];
    out[i]=out[i+N/2];
    out[i+N/2]=val;
    }
  
    constphase = cexp(I*distance*2.*M_PI/wvl);
    scaling = 1./(I*wvl*distance);


    
  for (i=0; i<n; i++)
    {
        y=(-N/2+i)*deltaf*wvl*distance;
        quadratic=cexp(I*M_PI*y*y/(wvl*distance));
        field.values[i]=scaling*constphase*quadratic*delta*out[i];

    }

    write_field_to_file(&field, fname);  

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

cleanup:
    if (field.size)   free(field.size);   field.size=NULL;
    if (field.values) free(field.values); field.values=NULL;*/
    return ret;
}

int read_field_from_file(struct field* field, const char* fname)
{
    int ret = 0;
    int i = 0;
    int n;
    double Re,Im;
    double val;				//Value of intensity or phase
    FILE* f = fopen(fname, "r");

    if (field == NULL)
    {
	fprintf(stderr, "error: field points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }
    if (field->dimensions < 1)
    {
	fprintf(stderr, "error: field has invalid dimensionality %d (%s:%d)\n",
		field->dimensions, __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    if (fname == NULL)
    {
	fprintf(stderr, "error: fname points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    f = fopen(fname, "r");
    if (f == NULL)
    {
	fprintf(stderr, "error: cannot open file [%s] (%s:%d) [%s]\n",
		fname, __FILE__, __LINE__, strerror(errno));
	ret = -1;
	goto cleanup;
    }

    /* read values from FILE f, set n */ 
    n=0;
    while (fscanf(f, "%lf %lf %lf \n", &Re, &Im, &val )!= EOF ){  /* get components */    
    n++;
    }
    rewind(f);		//bring back to beginning to start the proper reading

    field->components=n;
    field->values = malloc(n*sizeof(complex double));

    for (i=0; i<n; i++){ /* loop through input-reading */    
    fscanf(f, "%lf %lf %lf \n", &Re, &Im, &val);  
    field->values[i]=Re+I*Im;
    }

cleanup:
    if (f) fclose(f); f=NULL;
    return ret;
}

int write_field_to_file(struct field* field, const char* fname)
{
    int ret = 0;
    int n=field->components;
    int i;
    FILE* f = fopen(fname, "w");

    if (field == NULL)
    {
	fprintf(stderr, "error: field points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }
    if (field->dimensions < 1)
    {
	fprintf(stderr, "error: field has invalid dimensionality %d (%s:%d)\n",
		field->dimensions, __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }
    if (field->size == NULL)
    {
	fprintf(stderr, "error: field has invalid size array (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }
    if (field->values == NULL)
    {
	fprintf(stderr, "error: field has invalid values array (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    if (fname == NULL)
    {
	fprintf(stderr, "error: fname points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    f = fopen(fname, "w");
    if (f == NULL)
    {
	fprintf(stderr, "error: cannot open file [%s] (%s:%d) [%s]\n",
		fname, __FILE__, __LINE__, strerror(errno));
	ret = -1;
	goto cleanup;
    }


    /* TODO:
(done)     * write field to FILE* f
     */
    for (i=0;i<n;i++){
      if (print_intens) fprintf(f, "%f %f %f \n", creal(field->values[i]), cimag(field->values[i]), cabs(field->values[i]));
      else fprintf(f, "%f %f %f \n", creal(field->values[i]), cimag(field->values[i]), carg(field->values[i])); 
    }

cleanup:
    if (f) fclose(f); f=NULL;
    return ret;
}



double getPhase(double wvl, double dy, double dz){
  double r=sqrt(dy*dy+dz*dz);				//distance
  double opd=r-dz;					//optical path difference
  return fmod(opd*2.*M_PI/wvl,2*M_PI);			//might have problems if long int is not used
}

bool optimizeDelta(int* Nref, double* deltaref, double L, double wvl, double dz, double posy){
double deltastep=0.5e-6;
double dy;
double ph;
double prevph;
double phdif;
int N = *Nref;
double delta = *deltaref;
int i,j;
bool delta_is_set=false;
bool delta_reached_max=false;

while (!delta_is_set){
   j=0;
        while (delta>=delta_max){
        N=2*N;
        delta=L/N;
        }
   delta_reached_max=false;	  
   while(!delta_reached_max){
     /////////////////////////////phdif calculation//////////////////////////////////////////////////////
     phdif=2*M_PI;
     ph=0;
	  
     for(i=0;i<2;i++){
     prevph=ph;
     dy=delta*(N/2-i)+fabs(posy);
     ph=getPhase(wvl,dy,dz);
     }

      phdif=fabs(ph-prevph);
      phdif=(phdif<(2.*M_PI-phdif)) ? phdif : 2.*M_PI-phdif;  		
      ///////////////////////////////////////////////////////////////////////////////////////////////////
	  
      if(phdif<(0.1*M_PI)){
      delta=delta+deltastep;
      delta_reached_max=false; 		//keep iterating
      j++;		     		//suitable delta must have j!=0
      } 

      else if((phdif>=(0.1*M_PI))&&(j!=0)) {
      delta_reached_max=true;
      delta=L/N+(j-1)*deltastep;	//take previous value
      delta_is_set=(delta<delta_max);	//check, only successful case
      }
	
      else {
      N=2*N;
      delta=L/N;
      delta_reached_max=true;
      }	  
   }
}
    *deltaref=delta;
    *Nref=N;
    printf("Calculation parameters\n");
    printf("---------------------------\n");
    printf("N:             %5d\n",       N);
    printf("delta:         %8.2f µm\n",  delta*1e6);    
    printf("\n");
return true;
}
/*/Get field from file
get_field(struct &inside_crl, char* fname);    	
	
	    fprintf(stderr, "warning: field is read from file in %s.\n \
	    \t Problems (dim > 1) could appear \n\n", __FUNCTION__);
	
	    //Read values, get N and allocate memory (field.size) inside function
	    read_field_from_file(arg->field, fname);
	    n=field.components;
	
	    // retrieve size from total number of points in all dimensions 
	    for (i=0; i<field.dimensions; i++)
	    {
		//Problems with inverse process dim!=1. Usage of pow nth root is an option. 
	        //field.size[i] = N;	//as many components as dimensions
		//n *= field.size[i];
	      N=n;
	    field.size[i]=N;
	    }
    	
	    //Memory allocation
	    field.values = malloc(n*sizeof(complex double));
	    arg->field = malloc(sizeof(struct field));
	    arg->field->size = malloc(field.dimensions*sizeof(int));
	    arg->field->values = malloc(n*sizeof(complex double));
    	}*/

//Show field parameters
int show_field(struct field* argfield, char* message){
int ret=0;
int i;
printf("   ******************************\n");
printf("   Field in %s \n", message);
printf("   argfield.dimensions=%d \n", argfield->dimensions);
printf("   argfield.size[0]=%d \n", argfield->size[0]);
printf("   argfield.components=%d \n", argfield->components);
printf("      first 10 {Re} values: \n");
for(i=0;i<10;i++){
printf("   argfield.values[%d]=%f\n",i,creal(argfield->values[i]));
}
printf("\n");
return ret;
}

