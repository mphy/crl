#include "crlprop.h"


/* parse programm line arguments, set-up parameter, initialise and run the simulation */
int main(int argc, const char* argv[])
{
    int ret = 0;

    double energy_keV, distance1_m, f_m, distance2_m;
    struct parameters parameters;
    //argument structures for funtions and field memory allocation 
    struct s2c s2c; 			s2c.field = malloc(sizeof(struct field));
    struct insidecrl insidecrl; 	insidecrl.field = malloc(sizeof(struct field));
    struct c2f c2f;			c2f.field = malloc(sizeof(struct field));

    int tag = tag_val;           //tag determines whether initially create field or read from file
    int current=0;	         //int used to keep track of the step the simulation goes through (identifying propagation)
    double distance;		 //distance to perform fourier propagation
    char fnameread[20];		 //name of the file to be read
    char fnamewrite[20];	 //name of the file to wite to

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
    parameters.crl.number         =   lenses;
    parameters.crl.deltafactor    =   3.53e-6;
    parameters.crl.R		  =  f_m*2.*M_PI*parameters.crl.deltafactor;

    parameters.detector.distance  = distance2_m;// (40/39 is q for s=40 and f=1);
    parameters.detector.number    = 1;		//to be defined by field propagation
    parameters.detector.width     = 1;		//to be defined by field propagation
//    print_parameters(&parameters, stdout);


    /* If noread then copy parameters to s2c and generate field*/
    if (tag==-1){
    copy_xray(&parameters.xray,   &s2c.xray);
    copy_source(&parameters.source, &s2c.source);
    source_to_crl(&s2c, parameters.crl.aperture);
////    printf("FIRST STEP:  Create field tag=%d \n", tag);
    if (padf) pad_field(s2c.field);
    write_field_to_file(s2c.field, "entrance_plane.txt", parameters.crl.aperture);
    }			

    /* copy parameters to insidecrl structure*/
    copy_xray(&parameters.xray, &insidecrl.xray);
    copy_crl(&parameters.crl, &insidecrl.crl);

    /* copy parameters to c2f structure for the propagation*/
    copy_xray(&insidecrl.xray, &c2f.xray);
    copy_crl(&parameters.crl, &c2f.crl);
    copy_detector(&parameters.detector,&c2f.detector); 

    //Pre-process: Create/read field and allocate memory. Starting point will be after the lens
    //if field is coming from entrance_plane (both generated or read) allocate and execute PREVIOUS phase-shifting 
    if (tag>0){
	sprintf(fnameread,"crl_plane");
	sprintf(fnameread,"%s%d.txt", fnameread, tag);
	read_field_from_file(insidecrl.field, fnameread);
	tag=-1;//CHANGE (fourier propagation)
        } 

	//Get entrance plane field	  
	else{
	   if (tag==-1){	
	    //Getting from field already generated
            insidecrl.field->size = malloc(s2c.field->dimensions*sizeof(int));
            insidecrl.field->values = malloc(s2c.field->components*sizeof(complex double));
            copy_field(s2c.field, insidecrl.field);  
	    }

	    //Reading from entrance_plane
            else if (tag==0){
	    sprintf(fnameread,"entrance_plane");
	    sprintf(fnameread,"%s.txt", fnameread);
	    read_field_from_file(insidecrl.field, fnameread);
	    tag=-1;
 	    }
////        printf("CRL STEP:  Apply phase-shift tag=%d \n", tag);
     	sprintf(fnamewrite,"crl_plane");
	sprintf(fnamewrite,"%s1.txt", fnamewrite);
        crl_inside(&insidecrl, fnamewrite);
	} 

    //Once N is set, get detector parameters
    parameters.detector.number=pix_num; //previously parameters.detector.number=insidecrl.field->components;
//    parameters.detector.width=parameters.detector.number*c2f.xray.wavelength*c2f.detector.distance/parameters.crl.aperture;
    parameters.detector.width=0.01;
    parameters.detector.intensity = malloc(parameters.detector.number * sizeof(double)); 

    copy_detector(&parameters.detector,&c2f.detector); 
    print_parameters(&parameters, stderr);

    // Start lens loop: while current state (lens) has not reached total lens amount. (phshift) + /PROP+(PHSHIFT)/
    //										      (previous) +/LOOP+(post)/
	while(current!=lenses){

////	   printf("NEW STEP:  current state=%d , tag=%d \n", current, tag);

	       // get field from previous steps, allocate if first time
	       if (current==0){
	           c2f.field->size = malloc(insidecrl.field->dimensions*sizeof(int));
	           c2f.field->values = malloc(insidecrl.field->components*sizeof(complex double));
	           }
	       copy_field(insidecrl.field, c2f.field);


//RUN 	   //perform fourier propagation: from CRL to CRL or Detector (different methods)
////	   printf("Call propagation... \n");
   	   if (current==(lenses-1)){
	   distance=parameters.detector.distance;
//           get_spherical(&c2f, tag, distance);
	   if (padf) pad_field(c2f.field);

	   if (ray_add) ray_addition(&c2f, tag, distance);
           else crl_to_focus(&c2f, tag, distance);
	   }
	   else{
           distance=parameters.crl.separation;  
           crl_to_fourier(&c2f, tag, distance);
	   }
	  

////    	   printf("distance for propagation %d was %f \n", current, distance);
////	   printf("propagation successfully done. \n");
	   copy_field(c2f.field, insidecrl.field);

	   //if not reached detector then phase-shift before propagation 
	   if (current!=(lenses-1)){
	   sprintf(fnamewrite,"crl_plane");
	   sprintf(fnamewrite,"%s%d.txt", fnamewrite, current+2);
           crl_inside(&insidecrl, fnamewrite);
	   }

	   current=current+1;
      }
      write_field_to_file(c2f.field, "det_plane.txt", parameters.detector.width);

    
    /* call gnuplot script to generate plots */
    system ("gnuplot gpscript.gp");

cleanup:
    if (parameters.detector.intensity)   free(parameters.detector.intensity);   parameters.detector.intensity=NULL;
    return ret;
}//main


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
    char* fnamewrite="entrance_plane.txt";   //file name

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
  //optimizeDelta(&N, &delta, L, wvl, dz, posy);
  N=Nmin; delta=L/N;
  double phdif=fabs(getPhase(wvl,delta*(N/2),dz)-getPhase(wvl,delta*(N/2-1),dz));
  phdif=(phdif<(2.*M_PI-phdif)) ? phdif : 2.*M_PI-phdif; 
  if (phdif>=0.1*M_PI) fprintf(stderr, "warning: Wrong N, delta values phdif=%f \n", phdif);
      
////printf("N=%d  delta=%1.6fμ 	Manually set\n", N, delta*1e6);
//printf("Press 'Enter' to continue: ... "); while ( getchar() != '\n') ; 
   
 
    u = malloc(N*sizeof(complex double)); //dim 1
    if (u == NULL)
    {
	fprintf(stderr, "error: malloc failed (%s:%d) [%s]\n", __FILE__, __LINE__, strerror(errno));
	ret = -1;
	goto cleanup;
    }


	//Gaussian Amplitude
	double sigma=L/3./2.;
    for (i=0;i<N;i++)
    {
        dy=-0.5*L+i*delta+posy;
	ph=getPhase(wvl,dy,dz);

    A=exp(-dy*dy/2./sigma/sigma);

    Re=A*cos(ph);
    Im=A*sin(ph);
    u[i]=Re+I*Im;    
    }

    /* copy u to field; first get some memory */
    field.dimensions = field_dim;
    field.size = malloc(field.dimensions*sizeof(int));

    //Start allocation for s2c
    arg->field->size = malloc(field.dimensions*sizeof(int));


    /* calculate total number of points in all dimensions */
    for (i=0; i<field.dimensions; i++)
    {
	field.size[i] = N;	//as many components as dimensions
	n *= field.size[i];
    }

    /*copy n to field.components and allocate both for field and s2c.field*/ 
    field.components = n;
    field.values = malloc(n*sizeof(complex double));
    arg->field->values = malloc(n*sizeof(complex double));

    /* copy u to field.values */
    for (i=0; i<n;i++){
    field.values[i]=u[i];		
    }

    copy_field(&field, arg->field);

    /* now save the entrance field to a file */
    write_field_to_file(&field, fnamewrite, L);

cleanup:
    if (field.size)   free(field.size);   field.size=NULL;
    if (field.values) free(field.values); field.values=NULL;
    if (u) free(u); u=NULL;
    return ret;
}


/* propagate field through CRL device (apply phase-shift in real space)*/
//11111111111111111111111111111111111111111111111111111111111111111111111111111
int crl_inside(struct insidecrl* arg, char* fnamewrite){
    int ret = 0;


    double y;				           //postition on lens
    double wvl=arg->xray.wavelength;		   //wavelength
    double L = arg->crl.aperture;      		   //aperture
    //double posy = arg->crl.offset;      	   //lens off-axis
    double R = arg->crl.R;		      	   //radius  
    double deltafactor = arg->crl.deltafactor;     //delta factor (phase shift)
    int i; 					   //counter
    int N=Nmin;					   //Number of points 
    int n = 1;					   //dimensions (components)
    double* w;    				   //array holdind lens profile values
    double delta;
    double phshift;

    struct field field;
    field.size=NULL;
    field.values=NULL;

    if (arg == NULL)
    {
	fprintf(stderr, "error: arg points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }
    
    fprintf(stderr, "warning: equal CRL profiles are used and computed for each step.\n");
    fprintf(stderr, "warning: no scaling nor fft used in %s.\n", __FUNCTION__);
    
	  //Allocate memory for the field
          field.size = malloc(arg->field->dimensions*sizeof(int));
    	  field.values = malloc((arg->field->components)*sizeof(complex double)); //explicit allocation field
	  copy_field(arg->field, &field);

    n=field.components;

   //  retrieve size from total number of points in all dimensions 
    for (i=0; i<field.dimensions; i++)
    {
	//Problems with inverse process dim!=1. Usage of pow nth root is an option. 
        //field.size[i] = N;	//as many components as dimensions
	//n *= field.size[i];
      N=n; //dim 1
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

    // now save the crl field to a file and structure
    write_field_to_file(&field, fnamewrite, L);
    copy_field(&field, arg->field);
    //show_field(&field, N/2, 10, "in crl_inside after shifting");

cleanup:
    if (field.size)   free(field.size);   field.size=NULL;
    if (field.values) free(field.values); field.values=NULL;
    if (w) free(w); w=NULL; 
    return ret;
}


/* propagate field from CRL using Fourier Transform (to crl!) */
//22222222222222222222222222222222222222222222222222222222222222222222222222
int crl_to_fourier(struct c2f* arg, int tag, double distance){
    int ret=tag;

    
    double y;				           //position on lens
    double wvl=arg->xray.wavelength;		   //wavelength
    double L = arg->crl.aperture;      		   //aperture
    // double posy = arg->crl.offset;      	   //lens off-axis
    int i; 					   //counter
    int N=Nmin;					   //Number of points 
    int n = 1;					   //dimensions (components)
    double delta;				   //grid spacing on lens
    double deltaf;				   //freq grid spacing
    double complex quadratic;			   //quadratic phase factor
    double complex constphase;			   //constant phase
    double complex scaling;			   //normalisation (scaling) fft

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
 
   fprintf(stderr, "warning: scaling and correct use of fft to be reviewed in %s.\n", __FUNCTION__);

     field.size = malloc(arg->field->dimensions*sizeof(int));
     field.values = malloc((arg->field->components)*sizeof(complex double)); //explicit allocation field
     copy_field(arg->field, &field);
    
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
printf("CRL to CRL propagation N=%d delta=%fμ deltaf=%f \n", N, delta*1e6, deltaf); 
    in  = fftw_malloc(sizeof(fftw_complex) * n);
    out = fftw_malloc(sizeof(fftw_complex) * n);
    p   = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // fftw  and quadratic phase (y1)
    for (i=0; i<n; i++)
    {
        y=-0.5*L+delta*i;
        quadratic=cexp(I*M_PI*y*y/(wvl*distance));
	in[i]=quadratic*field.values[i];
    }

/*    //previous fftshift
    if (fftshift)
    for (i=0; i<N/2; i++){
    val=in[i];
    in[i]=in[i+N/2];
    in[i+N/2]=val;
    }
*/
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
    scaling=sqrt(1./(wvl*distance));
//    scaling = 1./(I*wvl*distance);
//    scaling = 1.;
//    scaling = cpow(1./(I*wvl*distance), 0.5);

    //apply quadratic phase (y2), scaling and constant phase
  for (i=0; i<n; i++)
    {
        y=(-0.5*N+i)*deltaf*wvl*distance;
        quadratic=cexp(I*M_PI*y*y/(wvl*distance));
        field.values[i]=scaling*constphase*quadratic*delta*out[i];

    }

//    write_field_to_file(&field, fnamewrite);  done in main
      copy_field(&field, arg->field);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

cleanup:
    if (field.size)   free(field.size);   field.size=NULL;
    if (field.values) free(field.values); field.values=NULL;
    return ret;
}

/* propagate field from CRL using Fourier Transform (to plane!) */
//33333333333333333333333333333333333333333333333333333333333333333333333333
int crl_to_focus(struct c2f* arg, int tag, double distance){
    int ret=tag;

    
    double y;				           //position on final lens/detector
    double wvl=arg->xray.wavelength;		   //wavelength
    double L = arg->crl.aperture;    		   //size of lens
    int i; 					   //counter
    int N=Nmin;					   //Number of points 
    int n = 1;					   //dimensions (components)
    double delta;				   //grid spacing on lens
    double deltaf;				   //freq grid spacing
    double complex quadratic;			   //quadratic phase factor
    double complex constphase;			   //constant phase
    double complex scaling;			   //normalisation (scaling) fft

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
   fprintf(stderr, "warning: detector members are not used.\n");
   fprintf(stderr, "warning: scaling and correct use of fft to be reviewed in %s.\n", __FUNCTION__);

     field.size = malloc(arg->field->dimensions*sizeof(int));
     field.values = malloc((arg->field->components)*sizeof(complex double)); //explicit allocation field
     copy_field(arg->field, &field);
    
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
printf("CRL to focus propagation N=%d delta=%fμ deltaf=%f \n", N, delta*1e6, deltaf); 
    in  = fftw_malloc(sizeof(fftw_complex) * n);
    out = fftw_malloc(sizeof(fftw_complex) * n);
    p   = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // fftw  and quadratic phase (y1)
    for (i=0; i<n; i++)
    {
        y=(-0.5*N+i)*delta;
        quadratic=cexp(I*M_PI*y*y/(wvl*distance));
	in[i]=quadratic*field.values[i];
    }

/*    //previous fftshift
    if (fftshift)
    for (i=0; i<N/2; i++){
    val=in[i];
    in[i]=in[i+N/2];
    in[i+N/2]=val;
    }
*/
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
      scaling=sqrt(1./(wvl*distance));
//    scaling = 1./(I*wvl*distance);
//    scaling = 1.;
//    scaling = cpow(1./(I*wvl*distance), 0.5);



  //apply quadratic phase (y2), scaling and constant phase
  for (i=0; i<n; i++)
    {
        y=(-0.5*N+i)*deltaf*wvl*distance;
        quadratic=cexp(I*M_PI*y*y/(wvl*distance));
        field.values[i]=scaling*constphase*quadratic*delta*out[i];

    }

//    write_field_to_file(&field, fnamewrite);  done in main
      copy_field(&field, arg->field);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

cleanup:
    if (field.size)   free(field.size);   field.size=NULL;
    if (field.values) free(field.values); field.values=NULL;
    return ret;
}


//ray addition			(field not used in recent implementation, changing)
int ray_addition(struct c2f* arg, int tag, double dz){
int ret=tag;
double wvl=arg->xray.wavelength;	   //wavelength
int i,j; 				   //counter
int N=arg->field->components;		   //Number of points
double L=arg->crl.aperture;		   //Lens aperure
int pix=arg->detector.number;		   //pixels in detector
double w=arg->detector.width;	           //detector width
double ph;				   //phase
double A;				   //amplitude
double y2;				   //position on detector
double y;				   //position on lens
complex double val;			   //complex value to add intensities

    struct field field;
    field.size=NULL;
    field.values=NULL;
    fprintf(stderr, "warning: using %s.\n", __FUNCTION__);

     field.size = malloc(arg->field->dimensions*sizeof(int));
     field.values = malloc((arg->field->components)*sizeof(complex double)); 
     copy_field(arg->field, &field);
////if (ray_add) printf("ray_add det.width=%fm pix=%d det_delta=%f \n ", w, pix, w/pix); //print

     N=field.components;	//     delta=L/N;     //     deltaf=1./(delta*N);


	for (j=0;j<pix;j++){
	val=0.;
	y2=(-0.5*pix+j)*w/pix;
	      for (i=0;i<N;i++){
	        A=cabs(arg->field->values[i]); 
	        y=((-0.5*N)+i)*L/N;
	        //d=taylor_sqrt_opd((y2-y),dz,5));   //d=sqrt((y1[j]-y)*(y1[j]-y)+dz*dz)-dz;
	        ph=getPhase(wvl,(y2-y),dz);
                //ADDITION
                val=val+A*cexp(I*ph); 

		}//loop on lens points
         arg->detector.intensity[j]=cabs(val)*cabs(val);


         }//loop on detector pix

copy_field(&field, arg->field);
write_intensity_file(arg->detector.intensity, "int.txt", pix, w);

if (field.size)   free(field.size);   field.size=NULL;
if (field.values) free(field.values); field.values=NULL;
return ret;
}


//get spherical
int get_spherical(struct c2f* arg, int tag, double dz){
    int ret=tag;
    int i;
    double wvl=arg->xray.wavelength;	   //wavelength
    double N=arg->field->components;
    double L=arg->crl.aperture;
    double y;
    double d;
    double A;
    fprintf(stderr, "warning: using %s.\n", __FUNCTION__);

    struct field field;
    field.size=NULL;
    field.values=NULL;

     field.size = malloc(arg->field->dimensions*sizeof(int));
     field.values = malloc((arg->field->components)*sizeof(complex double)); 
     copy_field(arg->field, &field);

	//Gaussian Amplitude
	double sigma=L/3./2.;
        for (i=0;i<N;i++){
        y=(-0.5*N+i)*L/N;
        d=taylor_sqrt_opd(y,dz,5);
        A=exp(-y*y/2./sigma/sigma);


        //field.values[i]=A*cexp(I*((d-dz)*2.*M_PI/wvl));
field.values[i]=A*cexp(I*fmod(d*2.*M_PI/wvl,2.*M_PI));

}
write_field_to_file(&field, "spherical.txt", arg->crl.aperture);

copy_field(&field, arg->field);
if (field.size)   free(field.size);   field.size=NULL;
if (field.values) free(field.values); field.values=NULL;
return ret;
}
