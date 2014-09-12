//Code including functions:
//read_field_from_file
//write_field_to_file
//optimizeDelta
//getPhase
#include "crlprop.h"

int read_field_from_file(struct field* field, const char* fname)
{
    int ret = 0;
    int i = 0;
    int n;
    double Re,Im,y;
    double val, phval;				//Value of intensity or phase
    FILE* f = fopen(fname, "r");

    if (field == NULL)
    {
	fprintf(stderr, "error: field points to NULL (%s:%d)\n", __FILE__, __LINE__);
	ret = -1;
	goto cleanup;
    }

    //Set dimensions and allocate size
    field->dimensions = field_dim;	//1 dimension
    field->size = malloc(field->dimensions*sizeof(int));
    fprintf(stderr, "warning: field is read from file in %s.\n \
    \t Problems (dim > 1) could appear \n\n", __FUNCTION__);
   
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
    while (fscanf(f, "%lf %lf %lf %lf %lf \n", &y, &Re, &Im, &val, &phval )!= EOF ){  /* get components */    
    n++;
    }
    rewind(f);		//bring back to beginning to start the proper reading

    //Set components and allocate values
    field->components=n;
    field->values = malloc(n*sizeof(complex double));

    for (i=0; i<n; i++){ /* loop through input-reading */    
    fscanf(f, "%lf %lf %lf %lf %lf \n", &y, &Re, &Im, &val, &phval);  
    field->values[i]=Re+I*Im;
    }

cleanup:
    if (f) fclose(f); f=NULL;
    return ret;
}

int write_field_to_file(struct field* field, const char* fname, double L)
{
    int ret = 0;
    int n=field->components;
    int i;
    FILE* f = fopen(fname, "w");
    double A;

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
      A=cabs(field->values[i]);
      //intensity/amplitude
      if (print_intens) fprintf(f, "%f %f %f %f %f \n", (-0.5+i*1./n)*L, creal(field->values[i]), cimag(field->values[i]), A*A, carg(field->values[i])); //A*A

      //phase
      else fprintf(f, "%f %f %f %f \n", (-0.5+i*1./n)*L, creal(field->values[i]), cimag(field->values[i]), carg(field->values[i])); 
    }

cleanup:
    if (f) fclose(f); f=NULL;
    return ret;
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

printf("####################\n");
printf("Delta-N optimization\n");
printf("N=%d delta=%8.2f µm\n",N, delta*1e6);

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
/*      delta=delta+deltastep;						AVOID APPLYING CHANGES (CONFLICT WITH PAD)
printf("add N=%d delta=%8.2f µm\n",N, delta*1e6);
      delta_reached_max=false; 		//keep iterating*/
delta_reached_max=true;
delta_is_set=true;
printf("WARNING: Optimize_delta changes not applied, delta<delta_max phdif<0.1pi\n");
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
printf("WARNING: WRONG N VALUES\n");
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
double getPhase(double wvl, double dy, double dz){
  double r=sqrt(dy*dy+dz*dz);				//distance
  double opd=r-dz;					//optical path difference
//printf("getPhase opd=%1.6fA ph=%e ph(2π)=%1.6f \n", opd*1.0e10, opd/wvl*2.*M_PI, fmod(opd*2.*M_PI/wvl,2*M_PI));
  return fmod(opd*2.*M_PI/wvl,2*M_PI);			//might have problems if long int is not used
}

int pad_field(struct field* arg){
int ret=0;
int i;
int N=arg->components;

for (i=0;i<N/4;i++){
arg->values[i]=0.;
}
for (i=3*N/4;i<N;i++){
arg->values[i]=0.;
}
return ret;
}

double getPhase_A(double wvl_A, double dy, double dz){
  double r=sqrt(dy*dy+dz*dz);				//distance
  double opd_A=(r-dz)*1.0e10;					//optical path difference
printf("getPhase_A opd=%1.6fA ph=%e ph(2π)=%1.6f \n", opd_A, opd_A/wvl_A*2.*M_PI, fmod(opd_A/wvl_A*2.*M_PI,2*M_PI));
  return fmod(opd_A/wvl_A*2.*M_PI,2*M_PI);			//might have problems if long int is not used
}


/*int get_field(struct field* field, const char* fnameread){
int ret=0;
  field->dimensions = field_dim;	//1 dimension
  field->size = malloc(field->dimensions*sizeof(int));
  fprintf(stderr, "warning: field is read from file in %s.\n \
  \t Problems (dim > 1) could appear \n\n", __FUNCTION__);
	
  //fnameread="entrance_plane.txt";
  //Read values, get N and implic. allocate memory (field.values) inside function field
  read_field_from_file(field, fnameread);
	    
  //arg->field=malloc(sizeof(struct field));
  //arg->field->size=malloc(field.dimensions*sizeof(int));
  //arg->field->values=malloc(field.components*sizeof(complex double));
	    
  //Once the field has been read, turn of reading TO BE DONE IN MAIN
  // (*tag)=-1;
  return ret;
  }*/
