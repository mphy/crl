//Code including functions:
//read_field_from_file
//write_field_to_file
//optimizeDelta
//getPhase
//write_intensity_file
//pad_field
//taylor_sqrt_opd

#include "crlprop.h"

int write_intensity_file(double* intensity, const char* fname, double pix, double width)
{
int ret=0;
FILE* f = fopen(fname, "w");
int i;
for (i=0;i<pix;i++){
fprintf(f, "%f %f \n", (-0.5+i*1./pix)*width, intensity[i]);
}
if (f) fclose(f); f=NULL;
return ret;
}

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
      //field values
      fprintf(f, "%f %f %f %f %f \n", (-0.5+i*1./n)*L, creal(field->values[i]), cimag(field->values[i]), A*A, carg(field->values[i])); 
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

////printf("####################\n");
////printf("Delta-N optimization\n");
////printf("N=%d delta=%8.2f µm\n",N, delta*1e6);

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
////    printf("Calculation parameters\n");
////    printf("---------------------------\n");
////    printf("N:             %5d\n",       N);
////    printf("delta:         %8.2f µm\n",  delta*1e6);    
////    printf("\n");
return true;
}

double getPhase(double wvl, double dy, double dz){
  double opd=taylor_sqrt_opd(dy,dz,5); 		//optical path difference
  return fmod(opd*2.*M_PI/wvl,2*M_PI);		//treat precision carefully
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

double taylor_sqrt_opd(double y, double z, int ord){
  double result=0.0;
  double prod=1.0;
  double sum=0.0;
  double x=1.0;
  int i;

    if (ord==0) result=0.;
     else{
	  for(i=1;i<=ord;i++){
  	  x=x*(y/z)*(y/z);
          prod=prod*(-1.)*(2.*i)*(2.*i-1.)/(4.*i*i);
	  sum=sum+prod*x/(1.-2.*i);
	  }  
	  result=sum;

     }
  return result*z;
}


