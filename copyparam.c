//Code including functions: 
//print_parameters
//copy_xray
//copy_crl
//copy_detector
//copy_field
//show_field
#include "crlprop.h"

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
    fprintf(f, "wavelength:    %8.4f Å\n",   para->xray.wavelength*1e10);
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
    fprintf(f, "width:         %8.4f m\n",  para->detector.width);
    fprintf(f, "pixel size:    %8.2f µm\n",  para->detector.width/para->detector.number*1e6);
    fprintf(f, "Airy spot r=1.22*wvl*2*z/L:    %8.2f µm\n",  1.22*para->xray.wavelength*para->detector.distance*2./para->crl.aperture*1e6);
    fprintf(f, "Airy spot r=1.22*wvl*z/L:    %8.2f µm\n",  1.22*para->xray.wavelength*para->detector.distance/para->crl.aperture*1e6);
//    fprintf(f, "number/pixels: to be set by field prop\n");
//    fprintf(f, "width: to be set by field prop\n");
//    fprintf(f, "pixel size: to be set by field prop\n");
    fprintf(f, "\n");

    fprintf(f, "Simulation N, delta \n");
    fprintf(f, "---------------------------\n");
    fprintf(f, "N defined(header): %5d\n",       Nmin);
    fprintf(f, "delta=%1.6fμm       \n",  para->crl.aperture/Nmin*1.e6);
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
//    out->size = malloc(in->dimensions * sizeof(int));
//    out->values = malloc(in->components * sizeof(complex double));
    memcpy(out->size, in->size, in->dimensions * sizeof(int));
    memcpy(out->values , in->values, in->components * sizeof(complex double));

cleanup:
    return ret;
}

//Show field parameters
int show_field(struct field* argfield, char* message){
int ret=0;
int i;
printf("   ******************************\n");
printf("   Field in %s \n", message);
printf("   argfield.dimensions=%d \n", argfield->dimensions);
printf("   argfield.size[0]=%d \n", argfield->size[0]);
printf("   argfield.components=%d \n", argfield->components);
printf("      first 10 {Re} {Im} values: \n");
for(i=0;i<10;i++){
printf("   argfield.values[%d]=%f,%f",i,creal(argfield->values[i]), cimag(argfield->values[i]));
printf("   \t %1.2f*exp(%1.4fi)\n",cabs(argfield->values[i]), carg(argfield->values[i]));
}
printf("\n");
return ret;
}

