#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <errno.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>

#define Nmin 1024
#define delta_max 10e-6
#define field_dim 1
#define fftshift true
#define print_intens true
#define lenses 1
#define tag_val -1	//tag_val to turn on/off getting field from file
			//Use: (-1 --> noread) (0 --> get from entrance_plane) (1 --> get from crl_1)...    


/* some general structures to store parameters */

/* x-ray energy / wavelength */
struct xray {
    double energy;     /* x-ray energy, in keV */
    double wavelength; /* x-ray wavelength, in metre */
    double wavenumber; /* x-ray wavenumber, in inverse metre */
};


/* source parameters */
struct source {
    double distance; /* distance from source to crl device (entrance), in metre */
};

/* crl device parameters */
struct crl {
    double f;            /* focal width of individual lens, in metre */
    double aperture;     /* lateral size of lenses, diametre, in metre */
    double separation;   /* separation distance between lenses, in metre */
    double offset;	 /* vertical position offset */
    int number;          /* number of lenses inside this crl device */
    double deltafactor;  /* delta factor (phase-shifting)*/
    double R;            /* radius */
};

/* detector parameters */
struct detector {
    double distance;   /* distance from crl device (exit), in metre */
    int number;        /* number of pixels */
    double width;      /* width / lateral size of detector, in metre */
    double* intensity; /* array to the propagated intensity */
};

/* collector for all relevant parameter structs */
struct parameters {
    struct xray xray;
    struct source source;
    struct crl crl;
    struct detector detector;
};


/* struct holding a complex-valued field */
struct field {
    int dimensions;         /* number of dimensions */
    int* size;              /* array holding number of points per dimension */
    complex double* values; /* array holding complex valued field */
    int components;         /* number of components */ 
};


/* parameters to propagate from point-source to crl device */
struct s2c {
    struct xray xray;
    struct source source;
    struct field* field;
};

/* parameters to the crl inside propagation*/
struct insidecrl{
    struct xray xray;
    struct crl crl;
    struct field* field;
};

/* parameters to the detector propagation*/
struct c2f{
    struct xray xray;
    struct crl crl;
    struct field* field;          
    struct detector detector;
};

/* parameters to the grid simulation*/
/*struct grid{
    int N;
    int dimensions;
    double delta;
}
*/
int source_to_crl(struct s2c* arg, double L);
int crl_inside(struct insidecrl* arg, char* fnamewrite);
int crl_to_fourier(struct c2f* arg, int tag, double distance);
int crl_to_focus(struct c2f* arg, int tag, double distance);

int print_parameters(struct parameters* para, FILE* f);
int write_field_to_file(struct field* field, const char* fname, double L);
int read_field_from_file(struct field* field, const char* fname);
double getPhase(double wvl, double dy, double dz);
bool optimizeDelta(int* N, double* delta, double L, double wvl, double dz, double posy);

int copy_xray(  struct xray* in,   struct xray* out);
int copy_source(struct source* in, struct source* out);
int copy_crl(struct crl* in, struct crl* out);
int copy_detector(struct detector* in, struct detector* out);

int copy_field(struct field* in, struct field* out);
int show_field(struct field* argfield, char* message);
//int get_field(struct field* field, const char* fnameread);
