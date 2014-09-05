//program genPlots to generate plots for crl propagation execution crlprop
//former bool genPlots(void)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(void){
char* plotfname="gpscript.gp";
char* fname;
char* data;
int i;
FILE* f = fopen(plotfname, "w");
for (i=0;i<3;i++){
fprintf(f,"################################################ \n");
if (i!=0) fprintf(f,"reset;\n");
fprintf(f,"set terminal pngcairo size 1200,480;\n");

if (i==0){
fname = "entrance_plane";
data = "(column(0)-N/2):1"; 
fprintf(f,"set title 'Entrance plane (Re)';\n");
}

else if (i==1){
fname = "crl_plane";
data = "(column(0)-N/2):1"; 
fprintf(f,"set title 'CRL plane (Re)';\n");
}

else if (i==2){
fname = "det_plane";
data = "(column(0)) < N/2 ? (column(0)):(column(0)-N):3"; 
fprintf(f,"set title 'Detector plane (Intensity)';\n");
}

fprintf(f,"set output '%s.png';\n", fname);
fprintf(f,"stats \"./%s.txt\" nooutput;\n", fname);
fprintf(f,"N=STATS_records;\n");
fprintf(f,"set grid;\n");
fprintf(f,"set xrange[-N/2:N/2];\n");
fprintf(f,"plot \"./%s.txt\" u %s w l;\n", fname, data);
fprintf(f,"################################################ \n");
}
if (f) fclose(f); f=NULL;
system("gnuplot gpscript.gp");
return 0;
}
