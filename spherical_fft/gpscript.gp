# 1:spatialcoord   2:Re   3:Im   4:intesity(amplitude)/phase


#.............................................................................. 
set terminal pngcairo size 1200,480;
set title 'Entrance plane (Re)';
set output 'entrance_plane.png';
stats "./entrance_plane.txt" nooutput;
N=STATS_records;
L=1e-3;
delta=L/N;
set grid;
#set xrange[-N/2:N/2];
plot "./entrance_plane.txt" u 1:2 w l;				#read coord
#plot "./entrance_plane.txt" u ((column(0)-N/2)*delta):2 w l;	#compute coord (defined)
#plot "./entrance_plane.txt" u (column(0)-N/2):2 w l;		#no coord


#.............................................................................. 
#.............................................................................. 
reset;
set terminal pngcairo size 1200,480;
set title 'CRL plane (Re)';
set output 'crl_plane.png';
stats "./crl_plane1.txt" nooutput;
N=STATS_records;
L=1e-3;
delta=L/N;
set grid;
#set xrange[-N/2:N/2];
plot "./crl_plane1.txt" u 1:2 w l;				#read coord
#plot "./crl_plane1.txt" u ((column(0)-N/2)*delta):2 w l;	#compute coord (defined)
#plot "./crl_plane1.txt" u (column(0)-N/2):2 w l;		#no coord

		
#.............................................................................. 
#.............................................................................. 
reset;
set terminal pngcairo size 1200,480;
set title 'Detector plane (Intensity)';
set output 'det_plane.png';
stats "./det_plane.txt" nooutput;
N=STATS_records;
set grid;
#set xrange[-N/2:N/2];
#set yrange[0:1000];
plot "./det_plane.txt" u 1:4 w l;						#read coord


#plot "./det_plane.txt" u (column(0)) < N/2 ? (column(0)):(column(0)-N):4 w l; 	#no coord + fftshift

#L=100e-6;
#delta=L/N;
#plot "./det_plane.txt" u ((column(0)-N/2)*delta):4 w l;				#compute coord (defined)
#.............................................................................. 
