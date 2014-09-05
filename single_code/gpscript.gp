################################################ 
set terminal pngcairo size 1200,480;
set title 'Entrance plane (Re)';
set output 'entrance_plane.png';
stats "./entrance_plane.txt" nooutput;
N=STATS_records;
set grid;
set xrange[-N/2:N/2];
plot "./entrance_plane.txt" u (column(0)-N/2):1 w l;
################################################ 
################################################ 
reset;
set terminal pngcairo size 1200,480;
set title 'CRL plane (Re)';
set output 'crl_plane.png';
stats "./crl_plane.txt" nooutput;
N=STATS_records;
set grid;
set xrange[-N/2:N/2];
plot "./crl_plane.txt" u (column(0)-N/2):1 w l;
################################################ 
################################################ 
reset;
set terminal pngcairo size 1200,480;
set title 'Detector plane (Intensity)';
set output 'det_plane.png';
stats "./det_plane.txt" nooutput;
N=STATS_records;
set grid;
set xrange[-N/2:N/2];
plot "./det_plane.txt" u (column(0)) < N/2 ? (column(0)):(column(0)-N):3 w l;
################################################ 
