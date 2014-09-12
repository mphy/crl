set title 'Detector plane (Intensity)';
stats "./det_plane.txt" nooutput;
N=STATS_records;
set grid;
#set xrange[-N/2:N/2];
#set yrange[0:1000];
plot "./det_plane.txt" u 1:4 ;						#read coord
#.............................................................................. 
