function diss = sampleSurface(flow, filt)
% This function takes a flowrate and filtration time and returns a
% dissolution time
% The surface is fitted to experimental data

a = 0.3403;
b = 13.92;
c = 2.41;
d = 3.056;


diss = filt.^a .* (b*exp(-c .* flow) + d);