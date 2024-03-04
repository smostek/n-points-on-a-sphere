function m = mag(C,dir)
% magnitude
% each row of C is a cartesian vector,
% m gives the magnitudes of those vectors
arguments
    C    double
    dir  double= 2
end
m = sqrt(sum(C.^2,dir));