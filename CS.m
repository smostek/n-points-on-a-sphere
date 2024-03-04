function [cs,T,m,R] = CS(C)
% Cross Sum, tangential force, magnitudes, raidal force
%
% C is an n-by-3 matrix of the cartesian coordinates of n-points
% THIS ASSUMES C IS COMPRISED OF UNIT VECTORS
% R(i) is the radial component of the force of the ith point
% T(i) is the non-radial (tangential) component of the force of the ith point

D = C - permute(C,[3 2 1]);
F = sum(D./mag(D).^3, 3, 'omitnan');
R = dot(F,C,2).*C;
T = F-R;
m = mag(T);
cs= sum(m);