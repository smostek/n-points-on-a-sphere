function [cs,T,m,R] = CS(C)
% Cross Sum - the average of the magnitudes of the tangential components of the forces
% Forces are calucluated using a normalized Coulomb's Law, with each vertex
% radiating an isotropic potential well of strength inversely proportional
% to distance
%
% C (3,:) - matrix of the cartesian coordinates of n-points
% THIS ASSUMES C IS COMPRISED OF UNIT VECTORS
%
% R(:,i) is the radial component of the force of the ith point
% T(:,i) is the non-radial (tangential) component of the force of the ith point
% m(i) is the magnitude of T(:,i)

D = C - permute(C,[1 3 2]);
F = sum(D./vecnorm(D).^3, 3, 'omitnan');
R = dot(F,C).*C;
T = F-R;
m = vecnorm(T);
cs= mean(m);