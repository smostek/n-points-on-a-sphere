function S = cart2sphere(C)
%Converts Cartesian to Spherical
%C is an n-by-3 matrix of the x,y,&z components of n points
%S is an n-by-3 matrix of the r,theta,&phi components of those same points
% also works with 3D matrices
r = sqrt(sum(C.^2,2));
S = [r, acos(C(:,3,:)./r), mod(atan2(C(:,2,:),C(:,1,:)),2*pi)];