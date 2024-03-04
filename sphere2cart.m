function C = sphere2cart(S)
% Converts from Spherical to Cartesian
% S is an n-by-3 matrix of r,theta,&phi components of n points
% C is an n-by-3 matrix of x,y,&z components for those same points
% also works with a 3D array.
C = S(:,1,:).*[sin(S(:,2,:)).*cos(S(:,3,:)), ...
               sin(S(:,2,:)).*sin(S(:,3,:)), ...
               cos(S(:,2,:))                ];