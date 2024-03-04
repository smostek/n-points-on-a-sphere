function [sphV] = cartvec2sphvec(S, cartV)
%converts a vector field from cartesian to spherical

%S is an n-by-3 matrix containing the r, theta, and phi components of n
%points
%cartV is an n-by-3 matrix with each row containing the ihat, jhat and khat
%components of the vector at the position corresponding to that row of P
%sphV is an n-by-3 matrix with each row containing the rhat, thetahat, and
%phihat components of the same vectors.

st = sin(S(:,2)); sp = sin(S(:,3));
ct = cos(S(:,2)); cp = cos(S(:,3));

vx = cartV(:,1); vy = cartV(:,2); vz = cartV(:,3);

sphV = [vx.*st.*cp + vy.*st.*sp + vz.*ct, ...
        vx.*ct.*cp + vy.*ct.*sp - vz.*st, ...
       -vx.*sp     + vy.*cp];