function cartV = sphvec2cartvec(C, sphV)
% Converts a vector field from cartesian to spherical
% C gives the cartesian coordinates of n points
% sphV gives the rhat, thetahat, and phihat components of the vector field
% at the coresponding row of P
x = C(:,1); y = C(:,2); z=C(:,3);
rho = sqrt(x.^2 + y.^2); rho(rho==0)=1e-16;
r   = sqrt(rho.^2+z.^2);

cartV = [sphV(:,1).*x./r + sphV(:,2).*x.*z./(r.*rho) - sphV(:,3).*y./rho, ...
         sphV(:,1).*y./r + sphV(:,2).*y.*z./(r.*rho) + sphV(:,3).*x./rho, ...
         sphV(:,1).*z./r - sphV(:,2).*rho./r];