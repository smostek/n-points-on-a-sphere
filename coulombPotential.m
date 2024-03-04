function u = coulombPotential(C, ref)
% returns the electrostatic potential at ref due to charges at C
% C contains row vectors giving the cartesian coordinates of point charges
% 
% if you want to know the potential at a point due to the charges at C,
% specify ref as the coordinates of that point.
%
% if you want to know the potential at the ith point of C due to all the
% others, specify ref as the index i.
%
% if you want to know the whole potential energy of the distribution, leave
% ref blank. THIS ASSUMES C IS COMPRISED OF UNIT VECTORS
if nargin==1
    [~,~,~,R]=CS(C);
    u = sum(mag(R));
    return
end

if isscalar(ref)
    i = ref;
    ref = C(i,:);
    C = C([1:i-1 i+1:end],:);
end
u = sum( mag(C-ref).^(-1) );