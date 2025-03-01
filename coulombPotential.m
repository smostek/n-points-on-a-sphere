function u = coulombPotential(C, ref)
% Returns the electrostatic potential at ref due to charges at C.
arguments
    C (3,:) {double}
    % cartesian coordinates of the point charges, 
    % IT IS ASSUMED C IS COMPRISED OF UNIT VECTORS

    ref = []
    % point at which the potential is to be determined, with three options:
    % 
    % 1) When unspecified, u will be the potential energy of the whole distribution
    % 
    % 2) When a scalar whole number, u will be the potential energy of just C(:,ref)
    % 
    % 3) When a 3-by-1 cartesian vector, u will be the potential energy at ref
end

if isempty(ref)
    u = sum(triu(1./permute(vecnorm(C-permute(C,[1,3,2])),[3 2 1]),1),'all');
    return
end

if isscalar(ref)
    i = ref;
    ref = C(:,i);
    C = C(:,[1:i-1 i+1:end]);
end
u = sum( vecnorm(C-ref).^(-1) );