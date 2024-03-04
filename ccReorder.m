function order = ccReorder(C)
% reorders a set of points on a sphere so that they cycle counter-clockwise
% when viewed from outside the sphere.
r=size(C,1);
if r>2
    N = cross(C(2,:)-C(1,:), C(3,:)-C(1,:));
    if sum((N-C(1,:)).^2) > sum((N+C(1,:)).^2)
        N = -N;
    end
    v = cartvec2sphvec(repmat(cart2sphere(N), [r,1]), C);
    [~,order] = sort(atan2(-v(:,2), v(:,3)));
else
    order=1:r;
end