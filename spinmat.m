function P_AB = spinmat(i,alpha)
% A compression and improvement of rotx, roty, and rotz functions
arguments
    i
    % As an integer, the index (1, 2, or 3) of the axis about which the rotation
    % takes place
    %
    % As a 3-by-1 unit vector, the axis or rotation. So, for instance
    % spinmat(3,angle) = spinmat([0;0;1],angle)
    % 
    % Additional columns indicate the indices of the successive rotations, whereby
    % spinmat([i(:,1),i(:,2)], [alpha(1),alpha(2)])
    % is the same as
    % spinmat(i(:,1),alpha(1))*spinmat(i(:,2),alpha(2))

    alpha
    % Angle(s) of rotation, in radians when specified as a number, but may also be a
    % symbolic variable.
    %
    % Each column in alpha corresponds to the axis specified in the same column of i,
    % so each input must have the same number of columns.
    % 
    % Additional rows of alpha generate additional pages in the output matrix, where
    % each page corresponds to the 3-by-3 rotation matrix produced by the sequence of
    % rotations in same row of alpha
    % 
    % For example, if i = [3,2] and alpha is 2-by-2 double, then
    % P_AB(:,:,1) = spinmat(i, alpha(1,:)),  and
    % P_AB(:,:,2) = spinmat(i, alpha(2,:))
    % 
    % Unfortunately, for symbolic alpha, multiple rows are not supported
end

if isa(alpha,'symfun')
    alpha = formula(alpha);
end
if isa(i,'symfun')
    i = formula(i);
end
[is1,is2] = size(i);
[as1,as2] = size(alpha);
if as2~=is2
    error('i and alpha must have same number of columns')
end
if as1>1
    if isa(alpha,'sym')
        error('symbolic manipulation must occur within a single page')
    end
    a = reshape(alpha(:,1),[1 1 as1]);
else
    a = alpha(1);
end

if is1==3
    % rotation about an arbitrary axis:
    n = i(:,1);
    P_AB = n*n.' + cos(a).*(eye(3)-n*n.') + sin(a).*skew(n);
elseif is1==1
    switch i(1)
        case 1 %rotation about x-axis
            P_AB = [ones(1,1,as1),zeros(1,2,as1); zeros(2,1,as1),[cos(a),-sin(a);sin(a),cos(a)]];
        case 2 %rotation about y-axis
            P_AB = [cos(a),zeros(1,1,as1),sin(a);zeros(1,1,as1),ones(1,1,as1),zeros(1,1,as1);-sin(a),zeros(1,1,as1),cos(a)];
        case 3 %rotation about z_axis
            P_AB = [[cos(a),-sin(a);sin(a),cos(a)],zeros(2,1,as1);zeros(1,2,as1),ones(1,1,as1)];
    end
else
    error('i must either be a row vector or a 3-by-n matrix')
end

if as2>1
    if as1>1
        P_AB = pagemtimes( P_AB, spinmat(i(2:is2),alpha(:,2:as2)) );
    else
        P_AB = P_AB*spinmat(i(:,2:is2),alpha(:,2:as2));
    end
end

end% spinmat

function S = skew(v)
% turns the cross product into matrix multiplication:
S = [ 0,  -v(3), v(2);
     v(3),  0,  -v(1);
    -v(2), v(1),  0];
end