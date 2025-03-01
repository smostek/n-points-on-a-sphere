function [x, history] = gradMin(f, x, tolerance, del)
% Use gradient descent to find a local minimum in arbitrary dimensions
arguments
    f function_handle 
    % Function to be minimized: takes in a row vector and returns a scalar

    x (1,:) 
    % Row vector of the initial positions; 
    % as an output, the coordinates of the local minimum

    tolerance = 1e-4  
    % Value of f at which minimization halts

    del = 0.1 
    % Initial step size

    % if two outputs are specified, the values of f at each step of the
    % minimization will be returned in column vector history
end

nDim = length(x);
history = x;
oldy = inf;
newy = f(x);
while all([newy>tolerance, del>tolerance, oldy>newy]) %[not minimized, taking meaninful steps, decreasing]
    dx = getDir(del);
    while sum(dx)==0 %dvec will be all zeros if stepping in any direction goes uphill
        del = del*0.5;
        dx = getDir(del);
    end
    x = x+dx;
    if nargout==2
        history = [history; x];
    end
    oldy = newy;
    newy = round(f(x), floor(-log10(oldy))+4);
end

    function dvar = getDir(step)
    delta = zeros(1,nDim);
    y = f(x);
    A = [zeros(nDim,nDim+1); [x y]];
    noChangeDims = false(1,nDim);
    p = floor(-log10(step))+4; %precision
    for d = 1:nDim
        dimys = [0 round(y,p) 0];
        for j = [-1,1]
            delta(d) = step*j;
            dimys(j+2)=round(f(x+delta), p);
        end
        [low, mini] = min(dimys);
        mini = mini-2+(dimys(1)==dimys(2)&&dimys(3)>=dimys(2));
        delta(d)=mini*step;
        A(d,:) = [x+delta, low];
        noChangeDims(d)= mini==0;
        delta(d)=0;
    end
    A([noChangeDims,false],:)=[]; A(:,[noChangeDims,false])=[];
    v = A\ones(length(A),1);
    v = v(1:end-1)/vecnorm(v(1:end-1));
    dvar = zeros(1,nDim);
    dvar(~noChangeDims)=step*v;
    end%GetDir

end