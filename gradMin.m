function [varvec, history] = gradMin(f, varvec)
%use gradient descent to find a minimum in arbitrary dimensions
%f takes in a column vector and returns a scalar
%varvec is a column vector of the initial positions
nDim = length(varvec);
history = varvec;
oldy = inf;
newy = f(varvec); dx=0.1;
while all([newy>1e-4 dx>1e-4 oldy>newy])
    dvec = getDir(dx);
    while sum(dvec)==0
        dx = dx*0.5; dvec = getDir(dx);
    end
    varvec = varvec+dvec;
    if nargout==2
        history = [history varvec];
    end
    oldy = newy;
    newy = round(f(varvec), floor(-log10(0.1*oldy))+3);
end

    function dvar = getDir(step)
    delta = zeros(nDim,1);
    y = f(varvec); A = [zeros(nDim+1,nDim), [varvec; y]];
    noChangeDims = false(nDim,1); p = round(abs(log10(step))*2);
    for d = 1:nDim
        dimys = [0 round(y,p) 0];
        for j = [-1,1]
            delta(d) = step*j;
            dimys(j+2)=round(f(varvec+delta), p);
        end
        [low, mini] = min(dimys);
        mini = mini-2+(dimys(1)==dimys(2)&&dimys(3)>=dimys(2));
        delta(d)=mini*step;
        A(:,d) = [varvec+delta; low];
        noChangeDims(d)= mini==0;
        delta(d)=0;
    end
    A([noChangeDims;false],:)=[]; A(:,[noChangeDims;false])=[];
    c = A'\ones(length(A),1);
    dvar = zeros(nDim,1);
    c = c(1:end-1)/c(end);
    c = c/sqrt(sum(c.^2));
    dvar(~noChangeDims)=step*c;
    end%GetDir

end