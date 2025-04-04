classdef shape < handle
properties
    cartCoor (3,:)
    % Cartesian coordinates of vertices

    minOrd (1,1)   
    % -log10(CS(self.cartCoor), the higher the value, the better the minimization (see CS)

    dual (3,:)
    % Cartesian coordinates of circumcenters of faces (see getSymmetry)

    faceOrder (:,1) cell 
    % Each cell corresponds to a face of the polyhedron, giving the indices of the
    % vertices that face neighbors (see getSymmetry)

    symmetry (5,1)
    % vector describing the symmetry group of the polyhedron and its orientation
    % (see getSymmetry)

    ratings 
    % vector giving the value of the distribution in each rating scheme (see rate)

end %properties

methods
function [s,coorf,varvec,nums] = shape(name,minimize,Opts)
% A polyhedron with a unit circumsphere
arguments 
    name = ''
    % char vector of numbers - each number creates a polygonal 'band' with that
    % many vertices. A band can be specified as being rotated by adding a tick (`)
    % after that number. For example, an Icosahedron can be constructed by passing 
    % any of '1 5 5` 1' , '3 3` 3 3`' , or '2 2` 4` 2` 2'
    % 
    % Alternatively, one can specify name as just a 3-by-n matrix representing the
    % cartesian coordinates of the vertices

    minimize logical = true 
    % When true, the heights of the bands will be adjusted using gradient descent
    % (see function gradMin) to miniize the forces on each vertex (see function CS)

    Opts.random logical = false
    % When true, the coordinates of the shape will be randomized. In this case, 
    % name must be specified as num2str(n) with n being the desired number of vertices
end

if ~isempty(name)
    if ischar(name)
        %get numbers from name
        if Opts.random
            n=str2double(name);
            if isnan(n)
                error('randomized distributions take only a single number name')
            end
            s.cartCoor=coorConvert([ones(1,n); asin(2*rand(1,n)-1); pi*(2*rand(1,n)-1)],'sph');
            return
        end
        [nums,bands,primes,NPole,SPole]=nameInterpreter(name);
        B = length(bands);
        n = sum(bands)+NPole+SPole;
        
        % assemble the shape
        phi = mod([zeros(1,NPole), cell2mat(arrayfun(@(b)(0:bands(b)-1)*2*pi/bands(b)+ ...
            primes(b)*pi/bands(b),1:B,'un',0)), zeros(1,SPole)],2*pi);
        vari= repelem(1:B,bands)';
    
        %minimize the setup:
        caps={[],-pi/2,[],pi/2};
        coorf=@(vars) coorConvert([ones(1,n); [caps{NPole+1},vars(vari)',caps{SPole+3}]; phi],'sph');
        varvec = asin((1/B)*(2*(1:B).'-(B+1)));
        if minimize
            varvec = gradMin(@(v)CS(coorf(v)),varvec);
            s.minOrd = -log10(CS(coorf(varvec)));
        end
        s.cartCoor = coorf(varvec);
    elseif isnumeric(name) && size(name,1)==3
        s.cartCoor = name;
    
    end
end
    
    function [nums,bands,primes,NPole,SPole] = nameInterpreter(txt)
        nums=textscan(txt, '%f', 'Delimiter', {' ', '` '});
        if length(nums)>1
            error('name must only contain integers, spaces, and ticks')
        end
        nums=nums{1}';
        
        NPole= nums(1)==1;
        SPole= nums(end)==1;
        bands= nums((1+NPole):(end-SPole));
        if any(nums<1) || any(nums~=round(nums))
            error('name must use only integers greater than zero')
        end
        if any(bands==1)
            error('cannot have 1s in middle of name vector')
        end
        if NPole&&txt(2)=='`' || SPole&&txt(end)=='`'
            error('cannot apply rotation to poles')
        end
        
        % get the ticks from name
        primes=false(1,length(bands));
        if NPole; txt([1,2])=    []; end
        if SPole; txt(end-1:end)=[]; end
        spaces= find(txt==' ');
        ticks = find(txt=='`');
        if max(ticks)>max(spaces)
            primes(end)=true;
            ticks(end)=[];
        end
        for t=1:length(ticks)
            primes(spaces== ticks(t)+1)= true;
        end
    end

end
%%%%%%%%%%%%%

function [tdata,rdata] = forceMin(self,mu,Opts)
% Minimize the Cross Sum of a distribution using a modified Newton's Method
arguments
    self shape

    mu double = 8
    % Target minimization order ( =-log10(CS) ), higher values correspond to greater
    % precision

    Opts.dispVal logical = true 
    % When true, mu will be displayed in the command window at the end of the process
end
N = length(self);
if N>1
    fprintf('Minimizing Distributions:')
    for a = 1:N
        self(a).forceMin(mu,'disp',false)
        if round( floor(10*a/N)-floor(10*((a-1)/N)) )
            fprintf(' %%')
        end
    end
    fprintf(' Done!\n')
    return
end

C = self.cartCoor;
f = @(c,p) p*round(CS(c)/p);

dt = 0.1;
[cs,T,M,R] = CS(C);
if nargout
    tdata = M;
    rdata = vecnorm(R);
end
nxt = @(c,t,f,m) cos(t*m).*c + sin(t*m).*f./m; %take a step along the sphere
ohno = false; %flag for if minimization is failing
tol = 10^(-mu);
while cs > tol && ~ohno
    p = min([0.1, cs*1e-3]); %precision increases as cs decreases
    cs = f(C,p);
    %shrink step-size until it's clear whether to go forwards or backwards:
    dir = [f(nxt(C,-dt,T,M),p), cs, f(nxt(C,dt,T,M),p)];
    [~,fb] = min(dir); fb=fb-2;
    while ( fb == 0 || any(diff(dir)==0) ) && ~ohno
        dt = dt*0.5;
        ohno = dt<tol;
        dir = [f(nxt(C,-dt,T,M),p), cs, f(nxt(C,dt,T,M),p)];
        [~,fb] = min(dir); fb = fb-2;
    end

    if ohno
        break
    end

    % take as many steps as keeps the score decreasing
    t = 0;
    steps = 0;
    while f(nxt(C,t+fb*dt,T,M),p) < f(nxt(C,t,T,M),p)
        t=t+fb*dt;
        steps = steps+1;
    end
    if steps>10; dt = dt*3; end
    C=nxt(C,t,T,M);
    [cs,T,M,R] = CS(C);
    if nargout
        tdata = [tdata;M];
        rdata = [rdata;vecnorm(R)];
    end
end
self.cartCoor = C./vecnorm(C);
self.minOrd = -log10(cs);
if Opts.dispVal
    fprintf('Minimized to mu = %.2f\n',self.minOrd)
end
% reevaluate any attributes that might have changed with the coordinates:
if self.symmetry(1)~=0
    self.getSymmetry
end
if ~isempty(self.ratings)
    self.rate
end

end%forceMin

function see(self,style, Opts)
arguments
    self (1,1) shape

    style char = 'normal'
    % Use to specify the visualization type; one of the following:
    %   normal      - Plots the points on a sphere
    %   solid       - Shows the polyhedron formed by making the points vertices
    %   force       - Shows the polyherdon with electrostatic force vectors on each
    %                 vertex
    %   flat        - South-pole based stereographic projection
    %   contourPlot - Indicates strength of the electrostatic potential on the 
    %                 surface of the sphere
    %   symmetry    - Highlight the symmetry group of the shape

    Opts.axes matlab.graphics.axis.Axes = gca
    % axes object on which the visuals will be displayed

    Opts.coor = self.cartCoor
    % Cartesian coordinates of vertices of the shape - use if you  want to visualize
    % coordinates without defining a new shape

    Opts.clear logical = true
    % When true, clears the target axis before drawing visuals

    Opts.mainColor = 'k'
    % Color character or RGB vector determining the color of shape's vertices and edges
end
validNames={'normal','solid','force','flat','contourPlot','symmetry'};
style = validatestring(style, validNames);

if isempty(gcf().Tag) % figure's tag can override the check to whether the axis a valie target
    if isempty(gca().Children) % newly created figure
        set(gcf,'Position',[80+10*rand 570+10*rand 560 420],'Color','w','UserData','standard')
    elseif ~strcmp(gcf().UserData,'standard') % existing special figure - create new one
        figure('Position',[80+10*rand 570+10*rand 560 420],'Color','w','UserData','standard');
        Opts.axes = axes;
    end
end
trgAx = Opts.axes;
view(trgAx, [0.05 -1 .2])
axis(trgAx, 'equal','off');
if Opts.clear; cla(trgAx); end
hold(trgAx, 'on')
feval(eval(['@' style]))
hold(trgAx, 'off')

function normal
    % Plots the points on a sphere
    [x,y,z]=sphere(15);
    surf(trgAx,x,y,z, 'EdgeColor',[0 0.9 0.9], 'FaceAlpha',0.2, 'FaceColor','c')
    C = Opts.coor;
    scatter3(trgAx, C(1,:), C(2,:), C(3,:), 50, Opts.mainColor, ...
        'filled','HitTest','off')
end
function solid
    % Shows the polyhedron defined by treating the electrons as vertices
    if isempty(self.faceOrder) || ~isequal(self.cartCoor,Opts.coor)
        self.getSymmetry('coor',Opts.coor,'dispSym',false);
    end
    f = @(a,i) Opts.coor(i,self.faceOrder{a});
    arrayfun(@(a) patch(trgAx,f(a,1),f(a,2),f(a,3),[0.9,0.9,0.9], ...
        'EdgeColor',Opts.mainColor), 1:size(self.dual,2));
end
function force
    % Shows the polyherdon with electrostatic force vectors on each vertex
    C = Opts.coor;
    [~,T] = CS(C);
    if isempty(self.dual) || ~isequal(self.cartCoor,C)
        self.getSymmetry('coor',C,'dispSym',false)
    end
    f = @(a,i) C(i,self.faceOrder{a});
    arrayfun(@(a) patch(trgAx,f(a,1),f(a,2),f(a,3),'w','EdgeColor',Opts.mainColor, ...
        'FaceAlpha',0.6, 'LineWidth',2), 1:size(self.dual,2));
    quiver3(C(1,:),C(2,:),C(3,:),T(1,:),T(2,:),T(3,:))
end
function flat
    % South-pole based stereographic projection
    if isempty(self.faceOrder) || ~isequal(self.cartCoor,Opts.coor)
        self.getSymmetry('coor',Opts.coor,'dispSym',false);
    end
    syms = self.symmetry;
    S=coorConvert(spinmat([3,2,3],[pi/2,0,0]-syms([5,4,3])')*Opts.coor,'cart');
    spi = find(pi/2-S(2,:)<0.2,1); %south pole index
    rho=2*tan((S(2,:)+pi/2)/2);
    if ~isempty(spi)
        % south pole goes to infinity in this projection, so it needs to be
        % handled differently.
        rho(spi)=NaN;
        x=rho.*cos(S(3,:));
        y=rho.*sin(S(3,:));
        % draw circle:
        t=linspace(0,2*pi,40);
        R=1.15*max(rho);
        plot(trgAx,R*cos(t),R*sin(t),'k')
        % find and draw the edges that connect to the south spole:
        e = cell2mat(cellfun(@(c) c(mod(find(c==spi,1),length(c))+1),self.faceOrder,'un',0)');
        arrayfun(@(a) plot(trgAx,x(e(a))*[1, R/rho(e(a))],y(e(a))*[1,R/rho(e(a))],'k'),1:length(e))
    else
        x=rho.*cos(S(3,:));
        y=rho.*sin(S(3,:));
    end
    scatter(trgAx,x,y,50,Opts.mainColor,'filled')
    arrayfun(@(a) patch(trgAx,x(self.faceOrder{a}),y(self.faceOrder{a}),'w', ...
        'EdgeColor',Opts.mainColor,'FaceAlpha',0), 1:size(self.dual,2));
    view(trgAx,[0 0 1])
end
function contourPlot
    % Indicates strength of the electrostatic potential on the surface of
    % the sphere
    res=100;
    th=asin((1/res)*(2*(1:res)-(res+1)));
    ph=pi*linspace(-1,1,res);
    [theta,phi] = meshgrid(th,ph);
    U = arrayfun(@(t,p) coulombPotential(Opts.coor,coorConvert([1;t;p],'sph')),theta,phi);
    [~,edges]=histcounts(U(U<1.5*size(Opts.coor,2)));
    M = contourc(th,ph,U,edges);
    i=1;
    while i<length(M)
        Ni=M(2,i);
        C=coorConvert([ones(1,Ni); M(:,i+1:i+Ni)],'sph');
        plot3(trgAx,C(1,:),C(2,:),C(3,:),'k')
        i = i+Ni+1;
    end
    [x,y,z]=sphere(res);
    surf(trgAx,x,y,z,'FaceColor','w','EdgeColor','none')
end
function symmetry
    % Highlight the symmetry group of the shape
    C = Opts.coor;
    if self.symmetry(1)==0 || ~isequal(self.cartCoor,C)
        self.getSymmetry('coor',Opts.coor,'dispSym',false);
    end

    syms = num2cell(self.symmetry);
    [type,k,psi,the,phi]=deal(syms{:});
    linCoor = {};
    linCol = []; %line colors
    linThk = []; %line widths
    surfCoor = {};
    surfCol = [];
    blu = [0.00,0.45,0.74];
    yel = [0.93,0.69,0.13]; 
    red = [0.64,0.08,0.18]; 
    % Draw the solid:
    [x,y,z]=sphere(40);
    surf(trgAx,0.99*x,0.99*y,0.99*z,'EdgeColor','none','FaceAlpha',0.7,'FaceColor',[.9,.9,.9])
    linstyl = ':';
    for a = 1:size(self.dual,2)
        line = sphereLine(C(:,self.faceOrder{a}([1:end,1])));
        plot3(trgAx,line(1,:),line(2,:),line(3,:),'Color',Opts.mainColor,'LineStyle',linstyl,'LineWidth',1)
    end
    scatter3(trgAx,C(1,:),C(2,:),C(3,:),30,Opts.mainColor,'filled')

    %take care of a couple special cases:
    if size(C,2)==1
        return
    end
    if size(C,2)==2
        C = C(:,1)+(C(:,2)-C(:,1))*[-0.5,1.5];
        plot3(trgAx,C(1,:),C(2,:),C(3,:),'Color',red)
        return
    end

    switch type
        case 1
            % C_k:
            linCoor = arrayfun(@(a)sphereLine(coorConvert([1,1;-1.57,1.57;[a,a]*2*pi/k],'sph')),1:k,'un',0);
            linCol = ones(length(linCoor),1)*blu;
            linThk = ones(length(linCoor),1)*1.5;

            a = reshape([0;0.2;0]+2*pi*(1:k)/k,1,[]);
            H = repmat([0,0,0.1],1,k);
            surfCoor = mat2cell([cos(a);sin(a);H],3,3*ones(1,k));
            surfCol = ones(k,1)*blu;
        case 2
            % C_kv:
            linCoor = arrayfun(@(a)sphereLine(coorConvert([1,1;-1.57,1.57;[a,a]*pi/k],'sph')),1:2*k,'un',0);
            linCol = repmat([blu;yel],k,1);
            linThk = ones(length(linCoor),1)*1.5;

            a = reshape([-0.1;0.1;0]+pi*(1:2*k)/k,1,[]);
            H = repmat(0.11*[-0.5,-0.5,1,0.5,0.5,-1],1,k);
            surfCoor = mat2cell([cos(a);sin(a);H],3,3*ones(1,2*k));
            surfCol = repmat([blu;yel],k,1);
        case 3
            % C_kh:
            t = linspace(0,2*pi,40);
            linCoor = arrayfun(@(a)sphereLine(coorConvert([1,1;-1.57,1.57;[a,a]*2*pi/k],'sph')),1:k,'un',0);
            linCoor{end+1} = [cos(t);sin(t);zeros(size(t))];
            linCol = [ones(length(linCoor)-1,1)*blu;red];
            linThk = ones(length(linCoor),1)*1.5;

            p = reshape([0;-0.2;0;0;-0.2;0]+2*pi*(1:k)/k,1,[]);
            p = [p, reshape([0;0.2;0]+pi*(1:2:2*k)/k,1,[])];
            t = repmat([-pi/4-[0,0,0.1],+pi/4+[0,0,0.1]],1,k);
            t = [t, repmat([0.1,0,-0.1],1,k)];
            surfCoor = mat2cell(coorConvert([ones(size(t));t;p],'sph'),3,3*ones(1,3*k));
            surfCol = [ones(2*k,1)*blu; ones(k,1)*red];
        case 4
            % R_k:
            linCoor = arrayfun(@(a)sphereLine(coorConvert([1,1;-1.57,1.57;[a,a]*pi/k],'sph')),1:2*k,'un',0);
            linCol = ones(2*k,1)*blu;
            linThk = ones(2*k,1)*1.5;

            a = reshape([0;-0.2;0]+pi*(1:2*k)/k,1,[]);
            H(3:3:6*k) = 0.1*(-1).^(0:2*k-1);
            surfCoor = mat2cell([cos(a);sin(a);H],3,3*ones(1,2*k));
            surfCol = ones(2*k,1)*blu;
        case 5
            % D_k:
            linCoor = arrayfun(@(a)sphereLine(coorConvert([1,1;-1.57,1.57;[a,a]*pi/k],'sph')),1:2*k,'un',0);
            linCol = repmat([blu;yel],k,1);
            linThk = ones(2*k,1)*1.5;

            a = reshape(0.15*[0;1;0;-1]+pi*(1:2*k)/k,1,[]);
            H = repmat(0.075*[1,-1,-1,1],1,2*k);
            surfCoor = mat2cell([cos(a);sin(a);H],3,4*ones(1,2*k));
            surfCol = repmat([blu;yel],k,1);
        case 6
            % D_kd:
            linCoor = arrayfun(@(a)sphereLine(coorConvert([1,1;-pi/2,0;[a,a]*pi/k],'sph')),1:2*k,'un',0);
            linCoor = [linCoor, arrayfun(@(a)sphereLine(coorConvert([1,1;pi/2,0;0,a*pi/k],'sph')),1:2*k,'un',0)];
            t = linspace(0,2*pi,40);
            linCoor{end+1} = [cos(t);sin(t);zeros(size(t))];
            linCol = [repmat([blu;yel],k,1);repmat([yel;blu],k,1);red];
            linThk = ones(length(linCoor),1)*1.5;

            p = reshape(0.15*[0;-1;1;0;-1;1]+pi*(1:2*k)/k,1,[]);
            t = repmat([-pi/4-0.1*[-1,1,1],+pi/4+0.1*[-1,1,1]],1,2*k);
            p = [p, reshape(0.075*[1;-1;-1;1]+pi*(1:2:4*k)/(2*k),1,[])];
            t = [t, repmat(0.15*[0,1,0,-1,0,-1,0,1],1,k)];
            surfCoor = mat2cell(coorConvert([ones(size(t));t;p],'sph'),3,[3*ones(1,4*k),4*ones(1,2*k)]);
            surfCol = [repmat([blu;yel;yel;blu],k,1); ones(2*k,1)*red];
        case 7
            % D_kh:
            linCoor = arrayfun(@(a)sphereLine(coorConvert([1,1;-1.57,1.57;[a,a]*pi/k],'sph')),1:2*k,'un',0);
            t = linspace(0,2*pi,40);
            linCoor{end+1} = [cos(t);sin(t);zeros(size(t))];
            linCol = [repmat([blu;yel],k,1);red];
            linThk = ones(length(linCoor),1)*1.5;

            p = reshape(0.15*[0;-1;1;0;-1;1]+pi*(1:2*k)/k,1,[]); %vert mirror triangles
            t = repmat([-pi/4-0.1*[-1,1,1],+pi/4+0.1*[-1,1,1]],1,2*k);
            p = [p, reshape([0;0.2;0;-0.2]+pi*(1:2:2*k)/k,1,[])]; %hori mirror triangles
            t = [t, repmat([0.1,0,-0.1,0],1,k)];
            surfCoor = mat2cell(coorConvert([ones(size(t));t;p],'sph'),3,[3*ones(1,4*k),4*ones(1,k)]);
            surfCol = [repmat([blu;blu;yel;yel],k,1); ones(2*k,1)*red];
        case {8,9,10}
            % T:
            % coordinates of a tetrahedron:
            C = coorConvert([ones(1,4);-pi/2,0.3398*ones(1,3);0,pi*([0,2,4])/3],'sph');
            % path between vertices that minimizes backtracking:
            linCoor = {sphereLine(C(:,[1 2 3 1 4 3 4 2]))};
            linCoor{2} = -linCoor{1}; %dual
            linCol = [blu;yel];
            if type==9
                % T_h:
                C = coorConvert([ones(1,6);0.6155*(-1).^(1:6);pi*(0:5)/3],'sph');
                linCoor(3:5)=arrayfun(@(a)sphereLine(C(:,[a a+1]),2*pi),1:2:5,'un',0);
                linCol = [linCol;red;red;red];
                linThk = [1;1;1.5;1.5;1.5];
            elseif type==10
                % T_d:
                C = [C,-C];
                linCoor{3} = sphereLine(C(:,[1 7 4 5 3 8 1 6 3 5 2 7 4 6 1 8 2]));
                linCol = [linCol; red];
                linThk = [1;1;1.5];
            else
                linThk = [1.5;1.5];
            end
        case {11,12}
            % O:
            t = linspace(0,2*pi);
            o=zeros(size(t));
            % draw an octahedron
            linCoor = {[cos(t);sin(t);o], [sin(t);o;cos(t)], [o;cos(t);sin(t)]};
            % coordinates of a cube:
            C = coorConvert([ones(1,8); asin(1/sqrt(3))*(-1).^(1:8); pi*([1,1,3,3,5,5,7,7])/4],'sph');
            % path between cube corners that minimizes backtracking
            linCoor{4} = sphereLine(C(:,[1 3 5 7 1 2 4 6 8 2 1 3 4 6 5 7 8]));
            linCol = [ones(3,1)*blu;yel];
            if type==12
                % O_h:
                % coordinates of rhombic dodecahedron:
                C = [C,[0 0 0 0 1 -1; 0 0 1 -1 0 0; 1 -1 0 0 0 0]];
                T = {[1 11 3 9 5 12 6 10 2 13 1 9 7 12 8 10 4 14 5],[3 14 6],[8 13 7],[2 11 4]};
                linCoor(5:8) = arrayfun(@(a)sphereLine(C(:,T{a})),1:4,'un',0);
                linCol = [linCol;ones(4,1)*red];
                linThk = [ones(4,1); ones(4,1)*1.5];
            else
                linThk = ones(4,1)*1.5;
            end
        case {13,14}
            % I:
            % coordinates of an icosahedron and a dodecahedron:
            I = coorConvert([ones(1,12); -pi/2,-0.4637*ones(1,5),+0.4636*ones(1,5),pi/2; 0,pi*(1:2:9)/5,pi*(0:2:8)/5,0],'sph');
            D = coorConvert([ones(1,20);0.9184*(-1).^(1:10),0.1887*(-1).^(1:10);pi*[0:9,0:9]/5],'sph');
            %draw the arcs in an order that avoids retracing paths:
            T = {[1 5 6 1 2 3 1 4 5 11 6 2 8 3 4 10 5],[12 10 9 12 8 7 12 11 7 2],[4 9 8],[3 9],[10 11],[6,7]};
            linCoor = arrayfun(@(a)sphereLine(I(:,T{a})),1:6,'un',0);
            T = {[1:2:9,1,11:20,10,2:2:10],[2,12],[3,13],[4,14],[5,15],[6,16],[7,17],[8,18],[9,19],[11,20]};
            linCoor(7:16) = arrayfun(@(a)sphereLine(D(:,T{a})),1:10,'un',0);
            linCol = [ones(6,1)*blu; ones(10,1)*yel];
            if type ==14
                % I_h:
                C = [I,D];
                % path between vertices of the rhombic triacontahedron:
                T = {[06 13 02 15 03 17 04 19 05 21 06 23 02 25 03 27 04 29 05 31 06],...
                     [07 24 08 26 09 28 10 30 11 32 07 14 08 16 09 18 10 20 11 22 07],...
                     [01 13], [01 15], [01 17], [01 19], [01 21],...
                     [12 14], [12 16], [12 18], [12 20], [12 22],...
                     [07 23], [02 24], [08 25], [03 26], [09 27],...
                     [04 28], [10 29], [05 30], [11 31], [06 32]};
                linCoor(17:38) = arrayfun(@(a)sphereLine(C(:,T{a})),1:22,'un',0);
                linCol = [linCol; ones(22,1)*red];
                linThk = [ones(16,1); ones(22,1)*1.5];
            else
                linThk = ones(16,1)*1.5;
            end
    end

    
    %Transformation matrix:
    P_AB = spinmat([3,2,3],[psi,the,phi]);
    for l = 1:length(linCoor)
        v = P_AB*linCoor{l};
        plot3(trgAx,v(1,:),v(2,:),v(3,:),'Color',linCol(l,:),'LineWidth',linThk(l))
    end
    for s = 1:length(surfCoor)
        v = P_AB*surfCoor{s};
        patch(trgAx,v(1,:),v(2,:),v(3,:),surfCol(s,:),'EdgeColor','none')
    end

        function line = sphereLine(C,t_f)
            % Gives the cartesian coordinates of a geodesic connectig points on a sphere
            % C gives the cartesian coordinates of each of the points being connected
            % t_f (optional) gives the arc length of each of the segments;
            % useful, for instance, if you want to draw the full great
            % circles containing each pair of points
            line = [];
            if nargin==1
                ang = @(j) acos(C(:,j).'*C(:,j+1));
            else
                ang = @(~) t_f;
            end
            for i = 1:size(C,2)-1
                proj = (eye(3)-C(:,i)*C(:,i).')*C(:,i+1);
                proj = proj/vecnorm(proj);
                t_f = ang(i);
                q = linspace(0,t_f,2*ceil(5*t_f));
                line = [line, cos(q).*C(:,i)+sin(q).*proj ];
            end
        end
end

end%see

function getSymmetry(self,tol,Opts)
% Determines the symmetry group of a distribution of points on the unit
% sphere given their coordinates
arguments
    self shape

    tol = 0.01
    % Determines, in a variety of checks, how close two elements should be to 
    % qualify as the same

    Opts.coor = self.cartCoor
    % Cartesian coordinates of shape whose symmetry is to be determined

    Opts.dispSym = 'true'
    % When true, displays the identified symmetry in the command window
end
% This function generates a 5-by-1 vector, syms, of the following structure:
% syms(1) is an integer from 1 to 14 specifying the type of symmetry in the
% order C_k, C_kv, C_kh, R_k*, D_k, D_kd, D_kh,
%       T, T_h, T_d, O, O_h, I, I_h
%                           *R_k is my notation for the rotoreflection
%                           symmetry S_2k
% syms(2) gives the k value for the first 7 of the groups above. If the
% shape has polyhedral symmetry, this value is arbitrary.
%
% syms(3:5) are psi, theta, and phi and are used to orient the principle axis
N = length(self);
if N>1
    fprintf('Evaluating Symmetry Groups:')
    for a=1:N
        try
        self(a).getSymmetry(tol,'disp',false)
        catch
            error(['error in shape #' num2str(a)])
        end
        if round( floor(10*a/N)-floor(10*((a-1)/N)) )
            fprintf(' %%')
        end
    end
    fprintf(' Done!\n')
    return
end

C=Opts.coor;
n = size(C,2);
symcheck = @(M,trans) all(ismember(tol*round(M'/tol),tol*round((trans*M).'/tol),'rows'));

%% Check if shape is a ploygon:
if n==1
    self.symmetry = [1;inf;atan2(C(2),C(1));acos(C(3));0];
    if Opts.dispSym
        fprintf('Symmetry Group: C_inf\n')
    end
    return
end
if n==2
    self.symmetry = [5;inf;atan2(C(2),C(1));acos(C(3));0];
    if Opts.dispSym
        fprintf('Symmetry Group: D_inf\n')
    end
    return
end
if iscoplanar(C, tol)
    D = cross(C(:,2)-C(:,1), C(:,3)-C(:,1));
    D = D/vecnorm(D);
    syms = zeros(5,1);
    syms([3 4]) = [atan2(D(2),D(1)); acos(D(3))]; % [psi; theta]
    C = spinmat([2 3],-syms([4,3])')*C;
    [~,ord] = sort(atan2(C(2,:),C(1,:)));
    C = C(:,ord);

    self.dual = D;
    self.faceOrder = {ord};

    % rotational symmetry:
    found = false;
    k = n;
    while k~=1 && ~found
        if symcheck(C,spinmat(3,2*pi/k))
            found = true;
        else
            fact=factor(k);
            k=k/fact(1);
        end
    end
    syms(2)=k;
    phi=0;
    % vertical mirror symmetry:
    if n>k && n>3
        C = [C, 0.5*(C+C(:,[2:end,1]))]-mean(C,2); % add midpoints
        bandphi = sort(atan2(C(2,:),C(1,:)));
        % iterate through a pi/k wide slice looking for a mirror plane
        j = 1;
        s_v=false;
        while bandphi(j) < (pi/k)*(1-k)+tol && ~s_v
            if abs(pi-(bandphi(j+n)-bandphi(j)))<tol
                % there is a point opposite the target
                c = cos(2*bandphi(j)); s = sin(2*bandphi(j));
                s_v = symcheck(C,[c,s,0; s,-c,0; 0,0,1]);
                phi = bandphi(j);
            end
            j = j+1;
        end
    else
        phi = atan2(C(2),C(1));
        s_v=symcheck(C,[cos(2*phi),sin(2*phi),0; sin(2*phi),-cos(2*phi),0; 0,0,1]);
    end
    syms(5) = phi;
    if s_v
        if k==1
            syms(1)=2;
        else
            syms(1) = 5;
        end
    else
        syms(1) = 1;
    end
    if Opts.dispSym
        nam = symmetryName(syms);
        fprintf(['Symmetry Group: ' nam '\n'])
    end
    self.symmetry = syms;
    return
end

    function copl = iscoplanar(C,tol)
        if size(C,2) <= 3
            copl = 1;
	        return
        end
        copl = rank(C(:,2:end)-C(:,1),tol) <= size(C,1)-1;
        % Based on code by Brett Shoelson, Ph.D.
        % brett.shoelson@mathworks.com
    end

%% Identify all Potential Axes
T=convhull(C');
% T gives the indices of each triangular face of C in a (:,3) array
% find circumcenters of triangular faces:
D = cross(C(:,T(:,2))-C(:,T(:,1)), C(:,T(:,3))-C(:,T(:,1)));
D = D./vecnorm(D);
% Remove Duplicates - coplanar faces (on a sphere) have same circumcenter
[~,i,grps] = unique(tol*round(D'/tol), 'rows');
D = D(:,i); % this way the values of D aren't the rounded values from the previous line
m = size(D,2);
self.dual = D;

fo = cell(m,1); % Face Order: essentially the same as T, but incorporates
                % non-triangular faces. The ith cell gives a list of vertex
                % indices that make up the ith face of the solid
E = zeros(3,m+n-2); % edge midpoint coordinates
ei = 1;
for i = 1:m
    IforT = grps==i;
    if sum(IforT)>1 % if more than one triangle from T makes up this face:
        verti = unique(T(IforT, :))'; %indices of all vertices that share a dual point
        % order the indices such that the vertices cycle
        % counterclockwise around the dual point:
        [~,sphV] = coorConvert(D(:,i),'cart',C(:,verti));
        [~,ord] = sort(atan2(sphV(3,:),sphV(2,:)));
        fo{i} = verti(ord);
    else
        fo{i} = T(IforT,:);
    end

    %find edge midpoints while we're here:
    l = length(fo{i});
    ord = [fo{i}; fo{i}(mod(1:l,l)+1)];
    ord = ord(:,find(ord(1,:)>ord(2,:))); % to make sure each edge is only counted once
    l = size(ord,2);
    E(:,ei:ei+l-1) = 0.5*(C(:,ord(1,:)) + C(:,ord(2,:)));
    ei = ei+l;
end
self.faceOrder = fo;
E = E./vecnorm(E);
    
F = [C,D,E]; %full set of axes
% We need to remove points y for which -y is not in F, because a symmetric axis
% should pass through two features. Meanwhile, we also need to remove one of each 
% pair {x,-x} that appears in F, because they both correspond to the same axis.
negs = squeeze(vecnorm(F+permute(F,[1 3 2])) < tol); % (i,j) is true if F(:,i)=-F(:,j)
[~,i] = find(triu(negs,1));
i = [i; find(sum(negs,2)==0)];
F(:,i)=[];


%% Handle low-symmetry case
o = size(F,2);
if o==0
    %shape has no symmetry axes; but there could be a single mirror plane,
    %in which case it would pass through the center of mass
    com = mean(C,2);
    psi = atan2(com(2),com(1));
    the = acos(com(3)/vecnorm(com));
    C2 = spinmat([2,3],-[the,psi])*C;
    % if this shape has a mirror plane, every point will either lie on it
    % or have a buddy with equal z-coordinate in this rotated basis
    z = round(C2(3,:),3);
    i = 1;
    buddy = i;
    while isscalar(buddy) && i < n
        buddy = find(z(i)==z);
        i = i+1;
    end
    if i==n
        % this shape has no symmetries
        phi = 0;
        type = 1;
        nam = 'C_1';
    else
        phi = atan2(C2(2,buddy(1))+C2(2,buddy(2)),C2(1,buddy(1))+C2(1,buddy(2)));
        if symcheck(C2,[cos(2*phi),sin(2*phi),0; sin(2*phi),-cos(2*phi),0; 0,0,1])
            type = 2;
            nam = 'C_s';
        else
            type = 1;
            nam = 'C_1';
        end
    end
    self.symmetry = [type;1;psi;the;phi];
    if Opts.dispSym
        fprintf(['Symmetry Group:' nam '\n'])
    end
    return
end

%% Calculate symmetries of each axis
syms = zeros(5,o);
for i = 1:o
    syms([3 4],i) = [atan2(F(2,i),F(1,i)); acos(F(3,i))]; % [psi; theta]
    C2 = spinmat([2,3],-syms([4,3],i)')*C; %rotate C so F(:,i) is on north pole

    % find horizontal mirror symmetry
    s_h = symcheck(C2,[1 0 0; 0 1 0; 0 0 -1]);

    % find rotational symmetry:
    % k is the order of the smallest set of points occupying the same
    % z-value, (ignoring poles)
    C3 = C2(:,find(vecnorm(C2([1 2],:))>tol)); %poles removed
    [~,~,grp] = unique(round(C3(3,:),2));
    [k,mini] = min(arrayfun(@(a) sum(grp==a),1:max(grp)));

    found = false;
    while k~=1 && ~found
        if symcheck(C2,spinmat(3,2*pi/k))
            found = true;
        else
            % if the band with the fewest points is not a regular
            % polygon, it might have a rotational symmetry with a lower
            % value of k:
            fact=factor(k);
            k=k/fact(1);
        end
    end
    syms(2,i) = k;

    % find vertical mirror symmetry:
    band = C3([1 2], find(grp==mini));
    b = size(band,2);
    if b>k % in the vast majority of cases this will be false
        % sort members by phi:
        bandphi(1:2:2*b) = sort(atan2(band(2,:),band(1,:)));
        % add midpoints:
        bandphi(2:2:2*b) = 0.5*(bandphi(1:2:2*b)+[bandphi(3:2:2*b),bandphi(1)+2*pi]);
        % iterate through a pi/k wide slice looking for a mirror plane
        j = 1;
        s_v=false;
        while bandphi(j) < bandphi(1)+pi/k+tol && ~s_v && j<=b
            if abs(pi-(bandphi(j+b)-bandphi(j)))<tol % there is a point opposite the target
                c = cos(2*bandphi(j)); s = sin(2*bandphi(j));
                s_v = symcheck(C2,[c,s,0; s,-c,0; 0,0,1]);
                phi = bandphi(j);
            end
            j = j+1;
        end
    else
        % if there the rotational symmetry matches the number of points
        % in the band (i.e. the band is a regular polygon), there will be
        % a mirror plane that passes through one of the vertices.
        phi = atan2(band(2),band(1));
        s_v=symcheck(C2,[cos(2*phi),sin(2*phi),0; sin(2*phi),-cos(2*phi),0; 0,0,1]);
    end
    syms(5,i) = phi;

    %perform any additional checks required to determine this axis'
    %symmetry:
    switch bin2dec(num2str([s_h s_v]))
        case 0
            % axis has neither vertical nor horizonal mirrors
            % options: C_k, R_k, D_k
            p_2k = pi/syms(2,i);
            if symcheck(C2,[cos(p_2k),-sin(p_2k),0; sin(p_2k),cos(p_2k),0; 0,0,-1])
                % axis has rotoreflection and so is R_k
                syms(1,i) = 4;
            else
                % Distinguish between C_k and D_k by searching for
                % perpendicular axes:
                perpi = find(abs(F(:,i).'*F)<tol);
                syms(1,i)=1;
                if length(perpi)>=k && k>1
                    thetahat = [cos(syms(3,i))*cos(syms(4,i)),sin(syms(3,i))*cos(syms(4,i)),-sin(syms(4,i))];
                    psihat = [-sin(syms(3,i)),cos(syms(3,i)),0];
                    phi = unique(sort(atan2(psihat*F(:,perpi),thetahat*F(:,perpi))));
                    phii=1;
                    found = false;
                    while phii<length(phi) && phi(phii)<2*pi/k+phi(1)+tol && ~found
                        found = symcheck(spinmat(3,-phi(phii))*C2,diag([1,-1,-1]));
                        phii=phii+1;
                    end
                    if found
                        syms(1,i)=5;
                        syms(5,i)=phi(phii-1);
                    end
                end
            end
        case 1
            % axis has vertical mirror but no horizontal mirror
            % options: C_kv, D_kd
            p_2k = pi/syms(2,i);
            if symcheck(C2,[cos(p_2k),-sin(p_2k),0; sin(p_2k),cos(p_2k),0; 0,0,-1])
                % axis has rotoreflection and so is D_kd
                syms(1,i) = 6;
            else
                % axis is C_kv
                syms(1,i) = 2;
            end
        case 2
            % axis has horizontal mirror but no vertical mirror
            % option: C_kh
            syms(1,i) = 3;
        case 3
            % axis has both horizontal and vertical mirrors
            % option: D_kh
            syms(1,i) = 7;
    end

end

%% Identify Symmetry Group from Axes
c1i = find(syms(2,:)==1); %indices of axes with no rotational symmetry
if length(c1i)==o
    if syms(1)==1
        scom = coorConvert(mean(C,2)); %Since it is technically possible for a symmetry
        % -less shape to have avoided the earlier check. Ask me how I know.
        self.symmetry = [1;1;scom(3);scom(2);0];
    else
        self.symmetry = syms(:,1);
    end
    if Opts.dispSym
        if syms(1) == 4
            fprintf('Symmetry Group: C_i\n')
        elseif syms(1)==2
            fprintf('Symmetry Group: C_s\n')
        else
            fprintf('Symmetry Group: C_1\n')
        end
    end
    return
end

syms(:,c1i) = [];
syms = sortrows(syms.',[2 1],'descend').';

reference = {[1 1 1 1; 3 3 3 3]; %T
             [4 4 4 4; 3 3 3 3]; %T_h
             [2 2 2 2; 3 3 3 3]; %T_d
             [5 5 5 5 5 5 5; 4 4 4 3 3 3 3]; %O
             [7 7 7 6 6 6 6; 4 4 4 3 3 3 3]; %O_h
             [5*ones(1,16); 5*ones(1,6) 3*ones(1,10)]; %I
             [6*ones(1,16); 5*ones(1,6) 3*ones(1,10)]};%I_h

poly = find(cellfun(@(c) isequal(c,syms([1,2],find(syms(2,:)>2))),reference),1);
if ~isempty(poly)
    % adjust phi value:
    % I don't know how to explain how, but this makes sure the visuals
    % are lined up correctly, so just trust it
    if dot([cos(syms(3))*sin(syms(4));sin(syms(3))*sin(syms(4));cos(syms(4))],...
        [cos(syms(8))*sin(syms(9));sin(syms(8))*sin(syms(9));cos(syms(9))])>0
        syms(8) = syms(8)+pi;
        syms(9) = pi-syms(9);
    end
    syms(5) = atan2(sin(syms(8)-syms(3)), cos(syms(4))*cos(syms(8)-syms(3))-sin(syms(4))*cot(syms(9)) );
    syms = [poly+7;0;syms(3:5,1)];
else
    syms = syms(:,1);
end

if Opts.dispSym
    nam = symmetryName(syms([1,2]));
    fprintf(['Symmetry Group: ' nam '\n'])
end
self.symmetry = syms;

end %getSymmetry

function rate(self,Opts)
% Evaluates the rating of each shape on a variety of different scales:
% Force Score - average of the magnitudes of the radial components of the
%               forces - unintuitively, this is the same as the total potential
%               energy divided by the number of points; smaller=better
% Distance Score - smallest distance between two vertices; larger=better
arguments
    self shape
    Opts.disp = true
    % When true, displays ratings in command window
end
N = length(self);
if N>1
    fprintf('Rating Distributions:')
    for a=1:N
        self(a).rate('disp',false)
        if round( floor(10*a/N)-floor(10*((a-1)/N)) )
            fprintf(' %%')
        end
    end
    fprintf(' Done!\n')
    return
end

C = self.cartCoor;
%calculate Force Score:
[~,~,~,R]=CS(C);
FScore = mean(vecnorm(R));

%calculate distance score
d = permute( vecnorm(C-permute(C,[1 3 2])), [3 2 1]);
d = d + 2*diag(ones(1,size(C,2))); %so it ignores that the distance to itself is 0
DScore = min(d,[],'all');

self.ratings = [FScore, DScore];
if Opts.disp
    fprintf('Force Score:%f\nDistance Score:%f\n',FScore,DScore)
end
end% rate

function tf = isCongruentTo(self,othr,tol)
% Returns true if the coordinates of each shape can be isometrically
% transformed into one another
arguments
    self (1,1) shape
    othr(1,1) shape
    tol {mustBeScalarOrEmpty} = 0.01
end
if isempty(self.dual)
    self.getSymmetry('dispSym',false)
end
if isempty(othr.dual)
    othr.getSymmetry('dispSym',false)
end
sym1 = self.symmetry;
sym2 = othr.symmetry;
tf = size(self.cartCoor,2)==size(othr.cartCoor,2) && isequal(sym1([1 2]),sym2([1 2]));
if tf
    hori = diag([1,1,-1]);
    vert = diag([1,-1,1]);
    I = eye(3);
    C1 = spinmat([3,2,3],-sym1([5,4,3])')*self.cartCoor;
    C2 = spinmat([3,2,3],-sym2([5,4,3])')*othr.cartCoor;
    symchk = @(M) all(ismember(tol*round(C1.'/tol), tol*round((M*C2).'/tol),'rows'));
    switch sym1(1)
        case 1 % C_k
            if sym1(2)==1
                % No symmetry
                [~,ord1] = sort(C1(3,:),'descend');
                [~,ord2] = sort(C2(3,:),'descend');
                if C1(3,ord1(1))>0.999
                    phi1=atan2(C1(2,ord1(2)),C1(1,ord1(2)));
                    phi2=atan2(C2(2,ord2(2)),C2(1,ord2(2)));
                else
                    phi1=atan2(C1(2,ord1(1)),C1(1,ord1(1)));
                    phi2=atan2(C2(2,ord2(1)),C2(1,ord2(1)));
                end
                M = [cos(2*phi1),sin(2*phi1),0; sin(2*phi1),-cos(2*phi1),0;0,0,1];
                tf = symchk(spinmat(3,phi1-phi2)) || symchk(M*spinmat(3,phi1-phi2));
            else
                tf = any([symchk(I),symchk(vert),symchk(hori),symchk(hori*vert)]);
            end
        case 2 %C_kv
            if sym1(2)==1
                %single mirror
                tf = symchk(I) || symchk(diag([-1,1,1]));
                %its a miracle to me that this check could be so simple
            else
                M = spinmat(3,pi/sym1(2));
                tf = any([symchk(I),symchk(hori),symchk(M),symchk(hori*M)]);
            end
        case 6 %D_kd
            tf = symchk(I) || symchk(hori);
        case 7 %D_kh
            if sym1(2)==2
                % three perpendicular mirrors, and any of their three lines
                % of intersection could be chosen as the principle axis
                M = spinmat(3,pi/2);
                tf = any([symchk(I),symchk(spinmat(2,pi/2)),symchk(spinmat(1,pi/2)), ...
                    symchk(M),symchk(M*spinmat(2,pi/2)),symchk(M*spinmat(1,pi/2))]);
            else
                tf = symchk(I) || symchk(spinmat(3,pi/sym1(2)));
            end
        case 8 %T
                M = hori*spinmat(3,pi/sym1(2));
            tf = any([symchk(I), symchk(M), symchk(vert), symchk(vert*M)]);
        case 10 %T_d
            tf = symchk(I) || symchk(hori*spinmat(3,pi/3));
        case {3,4,5,9,11,13} %C_kh, R_k, D_k, T_h, O or I
            tf = symchk(I) || symchk(vert);
            % it is plausible that for C_kh this is an insufficient check
        case {12,14} %O_h or I_h
            tf = symchk(I);
    end
end
end

function simulate(self, Opts)
% Simulates the motion of the points by treating each as an electron
arguments
    self (1,1) shape
    Opts.vel logical = false % When true, force changes points' velocities,
                             % otherwise, force changes poitns' positions.
    Opts.dt = 0.01           % Step size
    Opts.applyChanges logical = false % When true, the final coordinates of the
                             % simulation will be set to the shape's cartCoor
    Opts.tol = 1e-3          % value of Cross Sum at which the simulation halts
end

%set up axes:
if strcmp(gcf().UserData,'simulation')
    clf
    fig = gcf;
else
    fig = figure('Position',[80,430,850,420],'Color','w','UserData','simulation');
end
ax = axes(fig,'units', 'pixels', 'Position', [40,60,300,300]);
axis(ax,'equal','off')
hold(ax,"on")
ttl = uicontrol(fig,'style', 'text', 'Units', 'pixels', ... 
    'position', [540,230,110,20], 'String', 'Perturbation','BackgroundColor','w');
perSld = uicontrol(fig,'style','slider', 'Units', 'pixels','BackgroundColor','w', ...
    'position',[385,200,425,20],'min', 0, 'max', 0.2, 'Value', 0);
addlistener(perSld, 'Value', 'PostSet', @perSldcallback);
pb = uicontrol(fig,'Style','pushbutton', 'Units', 'pixels','BackgroundColor','w', ...
    'Position', [15,20,88,30], 'String', 'Simulate', 'Callback', @pbcallback);
% simSld = uicontrol(fig,'Style','slider','Units','pixels','BackgroundColor','w', ...
%     'Position',[110,20,360,30],'Enable','off');
% Frams = getframe(ax);

% draw initial configuration:
C = self.cartCoor;
S = coorConvert(C,'cart');
n = size(C,2);
points = scatter3(ax,C(1,:),C(2,:),C(3,:), 50, 'k', 'filled');
[x,y,z] = sphere;
sph = surf(ax,x,y,z,'EdgeColor',[0 0.9 0.9], 'FaceAlpha',0.1, 'FaceColor','c');
view(ax,[0.05 -1 .2])
gamma = rand(1,n)*2*pi;
delta = 0;


    function pbcallback(~, ~)
    S = S + [zeros(1,n); delta*[-cos(gamma); sin(gamma)]];
    C = coorConvert(S,'sph');

    % get rid of old structures and insert graph:
    delete([perSld, ttl])
    set(pb,'Style','togglebutton','String','stop','Value',0,'Callback',@(~,~)[])

    if Opts.vel
        program = {@updateWithVel,@drawPos,@plotNRG};
        vel = zeros(3,n);
        lname = 'Kinetic Energy';
    else
        program = {@updatePos,@drawPos,@plotCS};
        lname = 'log(CS)';
    end
    totalSteps = 1e4;
    div = 1;
    step = 0;
    t = 1:totalSteps/div;
    f = nan(size(t));
    g = axes(fig,'Units','pixels','Position',[390,120,420,200],'YGrid','on');
    hold(g,'on')
    l=plot(g,t,f,'b','DisplayName',lname);
    legend(g)
    [cs,T] = CS(C);
    
    % actually simulate:
    quitthis = false;
    while all([cs > Opts.tol, ~quitthis, step<=totalSteps])
        for i=1:3
            feval(program{i})
        end
        Frams(step) = getframe(ax);
        quitthis = get(pb, 'Value');
    end

    % post-simulation stuff:
    Tri = convhull(C');
    trisurf(Tri, C(1,:), C(2,:), C(3,:),ones(1,size(Tri,1)),'Parent',ax)
    delete(sph)
    set(pb,'Enable','off')
    % Frams(step+1) = getframe(ax);
    % set(simSld,'Enable','on','min',1,'max',step+1,'Value',step+1)
    % addlistener(simSld, 'Value', 'PostSet', @simSldcallback);
    if Opts.applyChanges
        self.cartCoor = C;
        self.minOrd   = -log10(cs);
    end

        function updateWithVel
            step = step+1;
            vel = vel + T*Opts.dt;
            C = C + vel*Opts.dt;
            C = C./vecnorm(C);
            [cs,T] = CS(C);
        end
        %%%%%%%%%%%%
        function updatePos
            step = step+1;
            C = C + T*Opts.dt;
            C = C./vecnorm(C);
            [cs,T] = CS(C);
        end
        %%%%%%%%%%
        function drawPos
            if mod(step,div)==0
                set(points,'XData',C(1,:),'YData',C(2,:),'ZData',C(3,:))
                drawnow
            end
        end
        %%%%%%%%%%%
        function plotCS
            if mod(step,div)==0
                l.YData(step/div)= log10(cs);
            end
        end
        %%%%%%%%%%%
        function plotNRG
            if mod(step,div)==0
                K = 0.5*sum(vel.^2,"all");
                l.YData(step/div)= K;
            end
        end
    end

    function perSldcallback(~,eventdata)
        delta = get(eventdata.AffectedObject, 'Value');
        C = coorConvert(S + [zeros(1,n); delta*[-cos(gamma); sin(gamma)]],'sph');
        points.XData = C(1,:); points.YData = C(2,:); points.ZData = C(3,:);
    end

    % function simSldcallback(~,eventdata)
    %     i = round(get(eventdata.AffectedObject, 'Value'));
    %     cla
    %     imshow(frame2im(Frams(i)),'Parent',ax)
    % end
end

function [self,counts]= removeDupes(self)
% Removes instances of congruent shapes from the shape vector self, leaving only the
% one instance with the lowest Cross Sum.
% counts(i) gives how many instances of self(i) (after duplicate removal) were
% present in self originally.
refi = 1;
while refi < length(self)
    counts(refi) = 1;
    chki = refi+1;
    while chki <= length(self)
        if self(refi).isCongruentTo(self(chki))
            if self(refi).minOrd < self(chki).minOrd
                self(refi) = self(chki).copySelf;
            end
            self(chki)=[];
            counts(refi)=counts(refi)+1;
        else
            chki = chki+1;
        end
    end
    refi = refi+1;
end
end

function other = copySelf(self)
% Sets all the properties of other equal to that of self.
% Avoids the strange reference issues that come with setting
% two objects equal to one anoter directly.
if length(self)>1
    other = arrayfun(@(a) self(a).copySelf, 1:length(self));
    return
end
other = shape();

other.cartCoor  = self.cartCoor;

other.minOrd    = self.minOrd;
other.dual      = self.dual;
other.faceOrder = self.faceOrder;
other.symmetry  = self.symmetry;
other.ratings   = self.ratings;
end

%%%%%%%%%%%%%%
end%methods
end%classdef