classdef shape < handle
properties
    name     = {}
    numVerts
    cartCoor
    minVal

    structure
    thresholds = [0.003, 0.01, 0.1]
    fingerprint = cell(3,1)
    ratings
end %properties

methods
function [s,coorf,varvec,nums] = shape(name,minimize,Opts)
arguments 
    name char = ''
    minimize logical = true
    Opts.random logical = false
end
if ~isempty(name)
    %get numbers from name
    s.name= {name, [0 0 1]};
    if Opts.random
        if length(name)>1
            error('randomized distributions take only a single number name')
        end
        n=str2double(name);
        s.numVerts=n;
        s.cartCoor=sphere2cart([ones(n,1) rand(n,1)*pi rand(n,1)*2*pi]);
        return
    end
    [nums,bands,primes,NPole,SPole]=nameInterpreter(name);
    B = length(bands);
    n = sum(bands)+NPole+SPole;
    s.numVerts=n;
    
    % assemble the shape
    phi = mod([zeros(NPole,1); cell2mat(arrayfun(@(b)(0:bands(b)-1)'*2*pi/bands(b)+ ...
        primes(b)*pi/bands(b),1:B,'un',0)'); zeros(SPole,1)],2*pi);
    vari= cell2mat(arrayfun(@(b)b*ones(bands(b),1), 1:B, 'un',0)');

    %minimize the setup:
    caps={[],0,[],pi};
    coorf=@(vars) sphere2cart([ones(n,1), [caps{NPole+1};vars(vari);caps{SPole+3}] phi]);
    varvec = (1:B)'*pi/(B+1);
    if minimize
        varvec = gradMin(@(v)CS(coorf(v)),varvec);
    end
    s.cartCoor = coorf(varvec);
    s.minVal = CS(s.cartCoor);
end
    
end
%%%%%%%%%%%%%

function [tdata,rdata] = forceMin(self, e, Opts)
arguments
    self shape
    e double = 1e-4
    Opts.panicMode = 0
    Opts.dispVal logical = true
end
N = length(self);
if N>1
    fprintf('Minimizing Distributions:')
    for a = 1:N
        self(a).forceMin(e,'panic',Opts.panicMode,'disp',Opts.dispVal)
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
    rdata = mag(R);
end
pm = ones(self.numVerts,1);
nxt = @(c,t,f,m) cos(t*m).*c + sin(t*m).*f./m; 
ohno = false;
while cs > e && ~ohno
    p = min([0.1, cs*1e-3]); %precision increases as cs decreases
    cs = f(C,p);
    dtCopy = dt;
    %shrink step-size until it's clear whether to go forwards or backwards
    dir = [f(nxt(C,-dt,T,M),p), cs, f(nxt(C,dt,T,M),p)];
    [~,fb] = min(dir); fb = fb-2;
    while ( fb == 0 || any(diff(dir)==0) ) && ~ohno
        dt = dt*0.5;
        ohno = dt<e;
        dir = [f(nxt(C,-dt,T,M),p), cs, f(nxt(C,dt,T,M),p)];
        [~,fb] = min(dir); fb = fb-2;
    end
    if ohno
        switch Opts.panicMode
            case 1
                dt = dtCopy;
                T = pm.*T; %undoes the redirection
                [pm, ohno] = adjustPM(C);
                fb = 1;
                T = pm.*T;
                % disp(pm')
            case 2
                f = @(c,p) p*round(coulombPotential(c)/p);
                dt = dtCopy;
                ohno = false;
            case 3
                f = @(c,p) p*round(mag(mean(c))/p);
                dt = dtCopy;
                ohno = false;
            otherwise
                break
        end
    end
    t = 0;
    % take as many steps as keeps the score decreasing
    steps = 0;
    while f(nxt(C,t+fb*dt,T,M),p) < f(nxt(C,t,T,M),p)
        t=t+fb*dt;
        steps = steps+1;
    end
    if steps>10; dt = dt*3; end
    C=nxt(C,t,T,M);
    [cs,T,M,R] = CS(C);
    if nargout
        tdata = [tdata M];
        rdata = [rdata mag(R)];
    end
    T = pm.*T;
end
self.cartCoor = C./mag(C);
self.minVal = cs;
if Opts.dispVal; disp(cs);end


    function [pm, ohno] = adjustPM(C)
        [~,ord] = sort(M);
        magCopy = M/M(ord(end));
        grps=diff([0 find(diff(magCopy(ord)')>0.1) size(C,1)]);
        G=length(grps);

        pm = double(dec2bin(0:2^G-1)=='1')'*2-1;
        pm = cell2mat(arrayfun(@(a) ones(grps(a),2^G).*pm(a,:), 1:G, 'un',0)');
        pm = pm(arrayfun(@(a) find(ord==a), 1:size(C,1)),:);
        [~,mini]=min(arrayfun(@(a) f(nxt(C,dt,pm(:,a).*T,M),p), 1:2^G));
        pm = pm(:,mini);
        ohno = f(C,p)<=f(nxt(C,dt,pm.*T,M),p);
    end

end%forceMin

function see(self,style, Opts)
% Input names: axes, style
arguments
    self (1,1) shape
    style char = 'normal'
    Opts.axes matlab.graphics.axis.Axes = gca
    Opts.coor = self.cartCoor
    Opts.clear logical = true
    Opts.mainColor = 'k'
    Opts.dualColor = 'r'
end
validNames={'normal','solid','area','highlight','bands', ...
    'force','potential','flat','contourPlot'};
style = validatestring(style, validNames);

if length(self) > 1
    viewManyShapes(self, "visType", style)
end
trgAx = Opts.axes;
view(trgAx, [0.05 -1 .2])
axis(trgAx, 'equal','off');
if Opts.clear; cla(trgAx); end
hold(trgAx, 'on')
feval(eval(['@' style]))
hold(trgAx, 'off')

function normal
    % normal Points-on-a-Sphere visual
    [x,y,z] = sphere;
    surf(trgAx, x,y,z, 'HitTest','off')
    C = Opts.coor;
    scatter3(trgAx, C(:,1), C(:,2), C(:,3), 100, Opts.mainColor, ...
        'filled','HitTest','off')
end
function solid
    if isempty(self.structure) || ~isequal(self.cartCoor,Opts.coor)
        self.getStructure('coor',Opts.coor);
    end
    if Opts.clear
        [x,y,z]=sphere(15);
        surf(trgAx,x,y,z, 'EdgeColor',[0 0.9 0.9], 'FaceAlpha',0.1, 'FaceColor','c')
    end
    D = self.structure{1};
    fo = self.structure{4};
    f = @(a,i) Opts.coor(fo{a},i);
    arrayfun(@(a) patch(trgAx,f(a,1),f(a,2),f(a,3),'w', ...
        'EdgeColor',Opts.mainColor), 1:size(D,1));
end
function area
    %circles bounded by the distributions' voronoi regions
    [x,y,z]=sphere(20);
    surf(x,y,z, 'FaceColor',[.7,.7 .7],'EdgeColor','none');
    [circles, regions] = self.rate;
    C = Opts.coor;
    if isempty(self.structure)
        self.getStructure('coor',Opts.coor)
    end
    D = self.structure{1};
    scatter3(C(:,1), C(:,2), C(:,3), 100, Opts.mainColor, 'filled')
    scatter3(D(:,1), D(:,2), D(:,3), 100, Opts.dualColor, 'filled')
    for i = 1:self.numVerts
        c=circles{i}; r = regions{i};
        plot3(c(:,1), c(:,2), c(:,3), 'b', 'LineWidth', 2)
        plot3(r(:,1), r(:,2), r(:,3), Opts.dualColor, 'LineWidth', 2)
    end
end
function highlight
    %highlight distinct faces and edges
    if isempty(self.structure)
        self.getStructure('coor',Opts.coor);
    end
    fo = self.structure{4};
    allUni = self.structure{5};
    uf = allUni{2};
    ue = allUni{3};
    C = Opts.coor;
    %Draw each set of similar faces in different colors:
    nU=length(uf);
    colr = (0:nU-1)/nU;
    for fi = 1:nU
        f = @(a,i) C(fo{uf{fi}(a)}, i);
        arrayfun(@(a) patch(trgAx,f(a,1),f(a,2),f(a,3),colr(fi), ...
            'EdgeColor',"none"), 1:length(uf{fi}));
    end
    
    %Draw each set of equal edges in different colors
    nU=length(ue);
    colr = 1000*(1-(1:nU)/nU)'; %makes sure the colors are spread out
    colr = round([colr/100 mod(colr,100)/10 mod(colr,10)])/10;
    edges = self.structure{2};
    uEdge = zeros(1,nU); d=cell(1,nU);
    for ei = 1:nU
        f = @(a,i) C(edges([1,2],ue{ei}(a)),i);
        pl = arrayfun(@(a) plot3(trgAx, f(a,1),f(a,2),f(a,3), 'Color',colr(ei,:), ...
            'LineWidth',4), 1:length(ue{ei}));
        uEdge(ei) = pl(1); 
        v = [f(1,1) f(1,2) f(1,3)];
        d{ei} = num2str( mag(v(1,:)-v(2,:)) );
    end
    S = cart2sphere(C);
    C=sphere2cart([1.1*ones(self.numVerts,1) S(:,[2 3])]);
    % text(C(:,1),C(:,2),C(:,3),cellstr(num2str((1:self.numVerts)')),"FontSize",15)
end
function bands
    %shows the parallel bands which make up a solution
    [x,y,z]=sphere(15);
    surf(trgAx,x,y,z, 'EdgeColor',[0 0.9 0.9], 'FaceAlpha',0.1, 'FaceColor','c')
    colr='rgby';
    txtbx=uicontrol('Style','popupmenu','String',self.name(:,1),'Units','normalized', ...
        'Position',[0.05 0.75 0.1 0.1], 'Value',1,'Callback',@dropdowncb);
    C = Opts.coor;
    h=scatter3(trgAx, C(:,1),C(:,2),C(:,3),100,Opts.mainColor,'filled');
    set(h, 'DeleteFcn', @(~,~) delete(txtbx))
    plots = drawplanes(1);
    
    function planes=drawplanes(i)
        [~,nums,~,~,bms,order]=bandReorder(C,self.name{i,2}, self.thresholds(3));
        p = @(b,k) C(order(bms(b)),k);
        planes = arrayfun(@(b) patch(trgAx,p(b,1),p(b,2),p(b,3),colr(i),'FaceAlpha',0.5, ...
                'EdgeColor',colr(i)), 1:length(nums) );
    end

    function dropdowncb(src,~)
        delete(plots)
        hold(trgAx,"on")
        plots=drawplanes(src.Value);
        hold(trgAx,"off")
    end
end
function force
    C = Opts.coor;
    [~,T] = CS(C);
    if isempty(self.structure)
        self.getStructure('coor',C)
    end
    D = self.structure{1};
    fo=self.structure{4};
    f = @(a,i) C(fo{a},i);
    arrayfun(@(a) patch(trgAx,f(a,1),f(a,2),f(a,3),'w','EdgeColor',Opts.mainColor, ...
        'FaceAlpha',0.6, 'LineWidth',2), 1:size(D,1));
    f = @(a,i) [C(a,i), C(a,i)+T(a,i)];
    arrayfun(@(a) plot3(trgAx, f(a,1),f(a,2),f(a,3), 'r','LineWidth',2), 1:self.numVerts)
end
function flat
    S=cart2sphere(Opts.coor);
    th = find(pi-S(:,2)<1e-2,1);
    SPole=~isempty(th);
    rho=2*tan(S(:,2)/2);
    if SPole
        rho(th)=NaN;
    end
    x=rho.*cos(S(:,3));
    y=rho.*sin(S(:,3));
    scatter(trgAx,x,y,50,Opts.mainColor,'filled')
    if isempty(self.structure)
        self.getStructure;
    end
    fo = self.structure{4};
    arrayfun(@(a) patch(trgAx,x(fo{a}),y(fo{a}),'w','EdgeColor',Opts.mainColor, ...
        'FaceAlpha',0), 1:size(self.structure{1},1));
    view(trgAx,[0 0 1])
    if SPole
        t=linspace(0,2*pi,40);
        r=1.15*max(rho);
        plot(trgAx,r*cos(t),r*sin(t),'k')
        edge = self.structure{2}';
        nb= [edge(find(edge(:,1)==th),2); edge(find(edge(:,2)==th),1)];
        arrayfun(@(a) plot(trgAx,x(nb(a))*[rho(nb(a)) r]/rho(nb(a)), ...
            y(nb(a))*[rho(nb(a)) r]/rho(nb(a)),'k'), 1:length(nb))
    end
end
function contourPlot
    res=100;
    th=linspace(0,pi,res+1);
    ph=linspace(0,2*pi,res+1);
    [theta,phi] = meshgrid(th,ph);
    U = arrayfun(@(t,p) coulombPotential(Opts.coor,sphere2cart([1,t,p])),theta,phi);
    [~,edges]=histcounts(U(U<1.5*self.numVerts));
    M = contourc(th,ph,U,edges);
    i=1;
    while i<length(M)
        Ni=M(2,i);
        C=sphere2cart([ones(Ni,1) M(:,i+1:i+Ni)']);
        plot3(trgAx,C(:,1),C(:,2),C(:,3),'k')
        i = i+Ni+1;
    end

    [x,y,z]=sphere(res);
    surf(trgAx,x,y,z,'FaceColor','w','EdgeColor','none')
end

end%see

function getStructure(self,dualThresh, uniqThresh,Opts)
arguments
    self
    dualThresh = 0.003
    uniqThresh = 0.01
    Opts.coor = self(1).cartCoor
end
% Structure is organized as:
% 1) Dual vertices
% 2) Edge indices
% 3) Voronoi Order
% 4) Face Orders
% 5) Unique vertices, faces, and edges
if length(self)>1
    for s=1:length(self)
        self(s).getStructure(dualThresh,uniqThresh)
    end
    return
end

self.thresholds(1) = dualThresh;
self.thresholds(2) = uniqThresh;
C = Opts.coor;
n = self.numVerts;
if iscoplanar(C, 0.01)
    C = C(ccReorder(C),:); self.cartCoor=C;
    T=[1:n; n:-1:1]; flat=true;
else
    T=convhull(C); flat=false;
end
% find circumcenters of of triangular faces of convhull
D1 = cross(C(T(:,2),:)-C(T(:,1),:), C(T(:,3),:)-C(T(:,1),:));
D1 = D1./mag(D1);
%Remove Duplicates
[~,~,grps] = unique(dualThresh*round(D1/dualThresh), 'rows');
m = max(grps);
D = cell2mat(arrayfun(@(g) D1(find(g==grps,1),:), 1:m, 'un',0)');

assocmat = false(n, m); %(i,j) will be true if the ith vert is part of the jth face
fo = cell(m, 1);% face order: lists the indices of the vertices which make up
                % each face, in an order that makes a counter-clockwise circle
for facei = 1:m
    IforT = grps==facei;
    if sum(IforT)>1 % if more than one triangle from T makes up this face:
        IofPinQ = unique(T(IforT, :),'stable'); %indices of points in question
        order=ccReorder(C(IofPinQ,:));
        fo{facei} = IofPinQ(order)';
    else
        fo{facei} = T(IforT,:);
    end
    assocmat(fo{facei},facei) = true;
end
%Now do the same thing, but backwards, finding the 'voronoi order' of the
%indices of the points of the dual that will draw a cycle around each vertex.
%Also take this time to identify which pairs of vertices make up the edges
vo = cell(n, 1);
if flat
    vo = mat2cell(repmat([1 2],[n,1]), ones(n,1),2);
    edge=[fo{1}; fo{1}(rem(1:length(fo{1}),length(fo{1}))+1); ones(1,m+n-2); 2*ones(1,m+n-2)];
else
    edge=zeros(4,m+n-2); e=1;
    for verti = 1:n
        IofPinQ = find(assocmat(verti,:)); %indices of points in question
        order = ccReorder(D(IofPinQ,:));
        vo{verti} = IofPinQ(order);
    
        edgei=find(sum(assocmat(:,IofPinQ),2)==2)'; %i's of verts connected to the target
        edgei(edgei<verti)=[]; l=length(edgei);
        edge([1,2],e:e+l-1)=[verti*ones(1,l); edgei];
        e=e+l;
    end
    % find which pair of faces are associated with each edge
    edge([3,4],:)=cell2mat(arrayfun(@(a) find( ...
        assocmat(edge(1,a),:)&assocmat(edge(2,a),:))',1:m+n-2,'un',0));
end
A = faceAreas(C,fo);
eLengths = arrayfun(@(a) mag(C(edge(1,a),:)-C(edge(2,a),:)), 1:m+n-2);
allProps = cell(3,1);

%edge properties: the length of each edge plus the sum of the areas of the
%faces it touches
allProps{3}=arrayfun(@(a) eLengths(a)+sum(A(edge([3,4],a))), 1:m+n-2);

%vert properties: the sum of the lengths of the edges which border it plus
%the areas of the faces it borders
allProps{1}=arrayfun(@(a) sum(eLengths(edge(1,:)==a | edge(2,:)==a)) ...
    + sum(A(assocmat(a,:))), 1:n);

%face properties: the sum of the lengths of the edges which border it plus
%its area plus the number of neighbors its vertices have
allProps{2}=arrayfun(@(a) sum(eLengths(edge(3,:)==a | edge(4,:)==a)) + A(a) ...
    +sum(assocmat(find(assocmat(:,a)),:),'all'), 1:m);


allUni  =cell(3,1);
for u =1:3
    [~,~,group]=unique(uniqThresh*round(allProps{u}/uniqThresh));
    i = 1:length(allProps{u});
    allUni{u}=arrayfun(@(a) i(group==a), 1:max(group),'un',0);
    self.fingerprint{u} = sort(cellfun(@length, allUni{u}));
end

self.structure = {D,edge, vo,fo,allUni};


    function copl = iscoplanar(C,tol)
    if size(C,1) <= 3
	    copl = 1;
	    return
    end
    rnk = rank(bsxfun(@minus,C(2:end,:),C(1,:)),tol);
    copl = rnk <= size(C,2) - 1;
    % Based on a function by Brett Shoelson, Ph.D.
    % brett.shoelson@mathworks.com
    end

    function a = faceAreas(C,fo)
        a=zeros(length(fo),1);
        for f = 1:length(fo)
            p=[C(fo{f},:); C(fo{f}(1),:)];
            i1 = 1:size(p,1)-1;
            %surprisingly neat formula for the area of a planar polygon
            a(f) = 0.5*sum(mag( cross(p(1,:)-p(i1,:),p(i1+1,:)-p(i1,:),2) ));
        end
    end
end

function [circles, regions] = rate(self)
if length(self) > 1
    for s = 1:length(self)
        self(s).rate;
    end
    return
end

C = self.cartCoor;
%calculate Force Score:
[~,~,~,R]=CS(C);
FScore = mean(mag(R));
    
%calculate Area Score
if isempty(self.structure)
    self.getStructure;
end
D = self.structure{1};
vo = self.structure{3};
n = self.numVerts;
circles = cell(1,n); regions = cell(1,n);
res = 20;
area = zeros(n,1);
for i = 1:n
    if size(vo{1},2) > 2
        numedges = length(vo{i});
        regions{i} = zeros(numedges*res, 3);
        for e = 1:numedges
            OP = D(vo{i}(e), :);
            OQ = D(vo{i}(rem(e,numedges)+1), :);
            aPOQ = acos(dot(OP, OQ));
            t = linspace(0, aPOQ, res)';
            regions{i}((e-1)*res+1:e*res, :) = (sin(aPOQ-t)/sin(aPOQ)).*OP + (sin(t)/sin(aPOQ)).*OQ;
        end
    else
        %for shapes on a plane
        t=linspace(0,pi,res)';
        mid = (C(i,:)+C(rem(i,n)+1,:))/2;
        mid = mid./mag(mid);
        regions{i} = D(1,:).*cos(t) + mid.*sin(t);
        numedges=1;
    end
    dist = acos(dot(repmat(C(i,:), [numedges*res,1]), regions{i}, 2));
    [~, npi] = min(dist); %index of the nearest point on the regions curve to vertex_i

    OP = C(i,:);
    OQ = regions{i}(npi,:);
    OA = dot(OP,OQ).*OP;
    AQ = OQ - OA;
    AB = cross(AQ, OP);
    t = linspace(0, 2*pi, 50)';
    circles{i} = OA + AQ.*cos(t)+AB.*sin(t);
    area(i) = 2*pi*(1 - dot(OP,OQ));
end
Ascore = sum(area)/(4*pi);

%calculate uniqueness score
allUni = self.structure{5};
UScore = length(allUni{1})+length(allUni{2})+length(allUni{3});
UScore = UScore/(2*n + 2*size(D,1)-2); %divide by total number of verts,faces,&edges

%calculate distance score
d = permute( mag( C-permute(C,[3 2 1]) ), [1 3 2]);
d = d + 2*diag(ones(1,n)); %so it ignores that the distance to itself is 0
DScore = min(d,[],'all');


self.ratings = [FScore, Ascore, UScore, DScore];
end

function tf = isSameAs(self, other)
%returns true if the two shapes have equivalent geometry
if isempty(self.structure)
    self.getStructure;
end
if isempty(other.structure)
    other.getStructure
end
tf = size(self.fingerprint,1)==size(other.fingerprint,1);
if tf
    tf = all(cellfun(@(c1,c2) isequal(c1,c2), self.fingerprint, other.fingerprint));
end
end

function simulate(self, Opts)
arguments
    self (1,1) shape
    Opts.mu = 0.5
    Opts.vel logical = false
    Opts.graph logical = true
    Opts.dt = 0.01
    Opts.track logical = false
    Opts.applyChanges logical = false
    Opts.potential logical = false
    Opts.lim = 1e-3
end

C = self.cartCoor;
S = cart2sphere(C);
figure('Position', [0 700 560 520])
ax = axes('units', 'pixels', 'Position', [130, 44, 300, 300]);
axis equal off
hold(ax,"on")
[x,y,z] = sphere;
sph = surf(ax,x,y,z);
view([0.05 -1 .2])
ttl = uicontrol('style', 'text', 'Units', 'normalized', ... 
    'position', [0.25, 0.85,  0.5, 0.05], 'String', 'Perturbation');
sld = uicontrol('style','slider', 'Units', 'normalized', ...
    'position',[0.25, 0.8, 0.5, 0.05],'min', 0, 'max', 0.2, 'Value', 0);
gamma = rand(self.numVerts,1)*2*pi;
delta = 0;
addlistener(sld, 'Value', 'PostSet', @sldcallback);

points = scatter3(C(:,1),C(:,2),C(:,3), 50, 'k', 'filled');

pb = uicontrol('Style','pushbutton', 'Units', 'normalized', ...
    'Position', [0.1, 0.8, 0.1, 0.05], 'String', 'simulate', 'Callback', @pbcallback);
n = self.numVerts;

    function pbcallback(~, ~)
    S = S + [zeros(n,1) delta*[-cos(gamma) sin(gamma)]];
    delete(pb); delete(sld)
    totalSteps = 1e4;
    div = 10;
    stopbutton = uicontrol('Style','togglebutton', 'Units', 'normalized', ...
         'Position', [0.1, 0.8, 0.1, 0.05], 'String', 'stop', 'Value', 0);
    vel = zeros(n, 3);
    step = 0;
    if Opts.vel
        program = {@updateWithVel,@drawPos};
    else
        program = {@updatePos,@drawPos};
    end
    [cs,T] = CS(C); T=cartvec2sphvec(S,T);
    if Opts.graph
        delete(ttl)
        t = 1:totalSteps/div;
        f = nan(size(t));
        g = axes('Units','pixels','Position',[150 375 260 125],'YGrid','on');
        hold(g,'on')
        l1=plot(g,t,f,'b');
        l2=plot(g,t,f,'r');
        l3=plot(g,t,f,'g');
        legend('cs','ts','com')
        program{3}=@plotCS;
    else
        ttl.String = num2str(cs);
        program{3} = @updateLabel;
    end
    if Opts.track
       history = scatter3(ax,C(:,1),C(:,2),C(:,3),20,'red');
       program{end+1}=@plotHistory;
    end

    quitthis = false;
    while all([cs > Opts.lim, ~quitthis, step<=totalSteps])
        for j=1:length(program)
            feval(program{j})
        end
        quitthis = get(stopbutton, 'Value');
    end
    axes(ax)
    trisurf(convhull(C), C(:,1), C(:,2), C(:,3))
    delete(sph); delete(stopbutton)
    hold off
    if Opts.applyChanges
        self.cartCoor = C;
        self.minVal   = csNew;
    end

        function updateWithVel
            step = step+1;
            vel = vel + Opts.mu*T*Opts.dt;
            S = S + vel*Opts.dt;
            C = sphere2cart(S);
        end
        %%%%%%%%%%%%
        function updatePos
            step = step+1;
            S = S+T*Opts.dt;
            C = sphere2cart(S);
        end
        %%%%%%%%%%
        function drawPos
            if mod(step,div)==0
                set(points,'XData',C(:,1),'YData',C(:,2),'ZData',C(:,3))
                [cs,T] = CS(C);
                T=cartvec2sphvec(S,T);
                drawnow limitrate
            end
        end
        %%%%%%%%%%%
        function plotCS
            if mod(step,div)==0
                ts = mag(sum(T));
                com = mag(mean(C,1));
                l1.YData(step/div)= log10(cs);
                l2.YData(step/div)= log10(ts);
                l3.YData(step/div)= log10(com);
            end
        end
        %%%%%%%%
        function updateLabel
            ttl.String = num2str(csNew);
            cs = csNew;
        end
        %%%%%%%
        function plotHistory
            if mod(step,100)==0
                history.XData = [history.XData C(:,1)'];
                history.YData = [history.YData C(:,2)'];
                history.ZData = [history.ZData C(:,3)'];
            end
        end
    end

    function sldcallback(~, eventdata)
    delta = get(eventdata.AffectedObject, 'Value');
    C = sphere2cart(S + [zeros(n,1) delta*[-cos(gamma) sin(gamma)]]);
    points.XData = C(:,1); points.YData = C(:,2); points.ZData = C(:,3);
    end
end

function clearHigherOrderAtts(self)
for s=1:length(self)
    self(s).structure=[];
    self(s).fingerprint=[];
    self(s).parents={};
    self(s).children={};
end
end

function rename(self, bandThresh)
% finds the names which best describe the shape's geometry
arguments
    self
    bandThresh = 0.1
end
if length(self)>1
    for s=1:length(self)
        self(s).rename(bandThresh)
    end
    return
end

if isempty(self.structure)
    self.getStructure;
end
D = self.structure{1};
edges = self.structure{2};
allUni = self.structure{5};
% The plan is to place put each vertex, face, and edge on the north pole
% and look for parallel bands of points. This section puts each of the
% reference points in a list.
uv = allUni{1}; uf = allUni{2}; ue = allUni{3};
C=self.cartCoor;
L=@(l) 1:length(l);
setupsToCheck= length([uv, uf, ue]);
names = cell(setupsToCheck,2);
names(:,2) = [arrayfun(@(v) C(uv{v}(1),:), L(uv), 'un',0)    ...
    arrayfun(@(f) D(uf{f}(1),:), L(uf), 'un',0)...
    arrayfun(@(e) mean(C(edges([1,2],ue{e}(1)),:),1), L(ue), 'un',0)];

numsSoFar=cell(setupsToCheck,1);
flags=false(setupsToCheck,1);
% namChk = @(nam,C) any(cellfun(@(c) isequal(nam,c), C));
for i = 1:setupsToCheck
    [C, nums, phi, bsi, bms,~,flpd] = bandReorder(self.cartCoor,names{i,2}, bandThresh);
    if flpd;names{i,2}=-names{i,2};end
    %get the midpoints of the edges of each band:
    mids=(phi+cell2mat(arrayfun(@(b) reCycle(phi(bms(b))), L(nums), 'un',0)'))/2;
    % Figure out which point within the reference band is the reference
    % point: do this by, for each point within the reference band, counting 
    % how many bands have a point or midpoint which lies near that point:
    fixi= banditofix(nums); % index of the band which is the reference
    f= @(k) phi(bsi(fixi)+k-1); %gives the value of the kth point in the ref band
    algn=arrayfun(@(a) sum(abs(fixPhi(phi,f(a))-f(a))<0.1 | ...
        abs(fixPhi(mids,f(a))-f(a))<0.1), 1:nums(fixi));
    [~,algn]=max(algn);
    phi=phi-phi(bsi(fixi)+algn-1); %rotates the shape so the ref point is at 0
    rotat=arrayfun(@(a) ~any(abs(fixPhi(phi(bms(a)),0))<0.1), L(nums));
    rotat(nums==1)=false;

    names{i,1}=nameConstructor(nums,rotat);

    % check whether this name should be thrown out:
    % flags the name if any of the bands have a center of mass far from 0, or
    % if any of the bands are rotated out of accordance with their primestatus, or 
    % if the shape name is symmetric yet the z-components aren't
    flags(i)= any([... 
        any(arrayfun(@(a) sum(mean(C(bms(a),1:2)).^2)>bandThresh^2 , L(nums))) ...
        any(arrayfun(@(a) rotatChk(phi(bms(a)),rotat(a),bandThresh), L(nums))) ...
        (isequal(nums(end:-1:1),nums) && abs(mean(C(:,3)))>bandThresh) ]);

    numsSoFar{i}=nums;
end
if all(flags)
    self.rename(bandThresh*1.5)
    return
end
names(flags,:)=[];
[~,i]=unique(names(:,1));
names=names(i,:);
%re-order the names in decreasing levels of rotational symmetry:
numsSoFar = numsSoFar(~flags);
[~,i]=sort(cellfun(@rotatSym, numsSoFar(i)),'descend');
self.name=names(i,:);
self.thresholds(3) = bandThresh;
%rotate everything for neatness
NPole = self.name{1,2};
self.cartCoor = north(self.cartCoor,NPole);
self.name(:,2) = arrayfun(@(a)north(self.name{a,2},NPole),(1:size(self.name,1))','un',0);
self.structure{1} = north(self.structure{1},NPole);

    function s = rotatSym(nam)
        % determines the s-fold rotational symmetry present in a shape
        % defined by nam
        nam(nam==1)=[];
        G = double(gcd(sym(nam)));
        if G == 1
            s= 1;
        elseif all(nam==G)
            s= max(nam);
        else
            s= max(nam/G);
        end
    end

    function tf = rotatChk(phis, prime, thresh)
        % determines whether the set of coordinates phis are rotated in
        % accordance with whether they have a prime or not
        % returns true if something is amiss
        if length(phis)==1; tf=false; return; end
        if prime
            phis = (phis+reCycle(phis))/2;
        end
        tf = ~any(abs(fixPhi(phis, 0))< thresh);
    end

    function phi = fixPhi(phi,ref)
        %returns values of phi such that the skip at each cycle is as far
        %away from ref as possible
        phi = mod(phi -ref+pi, 2*pi) +ref-pi;
        % phi(nums==1)=0.2;
    end

    function x = reCycle(x)
        x = x(mod(L(x),length(x))+1);
        x(end)=x(end)+2*pi;
    end
end

function [self,counts]= removeDupes(self)
% removes shapes that are geometrically identical
refi = 1;
while refi < length(self)
    counts(refi) = 1;
    chki = refi+1;
    while chki <= length(self)
        if self(refi).isSameAs(self(chki))
            if self(refi).minVal > self(chki).minVal
                self(refi) = self(chki);
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
    if length(self)>1
        other = arrayfun(@(a) self(a).copySelf, 1:length(self));
        return
    end
    other = shape();

    other.name      = self.name;
    other.numVerts  = self.numVerts;
    other.cartCoor  = self.cartCoor;
    other.minVal    = self.minVal;
    other.structure = self.structure;
    other.thresholds= self.thresholds;
    other.ratings   = self.ratings;
end

%%%%%%%%%%%%%%
end%methods
end%classdef