function symmetryContour(symGrp,F,boundaries,corners)
% see how the electrostatic potential is affected by symmetry

f = size(F,2);
b = length(nonzeros(boundaries));
c = sum(corners);

% set up axes:
follow = false;
if strcmp(gcf().UserData,'contour') % reusing existing figure
    clf
elseif ~isempty(gcf().UserData) % current fiugre is a non-contour figure
    figure
end
set(gcf,'Color','w','WindowButtonDownFcn',@changeFollowState, ...
    'WindowButtonMotionFcn',@movePoint,'UserData','contour')
ax = gca;
cla(ax)
set(ax,'Units','normalized','Position',[0.2,0,0.6,0.85])
[x,y,z] = sphere;
surf(ax,x,y,z,'FaceColor',[.7,.7,.7],'EdgeColor','none','HitTest','off');
set(ax,'XColor','w','YColor','w','ZColor','w','Color','w','Units','pixels')
axis(ax,'equal')
grid(ax,'off')
hold(ax,'on')
view(ax,[110,10])

% Draw the fundamental Domain
function h=Easyplot(S,colr,style)
    Cart = coorConvert(S,'sph');
    h=plot3(Cart(1,:),Cart(2,:),Cart(3,:),'Color',colr,'LineStyle',style,'LineWidth',2.5);
end
blu = [0.00,0.45,0.74];
yel = [0.93,0.69,0.13]; 
red = [0.64,0.08,0.18]; 
switch symGrp(1)
    case 1 %C_k
        t = linspace(-pi/2,pi/2,30);
        Easyplot([ones(size(t));t;zeros(size(t))],blu,'-');
        Easyplot([ones(size(t));t;2*pi*ones(size(t))/3],blu,':');
    case 2 %C_kv
        t = linspace(-pi/2,pi/2,30);
        Easyplot([ones(size(t));t;zeros(size(t))],blu,'--');
        Easyplot([ones(size(t));t;pi*ones(size(t))/3],yel,'--');
    case 3 %C_kh
        t = linspace(-pi/2,0,16);
        Easyplot([ones(size(t));t;zeros(size(t))],blu,'-');
        Easyplot([ones(size(t));t;2*pi*ones(size(t))/3],blu,':');
        p = linspace(0,2*pi/3,25);
        Easyplot([ones(size(p));zeros(size(p));p],red,'--');
    case 4 %R_k
        t = linspace(-pi/2,pi/2,30);
        Easyplot([ones(size(t));t;zeros(size(t))],blu,'-');
        Easyplot([ones(size(t));t;pi*ones(size(t))/3],blu,':');
    case 5 %D_k
        t = linspace(-pi/2,0,16);
        Easyplot([ones(size(t));t;zeros(size(t))],blu,'-');
        Easyplot([ones(size(t));t;pi*ones(size(t))/3],yel,'-');
        t = linspace(0,pi/2,16);
        Easyplot([ones(size(t));t;zeros(size(t))],blu,':');
        Easyplot([ones(size(t));t;pi*ones(size(t))/3],yel,':');
    case 6 %D_kd
        t = linspace(-pi/2,0,16);
        Easyplot([ones(size(t));t;zeros(size(t))],blu,'--');
        Easyplot([ones(size(t));t;pi*ones(size(t))/3],yel,'--');
        p = linspace(0,pi/6,7);
        Easyplot([ones(size(p));zeros(size(p));p],red,'-');
        p = linspace(pi/6,pi/3,7);
        Easyplot([ones(size(p));zeros(size(p));p],red,':');
    case 7 %D_kh
        t = linspace(-pi/2,0,16);
        Easyplot([ones(size(t));t;zeros(size(t))],blu,'--');
        Easyplot([ones(size(t));t;pi*ones(size(t))/3],yel,'--');
        p = linspace(0,pi/3,13);
        Easyplot([ones(size(p));zeros(size(p));p],red,'--');
    case 8 %T
        p = linspace(0,pi/4,15);
        Easyplot([ones(size(p));atan(-sin(pi/4-p));p],blu,'-');
        Easyplot([ones(size(p));atan(+sin(pi/4-p));p],yel,'-');
        p = linspace(-pi/4,0,15);
        Easyplot([ones(size(p));atan(-sin(pi/4+p));p],blu,':');
        Easyplot([ones(size(p));atan(+sin(pi/4+p));p],yel,':');
    case 9 %T_h
        p = linspace(0,pi/4,15);
        Easyplot([ones(size(p));atan(-sin(pi/4-p));p],blu,'-');
        p = linspace(-pi/4,0,15);
        Easyplot([ones(size(p));atan(-sin(pi/4+p));p],blu,':');
        p = linspace(-1,1,15)*pi/4;
        Easyplot([ones(size(p));zeros(size(p));p],red,'--');
    case 10 %T_d
        p = linspace(0,pi/4,15);
        Easyplot([ones(size(p));atan(-sin(pi/4-p));p],blu,'--');
        Easyplot([ones(size(p));atan(+sin(pi/4-p));p],yel,'--');
        t = linspace(-1,1,15)*atan(1/sqrt(2));
        Easyplot([ones(size(t));t;zeros(size(t))],red,'--');
    case 11 %O
        p = linspace(0,pi/6,15);
        Easyplot([ones(size(p));atan(-sin(pi/6-p))*sqrt(2);p],blu,'-');
        Easyplot([ones(size(p));atan(+sin(pi/6-p))/sqrt(2);p],yel,'-');
        p = linspace(-pi/6,0,15);
        Easyplot([ones(size(p));atan(-sin(pi/6+p))*sqrt(2);p],blu,':');
        Easyplot([ones(size(p));atan(+sin(pi/6+p))/sqrt(2);p],yel,':');
    case 12 %O_h
        p = linspace(0,pi/6,15);
        Easyplot([ones(size(p));atan(-sin(pi/6-p))*sqrt(2);p],blu,'-');
        Easyplot([ones(size(p));atan(+sin(pi/6-p))/sqrt(2);p],yel,'-');
        t = linspace(-atan(1/sqrt(2)),atan(1/sqrt(8)),15);
        Easyplot([ones(size(t));t;zeros(size(t))],red,'--');
    case 13 %I
        GR = (1+sqrt(5))/2; % golden ratio
        p = linspace(0,pi/10,15);
        Easyplot([ones(size(p));atan(-sin(pi/10-p))*GR;p],blu,'-');
        Easyplot([ones(size(p));atan(+sin(pi/10-p))/GR;p],yel,'-');
        p = linspace(-pi/10,0,15);
        Easyplot([ones(size(p));atan(-sin(pi/10+p))*GR;p],blu,':');
        Easyplot([ones(size(p));atan(+sin(pi/10+p))/GR;p],yel,':');
    case 14 %I_h
        GR = (1+sqrt(5))/2; % golden ratio
        p = linspace(0,pi/10,15);
        Easyplot([ones(size(p));atan(-sin(pi/10-p))*GR;p],blu,'-');
        Easyplot([ones(size(p));atan(+sin(pi/10-p))/GR;p],yel,'-');
        t = linspace(-atan(1/2),atan(0.5/GR^2),15);
        Easyplot([ones(size(t));t;zeros(size(t))],red,'--');
end

[~,G,Ccoor] = shapeWithSymmetry(symGrp,F,boundaries,corners);
Fcoor = reshape(pagemtimes(G,coorConvert([ones(1,size(F,2));F],'sph')),3,[]);
C = [Fcoor(:,1:f),Ccoor(:,1:c),Fcoor(:,f+1:end),Ccoor(:,c+1:end)];
n = size(C,2);

% display the cross sum:
ttl = title(ax,['Minimization Order: ' num2str(-log10(CS(C)))],' ');

% use the model of how U is distributed to determine the edges to use for
% the contour plot:
pd = makedist('GeneralizedExtremeValue','k',0.4,'sigma',n/8+1/2,'mu',0.9*n-0.7);
numLvls = 40;
edges = mean(pd) + std(pd)*linspace(-0.75,0.75,numLvls);
hCont = arrayfun(@(i)plot3(nan,nan,nan,'k'),1:numLvls);

% use these for sampling the contour plot:
res=75;
th=asin((1/res)*(2*(1:res)-(res+1)));
ph=linspace(0,2*pi,res);
[theta,phi] = meshgrid(th,ph);

scat = scatter3(ax,C(1,:),C(2,:),C(3,:),[],repelem([1;10],[f+b+c,n-(f+b+c)]),'filled');
I = 1; % index of the point being focused on
drawStuff;

function drawStuff
    set(scat,'XData',C(1,:),'YData',C(2,:),'ZData',C(3,:))
    set(ttl,'String',['Minimization Order: ' num2str(-log10(CS(C)))])

    U = arrayfun(@(t,p) coulombPotential(C(:,[1:I-1,I+1:n]),coorConvert([1;t;p],'sph')),theta,phi);
    M = contourc(th,ph,U,edges);
    i=1;
    val = M(1,1);
    j = find(round(edges,2)==round(val,2));
    set(hCont(1:j-1),'Visible','off')
    S=[];
    M = [M, [inf;0]];
    while i<=length(M)
        Ni=M(2,i);
        if M(1,i)-val > 0.05 % if we've reached the end of this set that share values
            Q=coorConvert(S,'sph');
            set(hCont(j),'XData',Q(1,:),'YData',Q(2,:),'ZData',Q(3,:),'Visible','on');
            S=[ones(1,Ni);M(:,i+1:i+Ni)];
            val = M(1,i);
            j=j+1;
        else
            S = [S,nan(3,1),[ones(1,Ni);M(:,i+1:i+Ni)]];
        end
        i = i+Ni+1;
    end
end

function movePoint(~,~)
    if follow && I <= f+b %don't move corner points
        r_ao = ax.CurrentPoint(1,:)';
        r_ba = ax.CurrentPoint(2,:)'-r_ao;
        d = vecnorm(r_ba);
        k = (1/d^2)*(-r_ao'*r_ba-sqrt((r_ao'*r_ba)^2-d^2*(r_ao'*r_ao-1)));
        if isreal(k)
            S = coorConvert(r_ao + k.*r_ba); %intersection point in spherical coors
            F(:,I) = S([2,3]);
            Fcoor = reshape(pagemtimes(G,coorConvert([ones(1,size(F,2));F],'sph')),3,[]);
            C = [Fcoor(:,1:f),Ccoor(:,1:c),Fcoor(:,f+1:end),Ccoor(:,c+1:end)];
            drawStuff;
            drawnow
        end
        % pause(0.05)
    end
end

function changeFollowState(~,~)
    if follow
        follow = false;
    else
        r_ao = ax.CurrentPoint(1,:)';
        r_ba = ax.CurrentPoint(2,:)'-r_ao;
        d = vecnorm(r_ba);
        k = (1/d^2)*(-r_ao'*r_ba-sqrt((r_ao'*r_ba)^2-d^2*(r_ao'*r_ao-1)));
        if isreal(k)
            r_po = r_ao + k.*r_ba;
            I_test = find(vecnorm(C-r_po)<=0.05,1);
            if isempty(I_test) || I > f+b+c
                follow = false;
            else
                I = I_test;
                follow = I <= f+b; % don't follow corners
                drawStuff;
                drawnow
            end
        end
    end
end

end