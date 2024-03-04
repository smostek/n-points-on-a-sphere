function minAnalyze(name)
arguments
    name char
end
close all

[self,coorf,vv,nums]=shape(name,false);
nxt = @(c,t,f,m) cos(t*m).*c + sin(t*m).*f./m;
C1 = self.cartCoor; n = self.numVerts;
C2 = C1;
minMode=1;
del = pi/4;
slot = @(list,i,val)[list(1:i-1,:); val; list(i+1:end,:)];
NPole = nums(1)==1;
I = 1+NPole;
bsi=arrayfun(@(a) sum(nums(1:a-1))+1, 1:length(nums));
bms=@(k) bsi(k):bsi(k)+nums(k)-1;
B = find(I<[bsi n],1)-1-NPole;
figure('Units','normalized','Position',[0.05 0.05 0.90 0.90])

%%%%%%%%%%%%%%%%%%%
% sphere & slider %
%%%%%%%%%%%%%%%%%%%
sa = axes('Units','normalized','Position',[0.11 0.55 0.20 0.30]);
[x,y,z]=sphere(16);
lims = vv(B)+[-del, del];
hold(sa,'on')
surf(sa,x,y,z,'FaceColor','c','FaceAlpha',0.1,'EdgeColor','c')
sScat = scatter3(sa,C1(:,1),C1(:,2),C1(:,3),50,'filled');
set(sScat,'CData',slot(zeros(n,3),I,[1 0 1]))
p = @(b,k) C1(bms(b),k);
planes = arrayfun(@(b) patch(sa,p(b,1),p(b,2),p(b,3),'b','FaceAlpha',0.5, ...
        'EdgeColor','b'), 1:length(nums) );
sSld = uicontrol('style','slider', 'Units', 'normalized', ...
    'position',[0.43, 0.45, 0.24, 0.03],...
    'min', lims(1), 'max', lims(2), 'Value', vv(B),'Callback',@sSldCB);

view(sa,[0.05 -1 .2])
axis(sa,'equal','off')

%%%%%%%%%%
% labels %
%%%%%%%%%%
S1 = cart2sphere(C1);
push = sphere2cart([1.1 1 1].*S1);
ltxt = text(sa,push(:,1),push(:,2),push(:,3),cellstr(num2str((1:n)')),'FontSize',15, ...
    'Color',[0.2 0.2 0.2]);

%%%%%%%%%%%%%%%%%%%%%
% Automation Button %
%%%%%%%%%%%%%%%%%%%%%
aBttn = uicontrol('Style','pushbutton','Units','normalized','Position', ...
    [0.15 0.90 0.11 0.03],'String','Go to force mode','Callback',@aBttnCB);

%%%%%%%%%%%%%%%%%%%%%
% Grad Min Selector %
%%%%%%%%%%%%%%%%%%%%%
gBttn = uicontrol('Style','radiobutton','String','include gradient descent', ...
    'Units','normalized','Position',[0.15 0.87 0.11 0.03],'Value',1);

%%%%%%%%%%%%%%%
% Zoom Slider %
%%%%%%%%%%%%%%%
zSld = uicontrol('Style','slider','Units','normalized','Position', ...
    [0.12, 0.45, 0.24, 0.03], 'Min',-5, 'Max',1,'Value',log10(del),...
    'Callback',@zSldCB);

%%%%%%%%%%%%%%%%%
% Force Vectors %
%%%%%%%%%%%%%%%%%
[~,T1,MAG1,R1] = CS(C1);
fFunc = @(c,F,i,a) [c(i,a) c(i,a)+F(i,a)];
fLine = arrayfun(@(i)plot3(sa,fFunc(C1,T1,i,1),fFunc(C1,T1,i,2), ...
    fFunc(C1,T1,i,3),'r','LineWidth',2),1:n);

%%%%%%%%%
% Title %
%%%%%%%%%
tText = uicontrol('Style','text','Units','normalized','Position', ...
    [0.48 0.88 0.13 0.07],'String', titleMaker(C1));

%%%%%%%%%%%%
% dropdown %
%%%%%%%%%%%%
dBox = uicontrol('Style','popupmenu','Units','normalized','Position', ...
    [0.05 0.75 0.04 0.02], 'String',cellstr(num2str((1:n)')), ...
    'Value',I, 'callback',@dBoxCB);

%%%%%%%%%%%%%%%%%
% Potental Well %
%%%%%%%%%%%%%%%%%
pa = axes('Units','normalized','Position',[0.08 0.07 0.3 0.35],'Color',[.94 .94 .94]);
xlabel('theta');ylabel('phi')
clim(pa,'manual')
view(pa,[32 37])
hold(pa,'on')
pPlot2 = well(C1,T1,I,0);
pPlot1 = surf(pa,[1 1;1 1]);

%%%%%%%%%%%%%%%%%%%%%%
% Direction Selector %
%%%%%%%%%%%%%%%%%%%%%%
dTbl = uitable('Position',[590 600 112 300],'Data',ones(n,1),'ColumnEditable',true,...
    'CellEditCallback',@dBoxCB);

%%%%%%%%%%%%%%%%%%%%%
% Simulation Button %
%%%%%%%%%%%%%%%%%%%%%
sBttn = uicontrol('Style','pushbutton','Units','normalized','Position', ...
    [0.03 0.90 0.11 0.03],'String','Simulate','Callback',@sBttnCB);

%%%%%%%%%%%%%%%%%
% history plots %
%%%%%%%%%%%%%%%%%
ha1 = axes('Units','normalized','Position',[0.7 0.55 0.20 0.30]);
title(ha1,'History of Nonradial Force (log plot)')
ha2 = axes('Units','normalized','Position',[0.7 0.08 0.20 0.30]);
title(ha2,'History of Radial Force')
hold([ha1 ha2],'on')
hLine1 = arrayfun(@(i) plot(ha1,1,1,'b'),1:n+1);
set(hLine1(end),'Color','g')
arrayfun(@(i)set(hLine1(i),'XData',[],'YData',[],'UserData',mod(i,n+1)),1:n+1)
hLine2 = arrayfun(@(i) plot(ha2,1,1,'b'),1:n+1);
set(hLine2(end),'Color','g')
arrayfun(@(i)set(hLine2(i),'XData',[],'YData',[],'UserData',mod(i,n+1)),1:n+1)


%%%%%%%%%%%%%%%%%%%
% magnitude plots %
%%%%%%%%%%%%%%%%%%%
ma1 = axes('Units','normalized','Position',[0.45 0.55 0.20 0.30]);
title(ma1,'Normalized Change in Nonradial Force')
ma2 = axes('Units','normalized','Position',[0.45 0.08 0.20 0.30]);
title(ma2,'Normalized Change in Radial Force','FontSize',10)
hold([ma1 ma2],'on')
mLine1 = plot(ma1,[0 2],[0 1],'r');
mLine2 = plot(ma1,[vv(B) vv(B)],[-1 1],'k');
mLine3 = plot(ma2,[0 2],[0 1],'r');
mLine4 = plot(ma2,[vv(B) vv(B)],[-1 1],'k');
plot(ma1, [-1 4],[0 0],'k')
plot(ma2, [-1 4],[0 0],'k')
ylim([ma1 ma2], [-1 1])

dBoxCB(dBox)
set(dTbl,'Enable','off')
dcm = datacursormode;
dcm.Enable = 'off';
dcm.UpdateFcn = @dcmCB;
set([ma1 ma2 ha1 ha2],'YGrid','on')


    function sSldCB(src,~)
        val = get(src,'Value');
        if minMode
            C2 = coorf(slot(vv,B,val));
            trgK = B+NPole;
            set(planes(trgK),'XData',C2(bms(trgK),1),'YData',C2(bms(trgK),2), ...
                             'ZData',C2(bms(trgK),3))
            if val ~= vv(B)
                set(zSld,'Enable','off')
            else
                set(zSld,'Enable','on')
            end
        else
            C2 = nxt(C1,val,T1,MAG1);
            axes(sa)
            f = @(a,i) C2(self.structure{4}{a},i);
            delete(planes)
            planes = arrayfun(@(a) patch(sa,f(a,1),f(a,2),f(a,3),'w', ...
                'EdgeColor','k'), 1:size(self.structure{1},1));
            if val ~=0
                set([zSld dTbl],'Enable','off')
            else
                set([zSld dTbl],'Enable','on')
            end
        end
        [~,T2] = CS(C2);
        T2 = dTbl.Data.*T2;
        set(sScat,'XData',C2(:,1),'YData',C2(:,2),'ZData',C2(:,3))
        set([mLine2 mLine4],'XData', [val val]);
        delete(fLine)
        fLine =  arrayfun(@(i)plot3(sa,fFunc(C2,T2,i,1),fFunc(C2,T2,i,2), ...
            fFunc(C2,T2,i,3), 'r','LineWidth',2),1:n);
        delete(pPlot2)
        S1 = cart2sphere(C1);
        push = sphere2cart([1.1 1 1].*S1);
        arrayfun(@(i)set(ltxt(i),'Position',push(i,:)),1:n)
        pPlot2 = well(C2,T2,I,0,S1(I,3));
        set(tText,'String', titleMaker(C2))
    end

    function dBoxCB(~,~)
        set([hLine1(I) hLine2(I)],'Color','b')
        I = get(dBox,'Value');
        set([hLine1(I) hLine2(I)],'Color','r')
        divs = 41;
        if minMode
            vv(B) = get(sSld,'Value');
            C1 = coorf(vv);
            [~,T1,MAG1,R1]=CS(C1);
            B = find(I<[bsi n+1],1)-1-NPole;
            if B==0||B>length(vv)
                set(sSld,'Enable','off')
                set([mLine1 mLine3],'Visible','off')
                return
            end
            set(sSld,'Enable','on')
            lims = vv(B)+[-del, del];
            set(sSld,'Min',lims(1),'Max',lims(2),'Value',vv(B))
            set([mLine2 mLine4],'XData', [vv(B) vv(B)]);
            xlim([ma1 ma2],lims)
            mx = linspace(lims(1),lims(2),divs);
            [MAG3,tmean,RF3,rmean] = predict(@(v) coorf(slot(vv,B,v)),mx);
        else
            C1 = nxt(C1,sSld.Value,T1,MAG1);
            [~,T1,MAG1,R1]=CS(C1);
            T1 = dTbl.Data.*T1;
            set(sSld,'Value',0)
            set([mLine2 mLine4],'XData', [0 0]);
            mx = linspace(-1,4,divs);
            [MAG3,tmean,RF3,rmean] = predict(@(v) nxt(C1,v,T1,MAG1),mx); 
        end
        clr = slot(repmat([0 0 1], [n 1]),I,[1 0 0]);
        ord = 1:n; ord(I)=[]; ord = [ord I];
        delete([mLine1,mLine3])
        %Nonradial Magnitude
        mLine1 = arrayfun(@(i) plot(ma1,mx,MAG3(i,:),'Color',clr(i,:)),ord);
        mLine1(n+1) = plot(ma1,mx,tmean,'g');
        arrayfun(@(i)set(mLine1(i),'UserData',mod(i,n+1)),[ord n+1])

        %Radial Magnitude
        mLine3 = arrayfun(@(i) plot(ma2,mx, RF3(i,:),'Color',clr(i,:)),ord);
        mLine3(n+1) = plot(ma2,mx,rmean,'g');
        arrayfun(@(i)set(mLine3(i),'UserData',mod(i,n+1)),[ord n+1])
        ylim(ma2, [-4 4].*max([rmean(end) rmean(1)]))

        %Potential Well
        set(pPlot2,'Visible','off')
        delete(pPlot1)
        pPlot1 = well(C1,T1./dTbl.Data,I,1);
        set(sScat,'CData',slot(zeros(n,3),I,[1 0 1]))
        set([zSld dTbl],'Enable','on')

        %History
        RF1 = mag(R1);
        plotAppend(hLine1, length(hLine1(1).XData)+1, log10([MAG1;mean(MAG1)]))
        plotAppend(hLine2, length(hLine2(1).XData)+1, [RF1; mean(RF1)])

    end

    function aBttnCB(~,~)
        if minMode
            if gBttn.Value
                [vv,hst] = gradMin(@(v)CS(coorf(v)), vv);
                [~,~,~,~,tdata,rdata] = predict(@(v)coorf(hst(:,v)),1:size(hst,2));
                tdata = tdata(:,2:end); rdata = rdata(:,2:end);
            else
                vv(B)=get(sSld,'Value');
                tdata = 1; rdata=1;
            end
            minMode = 0;
            C1 = coorf(vv);
            self.cartCoor = C1;
            self.getStructure;
            [~,T1,MAG1,R1]=CS(C1);
            xlim([ma1 ma2],[-1 4])
            set(sSld,'Min',-1,'Value',0,'Max',4)
            set(aBttn,'String','Force Min')
            set(dTbl,'Enable','on')
            delete(gBttn)
        else
            self.cartCoor = nxt(C1,sSld.Value,T1,MAG1);
            [tdata,rdata]=self.forceMin('dispVal',false);
            C1 = self.cartCoor;
            self.getStructure
            [~,T1,MAG1,R1]=CS(C1);
            set(aBttn,'Enable','off')
        end
        tdata=tdata(:,1:end-1); rdata=rdata(:,1:end-1);
        l = 1:size(tdata,2);
        plotAppend(hLine1, length(hLine1(1).XData)+l, log10([tdata; mean(tdata,1)]))
        plotAppend(hLine2, length(hLine2(1).XData)+l, [rdata; mean(rdata,1)])
        dBoxCB(dBox)
        sSldCB(sSld)
        set(pPlot2,'Visible','off')
    end

    function zSldCB(src,~)
        del = 10^get(src,'Value');
        if minMode
            if abs(sSld.Value-vv(B)) >= del
                sSld.Value = del*(sign(sSld.Value-vv(B))+1)+vv(B)-1.01*del;
                sSldCB(sSld)
            end
            set(sSld,'Min',vv(B)-del,'Max',vv(B)+del)
            xlim(vv(B)+[-del del])
        end
        % ylim([ma1 ma3],[0,2*del])
        delete(pPlot1)
        pPlot1 = well(C1,T1,I,1);
        % ylim([ma1 ma2 ma3],diff(zlim(pa))*[-0.5 0.5])
        % set([mLine2 mLine4 mLine6],'YData',ylim(ma2))
    end

    function [t_adj,mean1,r_adj,mean2,tm,rm] = predict(coorf, x)
        [cs,~,~,R]=CS(C1);
        R = mag(R);
        tm = cell2mat(arrayfun(@(v) tanMag(coorf(v)),x, 'un',0));
        mean1 = (sum(tm,1)-cs)./cs;
        t_adj = (tm - tanMag(C1))./tanMag(C1);
        rm = cell2mat(arrayfun(@(v) radMag(coorf(v)),x, 'un',0));
        mean2 = (sum(rm,1)-sum(R))./sum(R);
        r_adj = (rm - R)./R;

        function tm = tanMag(C)
            [~,~,tm] = CS(C);
        end
        function rm = radMag(C)
            [~,~,~,r] = CS(C);
            rm = mag(r);
        end
    end

    function w = well(C,T,i,colr,ref)
        L = (-del/2:del/16:del/2);
        other = C([1:i-1 i+1:end], :);
        c1 = C(i,:);
        s1 = cart2sphere(c1);
        if nargin == 5
            s1(3) = mod(s1(3)-ref+pi,2*pi)+ref-pi;
        end
        [theta,phi] = meshgrid(s1(2)+L, s1(3)+L);
        uSet = arrayfun(@(t,p) coulombPotential(other,sphere2cart([1,t,p])), theta, phi);
        if colr
            colr = uSet;
            adjclr = true;
        else
            colr = max(clim)*ones(size(uSet));
            adjclr = false;
        end
        t = T(i,:);
        c2 = nxt(c1,1,t,mag(t));
        s2 = cart2sphere(c2);
        %make sure it projects to the right place in spherical coordinates
        s2(3) = mod(s2(3)-s1(3)+pi,2*pi)+s1(3)-pi;
        w(1) = surf(pa, theta,phi,uSet,colr, 'FaceAlpha',0.5);
        v = [s1(2) s1(3) coulombPotential(C,i);
             s2(2) s2(3) coulombPotential(other,c2)];
        w(2) = scatter3(pa,v(1,1),v(1,2),v(1,3),30,'k','filled');
        w(3) = plot3(pa,v(:,1),v(:,2),v(:,3),'k','LineWidth',3);
        if adjclr; clim(pa,zlim(pa));end
    end

    function txt = dcmCB(~,info)
        i = info.Target.UserData;
        if isempty(i)
            Pos = info.Position;
            txt = ['(' num2str(Pos(1),2) ', ' num2str(Pos(2),2) ', ' num2str(Pos(3),2) ')'];
        elseif i==0
            txt = 'average';
        else
            txt = ['i: ' num2str(i)];
        end
    end

    function plotAppend(h,x,y)
        if ~isempty(x)
            arrayfun(@(i)set(h(i),'XData',[h(i).XData log10(x)],'YData',[h(i).YData y(i,:)]), 1:length(h))
        end
    end

    function sBttnCB(~,~)
        set([sBttn aBttn zSld sSld],'Enable','off')
        set([planes fLine pPlot1 pPlot2 mLine1 mLine3],'Visible','off')
        C = C2;
        [cs3,T3,mag3,R3] = CS(C);
        dt = 1/15;
        while cs3 > 1e-4
            C = C+T3*dt;
            RF = mag(R3);
            plotAppend(hLine1, length(hLine1(1).XData)+1,log10([mag3;mean(mag3)]))
            plotAppend(hLine2, length(hLine2(1).XData)+1,[RF;mean(RF)])
            set(sScat,'XData',C(:,1),'YData',C(:,2),'ZData',C(:,3))
            drawnow limitrate
            [cs3,T3,mag3,R3] = CS(C);
        end
        C1 = C;
        self.cartCoor = C1;
        minMode = 0;
        xlim([ma1 ma2],[-1 4])
        set(sSld,'min',-1,'Value',0,'max',4,'Enable','on')
        self.getStructure;
        dBoxCB(dBox)
        sSldCB(sSld)
        set(ltxt,'Visible','on')
        set(pPlot2,'Visible','off')
    end

    function t = titleMaker(C)
        [cs,T]=CS(C);
        t= {['Cross Sum = ' num2str(cs)], ...
            ['Potential Energy = ' num2str(coulombPotential(C),8)], ...
            ['Tan Sum Mag = ' num2str(mag(sum(T)))], ...
            ['CoM Magnitude = ' num2str(mag(mean(C)))]   };
    end

end