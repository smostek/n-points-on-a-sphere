function [s,G,Ccoor] = shapeWithSymmetry(grp,F,boundaries,corners)
% Use fundamental domains to efficiently determine the coordinates of a shape
% given its symmetry
arguments
    grp
    % one- or two-element vector of whole numbers. First number determines 
    % symmetry group as an index of the following list:
    % C_k, C_kv, C_kh, R_k, D_k, D_kd, D_kh,
    % T,   T_h,  T_d,  O,   O_h, I,    I_h
    % for the first 7 groups, a second number is used to determine k
    % So, for example, [5 3] specifies D_3 symmetry, while [13] specifies
    % Icosahedral symmetry

    F (2,:) double = []
    % Spherical coordinates [theta;phi] of points within the boundary of
    % the fundamental domain (FD)

    boundaries double = []
    % Distance of points along the FD boundaries. 
    % Each column corresponds to one set of boundary values If two points are
    % wanted on one boundary line, and no others, then two columns are needed;
    % all the other matrix values would be zero

    corners (1,:) logical = []
    % Determines whether to add points at the corners of the FD
end

if isscalar(grp)
    if grp < 7
        error('axial groups require an index')
    else
        k = subsref([3 3 3 4 4 5 5], struct('type','()','subs',{{grp-7}}));
    end
else
    k = grp(2);
    grp = grp(1);
end
%% Locate Principal Axes
switch(grp)
    case {1,2,3,4,5,6,7}
        ax = [0;0;1];
        instr = [1;1];
    case {8,9,10}
        ax = coorConvert([ones(1,4);atan(1/sqrt(2))*repmat([-1,1],[1,2]);(pi/2)*(0:3)],'sph');
        instr = [1:4; ones(1,4)];
    case {11,12}
        ax = coorConvert([ones(1,6);atan(1/sqrt(2))*repmat([-1,1],[1,3]);(pi/3)*(0:5)],'sph');
        instr = [1,2,3,5,6,6; 1,1,1,1,1,2];
    case {13,14}
        ax = coorConvert([ones(1,12);repmat(atan(1/2)*[-1,1],[1,5]),-pi/2,pi/2; (pi/5)*(0:9),0,0],'sph');
        instr = [1,2,2,2,2,3,11,10,9,9,9,1; 1,1,2,3,4,1,2,1,1,2,3,1];
        GR = (1+sqrt(5))/2; %golden ratio
end
%% Add Corners
Ccoor = [];
if ~isempty(corners)
    switch grp
        case {1,2}
            % C_k or C_kv
            if corners(1)
                % Pn
                Ccoor = ax;
            end
            if corners(2)
                % Ps
                Ccoor = [Ccoor -ax];
            end
        case {3,4}
            % C_kh or R_k
            if corners
                % Pn
                Ccoor = [ax, -ax];
            end
        case {5,7}
            % D_k or D_kh
            if corners(1)
                % Pn
                Ccoor = [ax,-ax];
            end
            if corners(2)
                % A1
                Ccoor = [Ccoor [cos(2*pi*(1:k)/k);sin(2*pi*(1:k)/k);zeros(1,k)]];
            end
            if corners(3)
                % A2
                Ccoor = [Ccoor [cos(pi*(1:2:2*k)/k);sin(pi*(1:2:2*k)/k);zeros(1,k)]];
            end
        case 6
            % D_kd
            if corners(1)
                % Pn
                Ccoor = [ax,-ax];
            end
            if corners(2)
                % A1
                Ccoor = [Ccoor [cos(pi*(1:2*k)/k);sin(pi*(1:2*k)/k);zeros(1,2*k)]];
            end
            if corners(3)
                % A3
                Ccoor = [Ccoor [cos(0.5*pi*(1:2:4*k)/k);sin(0.5*pi*(1:2:4*k)/k);zeros(1,2*k)]];
            end
        case {8,10}
            % T or T_d
            if corners(1)
                % T1
                Ccoor = ax;
            end
            if corners(2)
                % T2
                Ccoor = [Ccoor -ax];
            end
            if corners(3)
                % Ot
                Ccoor = [Ccoor coorConvert([ones(1,6);-pi/2,zeros(1,4),pi/2;0,(pi/4)*(1:2:7),0],'sph')];
            end
        case 9
            % T_h
            if corners(1)
                % T1
                Ccoor = [ax, -ax];
            end
            if corners(2)
                % Ot
                Ccoor = [Ccoor [0,0;0,0;1,-1], [[1,-1,-1,1;1,1,-1,-1]/sqrt(2);zeros(1,4)]];
            end
        case {11,12}
            % O
            if corners(1)
                % Oo
                Ccoor = ax;
            end
            if corners(2)
                % C
                Ccoor = [Ccoor coorConvert([ones(1,8);-pi/2,repmat(atan(1/sqrt(8))*[-1,1],[1,3]),pi/2;0,(pi/3)*(1:6),0],'sph')];
            end
            if corners(3)
                % CO
                Ccoor = [Ccoor coorConvert([ones(1,12); ...
                    repmat(atan(sqrt(2))*[-1,1],[1,3]),zeros(1,6); ...
                    (pi/3)*(1:6),(pi/6)*(1:2:11)],'sph')];
            end
        case {13,14}
            % I or I_h
            if corners(1)
                % I
                Ccoor = ax;
            end
            if corners(2)
                % D
                Ccoor = [Ccoor coorConvert([ones(1,20); ...
                    repelem(atan([-GR^2,GR^2,-1/GR^2,1/GR^2]/2),5); ...
                    repmat([1:2:9,0:2:8]*(pi/5),[1,2])],'sph')];
            end
            if corners(3)
                % ID
                Ccoor = [Ccoor coorConvert([ones(1,30); ...
                    repmat([-1,1]*atan(GR),[1,5]),repmat([-1,1]*atan(1/GR),[1,5]),zeros(1,10); ...
                    pi*(0:9)/5,pi*(1:10)/5, pi*(1:2:19)/10],'sph')];
            end
    end
end

%% Add interior
if isempty(F)
    Fcoor = [];
    if ~isempty(boundaries) && sum(boundaries,'all')~=0
        % (ignoring corners) shape only has points on mirror planes, so
        % make the transMats for the group corresponding to if you removed
        % all the mirror symmetries from the target group (7 is a special
        % case)
        G = transMats(subsref([1 1 1 4 5 4 3 8 8 8 11 11 13 13], struct('type','()','subs',{{grp}})),k);
    end
else
    G = transMats(grp,k);
    Fcoor = reshape(pagemtimes(G,coorConvert([ones(1,size(F,2));F],'sph')),3,[]);
    if ~isempty(boundaries) && sum(boundaries,'all')~=0
        % remove the mirror transformations
        G = G(:,:,1:end/2); 
    end
end

%% Add Boundaries
Bcoor = zeros(3,0);
flag = false;
if ~isempty(boundaries) && sum(boundaries,'all')~=0
    switch grp
        case 2
            if any(boundaries(1,:))
                % V1
                p = nonzeros(boundaries(1,:))';
                l = length(p);
                flag = ~all(inRange(p,-pi/2,pi/2));
                Bcoor = [cos(p);zeros(1,l);-sin(p)];
            end
            if any(boundaries(2,:))
                % V2
                p = nonzeros(boundaries(2,:))';
                l = length(p);
                flag = flag || ~all(inRange(p,-pi/2,pi/2));
                Bcoor = [Bcoor, coorConvert([ones(1,l);p;pi*ones(1,l)/k],'sph')];
            end
        case 3
            % H1
            flag = any(boundaries<0 | boundaries>=2*pi/k);
            Bcoor = [cos(boundaries);sin(boundaries);zeros(size(boundaries))];
        case 6
            if any(boundaries(1,:))
                % V_3
                p = nonzeros(boundaries(1,:))';
                l = length(p);
                flag = ~all(inRange(p,-pi/2,0));
                Bcoor = [cos(p);zeros(1,l);-sin(p)];
            end
            if any(boundaries(2,:))
                % V_4
                p = nonzeros(boundaries(2,:))';
                l = length(p);
                flag = flag || ~all(inRange(p,-pi/2,0));
                Bcoor = [Bcoor, coorConvert([ones(1,l);p;pi*ones(1,l)/k],'sph')];
            end
        case 7
            % D_kh - the most annoying case
            % points on the vertical mirrors need to have C_kh applied to them, while
            % points on the horizontal mirror need to have C_kv applied, so we'll
            % deal with the latter set later
            if any(boundaries(1,:))
                %V3
                p = nonzeros(boundaries(1,:))';
                l = length(p);
                flag = ~all(inRange(p,-pi/2,0));
                Bcoor = [cos(p);zeros(1,l);-sin(p)];
            end
            if any(boundaries(2,:))
                %V4
                p = nonzeros(boundaries(2,:))';
                l = length(p);
                flag = flag || ~all(inRange(p,-pi/2,0));
                Bcoor = [Bcoor, coorConvert([ones(1,l);p;pi*ones(1,l)/k],'sph')];
            end
        case 9
            % M_Th
            flag = any(abs(boundaries)>pi/4);
            Bcoor = coorConvert([ones(size(boundaries));zeros(size(boundaries));boundaries],'sph');
        case 10
            if any(boundaries(1,:))
                % B_T1O
                p = nonzeros(boundaries(1,:))';
                l = length(p);
                flag = flag || ~all(inRange(p,-atan(1/sqrt(2)),0));
                Bcoor = coorConvert([ones(1,l);p;pi/4+asin(tan(p))],'sph');
            end
            if any(boundaries(2,:))
                % B_OT2
                p = nonzeros(boundaries(2,:))';
                l = length(p);
                flag = flag || ~all(inRange(p,0,atan(1/sqrt(2))));
                Bcoor = [Bcoor coorConvert([ones(1,l);p;pi/4-asin(tan(p))],'sph')];
            end
            if any(boundaries(3,:))
                % M_Td
                p = nonzeros(boundaries(3,:))';
                l = length(p);
                flag = flag || any(abs(p)>atan(1/sqrt(2)));
                Bcoor = [Bcoor coorConvert([ones(1,l);p;zeros(1,l)],'sph')];
            end
        case 12
            if any(boundaries(1,:))
                % B_OCO
                p = nonzeros(boundaries(1,:))';
                l = length(p);
                flag = ~all(inRange(p,-atan(1/sqrt(2)),0));
                Bcoor = coorConvert([ones(1,l);p;pi/6+asin(tan(p)/sqrt(2))],'sph');
            end
            if any(boundaries(2,:))
                % B_COC
                p = nonzeros(boundaries(2,:))';
                l = length(p);
                flag = flag || ~all(inRange(p,0,atan(1/sqrt(8))));
                Bcoor = [Bcoor coorConvert([ones(1,l);p;pi/6-asin(tan(p)*sqrt(2))],'sph')];
            end
            if any(boundaries(3,:))
                % M_O
                p = nonzeros(boundaries(3,:))';
                l = length(p);
                flag = ~all(inRange(p,-atan(1/sqrt(2)),atan(1/sqrt(8))));
                Bcoor = [Bcoor coorConvert([ones(1,l);p;zeros(1,l)],'sph')];
            end
        case 14
            if any(boundaries(1,:))
                % B_IID
                p = nonzeros(boundaries(1,:))';
                l = length(p);
                flag = ~all(inRange(p,-atan(1/2),0));
                Bcoor = coorConvert([ones(1,l);p;pi/10+asin(tan(p)/GR)],'sph');
            end
            if any(boundaries(2,:))
                % B_IDD
                p = nonzeros(boundaries(2,:))';
                l = length(p);
                flag = flag || ~all(inRange(p,0,atan(1/(2*GR^2))));
                Bcoor = [Bcoor coorConvert([ones(1,l);p;pi/10-asin(tan(p)*GR)],'sph')];
            end
            if any(boundaries(3,:))
                % M_I
                p = nonzeros(boundaries(3,:))';
                l = length(p);
                flag = flag || ~all(inRange(p,-atan(1/2),atan(1/(2*GR^2))));
                Bcoor = [Bcoor coorConvert([ones(1,l);p;zeros(1,l)],'sph')];
            end
        otherwise
            warning('the selected symmetry group does not have mirror planes')
    end
    if flag
        warning('boundary value outside acceptable range')
    else
        Bcoor = reshape(pagemtimes(G,Bcoor),3,[]);
    end
    if grp==7 && any(boundaries(3,:))
        % fill in the horizontal mirror skipped earlier
        p = nonzeros(boundaries(3,:));
        p = reshape([p;-p]+(2*pi/k)*(0:k-1),1,[]);
        Bcoor = [Bcoor [cos([p,-p]);sin([p,-p]);zeros(size([p,p]))]];
    end
end

s=shape([Fcoor Bcoor Ccoor]);

% helper functions:

    function G = transMats(type,k)
        G = spinmat(ax(:,1),(2*pi/k)*(0:k-1)');
        for i = 2:size(ax,2)
            if i==12
                %unfortunately I symmetry is dumb and needs an exception
                G = cat(3,G,pagemtimes(spinmat(ax(:,[7,11]),(2*pi/k)*[1,2]),G(:,:,1:k)));
            else
                G = cat(3,G,pagemtimes(spinmat(ax(:,instr(1,i)),(2*pi/k)*instr(2,i)),G(:,:,1:k)));
            end
        end
    
        switch type
            case {1,8,11,13} % no mirrors
                sigma = [];
            case {2,7,10,12,14} % vertical mirrors
                if type == 7
                    G = cat(3,G,pagemtimes(G,diag([1,1,-1])));
                end
                sigma = diag([1,-1,1]);
            case {3,9} % horizontal mirrors
                sigma = diag([1,1,-1]);
            case {4,6} % rotoreflections
                G = cat(3,G,pagemtimes(G,diag([1,1,-1])*spinmat(3,pi/k)));
                if type==6
                    sigma = diag([1,-1,1]);
                else
                    sigma = [];
                end
            case 5 % perpendicular rotation
                sigma = diag([1,-1,-1]);
        end
        
        if ~isempty(sigma)
            G = cat(3,G,pagemtimes(G,sigma));
        end
        
    end %transMats
    
    function tf = inRange(n,a,b)
        tf = n>a & n<b; 
    end

end % main function
