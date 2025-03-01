function [newC, newF] = coorConvert(oldC,ctype,oldF,ftype)
% Converts between cartesian and spherical coordinates, for both points and vectors
arguments
    oldC (3,:) double
    % Either [x;y;z] or [r;theta;phi] - theta is elevation (-pi/2 at north pole,
    % +pi/2 at south pole); phi is azimuth (-pi to pi, 0 is along positive x)

    ctype = 'Cartesian'
    % use to specify the coordinate system of oldC - convert FROM ctype

    oldF (3,:) double = [];
    % vector components in either cartesian [xhat;yhat;zhat] or spherical
    % [thetahat;phihat;rhat] - note that the converted values are the second output

    ftype = 'Cartesian'
    % use to specify the coordinate systme of oldF - convert FROM ftype
end

if strcmp(validatestring(ctype,{'Cartesian','Spherical'}),'Cartesian')
    % convert x,y,z to r,theta,phi
    r = vecnorm(oldC);
    newC = [r; -asin(oldC(3,:)./r); atan2(oldC(2,:),oldC(1,:))];
    S = newC;
else
    % convert r,theta phi to x,y,z
    S = oldC;

    newC = S(1,:).*[cos(oldC(3,:)).*cos(oldC(2,:));
                    sin(oldC(3,:)).*cos(oldC(2,:));
                   -sin(oldC(2,:))];
end

if ~isempty(oldF)
    S = permute(S,[1 3 2]);
    st = sin(S(2,1,:));
    ct = cos(S(2,1,:));
    sp = sin(S(3,1,:));
    cp = cos(S(3,1,:));
    o = zeros(size(st));
    if strcmp(validatestring(ftype,{'Cartesian','Spherical'}),'Cartesian')
        %convert cartesian vectors to spherical
        P = [cp.*ct, sp.*ct, -st;
              -sp,     cp,    o;
             cp.*st, sp.*st,  ct];
    else
        %convert spherical vectors to cartesian
        P = [ct.*cp, -sp, cp.*st;
             ct.*sp,  cp, sp.*st;
              -st,    o,    ct];
    end
    newF = permute(pagemtimes(P,permute(oldF,[1 3 2])), [1 3 2]);
end