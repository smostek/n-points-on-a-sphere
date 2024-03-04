function P=toPhaseSpace(C)
    %P is a vector giving phi_3 to phi_n followed by theta_2 to theta_n
    %rotate a well-chosen point to the north pole:
    d = mag(C-permute(C,[3 2 1]));
    [~,top] = min(mean(d));
    C = C([top (1:top-1) (top+1:end)],:); %lowest avg dist is #1
    C = north(C);
    %convert to spherical and remove three zeros
    S = cart2sphere(C(2:end,:));
    [~,top] = min(S(:,2));
    S = S([top (1:top-1) (top+1:end)],[2 3]); %lowest theta is #2
    S(:,2)=mod(S(:,2)-S(1,2),2*pi); %rotate so phi_2=0;
    S(find(pi-S(:,1)<0.01,1),2) = 2*pi; %fix south pole phi value
    P = sortrows(S(:,[2,1])); %sort by phi 
    P = [P(2:end,1);P(:,2)];
end