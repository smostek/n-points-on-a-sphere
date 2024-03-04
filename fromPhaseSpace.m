function C=fromPhaseSpace(P)
    S(3:n,2) = P(1:n-2);
    S(2:n,1) = P(n-1:end);
    C = sphere2cart([ones(n,1) S]);
end