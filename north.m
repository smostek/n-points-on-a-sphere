function C = north(C,ref)
% rotates C such that ref goes to the north pole
% C and ref are given in cartesian
arguments
    C   (:,3) double
    ref (1,3) double = C(1,:)
end
ref = ref/mag(ref);
x=ref(1);
y=ref(2);
z=ref(3);
p=sqrt(x^2+y^2);
r=sqrt(p^2+z^2);

if p ==0
    C=C*sign(z);
    return
end

P = [z*x/(p*r) z*y/(p*r) -p/r; -y/p x/p 0; x/r y/r z/r];
% thetahat', phihat', rhat' -> ref goes to [0;0;1]

C = (P*(C'))';