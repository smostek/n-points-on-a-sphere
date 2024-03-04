function [C,nums,phi,bsi,bms,order,flipped] = bandReorder(C, ref, thresh)
if nargin==2
    thresh=0.1;
end

C = north(C, ref);
[z,order] = sort(C(:,3), 'descend');
nums = diff([0 find(diff(z')<-thresh) size(C,1)]);

numVal = @(nam) sum((1:length(nam)).*nam);
flipped = false;
if numVal(nums(end:-1:1)) > numVal(nums)
    flipped = true;
    nums = nums(end:-1:1);
    order=order(end:-1:1);
    C = -C;
end
bsi=arrayfun(@(a) sum(nums(1:a-1))+1, 1:length(nums));
bms=@(k) bsi(k):bsi(k)+nums(k)-1;
C= C(order,:);
phi = mod(atan2(C(:,2), C(:,1)), 2*pi);
for b = 1:length(nums)
    k = bms(b);
    [phi(k),reorder] = sort(phi(k));
    C(k,:)=C(k(reorder),:);
    order(k)=order(k(reorder));
end