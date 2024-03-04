function [nums,bands,primes,NPole,SPole] = nameInterpreter(txt)
nums=textscan(txt, '%f', 'Delimiter', {' ', '` '});
if length(nums)>1
    error('name must only contain integers, spaces, and ticks')
end
nums=nums{1}';

NPole= nums(1)==1;
SPole= nums(end)==1;
bands= nums((1+NPole):(end-SPole));
B = length(bands);
if any(nums<1) || any(nums~=round(nums))
    error('name must use only integers greater than zero')
end
if any(bands==1)
    error('cannot have 1s in middle of name vector')
end

% get the ticks from name
primes=false(1,B);
if NPole&&txt(2)=='`' || SPole&&txt(end)=='`'
    error('cannot apply rotation to poles')
end
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