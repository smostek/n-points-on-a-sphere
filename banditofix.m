function i = banditofix(bands)
np=bands(1)==1;
bands(bands==1)=[];
if isempty(bands); i=1; return; end

twos= bands==2; 
odds= logical(mod(bands,2));
bandsCopy = bands;
while ~any(odds) && ~any(twos)
    bandsCopy=bandsCopy/2;
    twos= bandsCopy==2;
    odds= logical(mod(bandsCopy,2));
end
if any(odds)
    i = find(odds,1);
else
    i = find(twos,1);
end
i=i+np;