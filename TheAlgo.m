function [success, fail] = TheAlgo(n)

clc
fprintf('Generating Starting Positions')
s=initPos(n,false);
fprintf(': %d Created\n',length(s))
s.forceMin(1e-8,'disp',false)

fprintf('Removing Duplicates')
minned = arrayfun(@(a)s(a).minVal<1e-5, 1:length(s));
fail = s(~minned);
success=s(minned).removeDupes;
fprintf(': %d Unique Solutions Found\nAsigning Names', length(success))
success.rename(0.1)
fprintf(': Done!\n')
viewManyShapes(success, 'de',{'ba','hi'}, 'ra',[1 2], 'title','name')