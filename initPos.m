function allS = initPos(n,min)
if nargin==1
    min = true;
end

p = partitions(n);
p = p(p(:,1)<=2,:); %delete partitions with more than 2 ones
allS = [];

for i = 1:length(p(:,1))
    % Find all the unique ways of arranging the bands defined by the
    % partition
    NPole = p(i,1)>0;
    SPole = p(i,1)>1;
    bandOrders = partDecoder(p(i,2:end), 2:length(p(i,:)));
    B = length(bandOrders);
    bandOrders = unique(perms(bandOrders), 'rows');
    if (SPole||~NPole)
        bandOrders=removeMirrors(bandOrders);
    end

    % convert partitions into names
    names = ones(size(bandOrders,1), B+NPole+SPole);
    names(:,NPole+1:NPole+B) = bandOrders;
    % add each name set to allS
    biny= dec2bin(0:2^(B-1)-1)=='1';
    for j = 1:size(bandOrders,1)
        r=2^(B-1);
        s=shape();
        bi = banditofix(bandOrders(j,:));
        ticks= [biny(:,1:bi-1), false(r,1), biny(:,bi:end)];
        if isequal(names(j,:), names(j,end:-1:1))
            ticks=removeMirrors(ticks);
            if all(diff(bandOrders(j,:))==0)
                ticks=removeMirrors(ticks,~ticks);
            end
            r=size(ticks,1);
        end
        ticks=[false(r,NPole) ticks false(r,SPole)];
        s(r)=shape();
        for si=1:r
            s(si)=shape( nameConstructor(names(j,:),ticks(si,:)),min );
        end
        allS=[allS, s];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function l = partDecoder(part, i)
    if length(part)==1
        l=i(1)*ones(1,part);
    else
        l=[i(1)*ones(1,part(1)), partDecoder(part(2:end),i(2:end))];
    end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function l1 = removeMirrors(l1,l2)
    if nargin==1; l2=l1; end
    l2=l2(:,end:-1:1);
    k=1;
    while k <= length(l1(:,1))
        di=find(all( l1(k,:)==[l2(1:k-1,:);l2(k+1:end,:)] ,2));
        if di>=k; di=di+1;end
        l1(di,:)=[]; l2(di,:)=[];
        k = k+1;
    end
    end
end