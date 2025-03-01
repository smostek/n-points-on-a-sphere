function s = TheAlgo(n)
% Attempts to find every equilibrium distribution of n electrons on a sphere
if mod(n,1)~=0 || n<1
    error('n must be a whole number')
end
if n<3
    if n==1
        s=shape([0;0;1]);
    else
        s=shape([0,0;0,0;1,-1]);
    end
    s.minOrd=inf;
    s.getSymmetry
    s.rate
    s.see
    return
end

clc
tic

fprintf('Generating Starting Positions')
s=initPos(n,false);
N1 = length(s);
fprintf(': %d Created\n',N1)
s.forceMin(12)
s = s(arrayfun(@(a)s(a).minOrd>8, 1:N1)); %remove non-minimized solutions
fprintf('Minimization Success Ratio = %.4f\n', length(s)/N1)
s.getSymmetry
fprintf('Removing Duplicates')
s=s.removeDupes;
fprintf(': %d Unique Solutions Found\n', length(s))
s.rate

toc;
viewManyShapes(s)


function allS = initPos(n,min)
    if nargin==1
        min = true;
    end
    
    p = partitions(n);
    % Each row of p describes a partition. It has n columns where the value
    % of the ith column gives how many of i contribute to the partition.
    % For instance, one row of the partitions of 7 would be
    % 2 1 1 0 0 0 0, meaning 7 = 1+1+2+3

    p = p(p(:,1)<=2,:); %delete partitions with more than 2 ones
    allS = [];
    
    for i = 1:length(p(:,1))
        % Find all the unique ways of arranging the bands defined by the
        % partition
        NPole = p(i,1)>0;
        SPole = p(i,1)>1;
        % 'bands' are values in the partition greater than one; ones become
        % the poles
        bandOrders = partDecoder(p(i,2:end), 2:length(p(i,:)));
        B = length(bandOrders);
        bandOrders = unique(perms(bandOrders), 'rows');
        if (SPole||~NPole)
            bandOrders=removeMirrors(bandOrders);
        end
    
        % convert partitions into names
        names = ones(size(bandOrders,1), B+NPole+SPole);
        names(:,NPole+1:NPole+B) = bandOrders;
        % a given order of the bands can have no more than 2^(B-1) distinct
        % rotations
        biny= dec2bin(0:2^(B-1)-1)=='1';
        % add each name set to allS
        for j = 1:size(bandOrders,1)
            r=2^(B-1);
            s=shape();
            bi = banditofix(bandOrders(j,:));
            ticks= [biny(:,1:bi-1), false(r,1), biny(:,bi:end)];
            if isequal(names(j,:), names(j,end:-1:1))
                % palindromic name numbers allow for more reduction
                ticks=removeMirrors(ticks);
                if all(diff(bandOrders(j,:))==0)
                    % uniform name numbers have even more reduction
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
        % converts a partition from the format given by the partition
        % function to one in which each of the nonzero elements become
        % their own entry in a vector.
        % For example, [2 1 1 0 0 0 0] becomes [1 1 2 3]
        if length(part)==1
            l=i(1)*ones(1,part);
        else
            l=[i(1)*ones(1,part(1)), partDecoder(part(2:end),i(2:end))];
        end
        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function l1 = removeMirrors(l1,l2)
        % removes all rows of array l1 that appear in reverse order in l2
        % preserves palindromes
        if nargin==1
            l2=l1;
        end
        l2=l2(:,end:-1:1);
        k=1;
        while k <= size(l1,1)
            di=find(all( l1(k,:)==[l2(1:k-1,:);l2(k+1:end,:)] ,2));
            if di>=k; di=di+1;end
            l1(di,:)=[]; l2(di,:)=[];
            k = k+1;
        end
        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function nam = nameConstructor(nums,prims)
        if length(nums)==1
            if prims(1)
                nam=[num2str(nums) '`'];
            else
                nam= num2str(nums);
            end
        else
            if prims(1)
                nam = [num2str(nums(1)) '` ' nameConstructor(nums(2:end),prims(2:end))];
            else
                nam = [num2str(nums(1)) ' '  nameConstructor(nums(2:end),prims(2:end))];
            end
        end
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function i = banditofix(bands)
        % determines the index of the band that should be fixed in place, because
        % rotating that band doesn't change the shape geometrically.
        np=bands(1)==1; %north pole
        bands(bands==1)=[];
        if isempty(bands)
            i=1;
            return
        end
        
        twos= bands==2; 
        odds= logical(mod(bands,2));
        bandsCopy = bands;
        % if any of the bands are odd, any of the odd bands can be fixed
        % failing that, if any of the bands are two, any such band can be fixed
        % if all bands are even and greater than two, divide by two until either
        % the above conditions are met
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
        end
    
end

function plist = partitions(total_sum,candidate_set,max_count,fixed_count)
% extracts the list of all partitions of a number as integer sums of a list of candidates
% usage: plist = partitions(total_sum,candidate_set)
% usage: plist = partitions(total_sum,candidate_set,max_count,fixed_count)
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 2
% Release date: 7/15/08
    
    % default for candidate_set
    if (nargin<2) || isempty(candidate_set)
      candidate_set = 1:total_sum;
    end
    
    % how many candidates are there
    N = length(candidate_set);
    
    % error checks
    if any(candidate_set<0)
      error('All members of candidate_set must be >= 0')
    end
    % candidates must be sorted in increasng order
    if any(diff(candidate_set)<0)
      error('Efficiency requires that candidate_set be sorted')
    end
    
    % check for a max_count. do we supply a default?
    if (nargin<3) || isempty(max_count)
      % how high do we need look?
      max_count = floor(total_sum./candidate_set);
    elseif length(max_count)==1
      % if a scalar was provided, then turn it into a vector
      max_count = repmat(max_count,1,N);
    end
    
    % check for a fixed_count
    if (nargin<4) || isempty(fixed_count)
      fixed_count = [];
    elseif (fixed_count<0) || (fixed_count~=round(fixed_count))
      error('fixed_count must be a positive integer if supplied')
    end
    
    % check for degenerate cases
    if isempty(fixed_count)
      if total_sum == 0
        plist = zeros(1,N);
        return
      elseif (N == 0)
        plist = [];
        return
      elseif (N == 1)
        % only one element in the set. can we form
        % total_sum from it as an integer multiple?
        p = total_sum/candidate_set;
        if (p==fix(p)) && (p<=max_count)
          plist = p;
        else
          plist = [];
        end
        return
     end
    else
      % there was a fixed_count supplied
      if (total_sum == 0) && (fixed_count == 0)
        plist = zeros(1,N);
        return
      elseif (N == 0) || (fixed_count <= 0)
        plist = [];
        return
      elseif (N==1)
        % there must be a non-zero fixed_count, since
        % we did not trip the last test. since there
        % is only one candidate in the set, will it work?
        if ((fixed_count*candidate_set) == total_sum) && (fixed_count <= max_count)
          plist = fixed_count;
        else
          plist = [];
        end
        return
      end
    end
    
    % finally, we can do some work. start with the
    % largest element and work backwards
    m = max_count(end);
    % do we need to back off on m?
    c = candidate_set(end);
    m = min([m,floor(total_sum/c),fixed_count]);
    
    plist = zeros(0,N);
    for i = 0:m
      temp = partitions(total_sum - i*c, ...
          candidate_set(1:(end-1)), ...
          max_count(1:(end-1)),fixed_count-i);
      plist = [plist;[temp,repmat(i,size(temp,1),1)]];  %ok
    end
end

end