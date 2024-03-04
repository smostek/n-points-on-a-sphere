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