if ~exist('type_mask', 'var')
    type_mask=1;
end
if type_mask == 1
    names = dir('*cells.mat');
else
    names = dir('*smask_cells.mat');
end
path = pwd;
if start==1
    finalcelllist = [];
end
for i = start:size(names,1)
% for i = 3
    names(i).name
    findjumps_multijumps(pwd, names(i).name);
    load(names(i).name, 'included_cells');
    pass = find(included_cells);
    if ~isempty(pass)
        for j = 1:size(pass,1)
            k = size(finalcelllist,1);
            finalcelllist(k+1,1).name = names(i).name;
            finalcelllist(k+1,1).index = included_cells(pass(j));
        end 
    end        
end
save('finallist.mat', 'finalcelllist');