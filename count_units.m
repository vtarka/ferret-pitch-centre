function nUnits = count_units(unit_list,column)

if ~exist('column','var')
    column = 3;
end

nUnits = 0;
for pen = 1:length(unit_list)

    pen_units = unit_list{pen,column};

    if ~isempty(pen_units)
        nUnits = nUnits + length(unique(pen_units));
    end
        
end
end