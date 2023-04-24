function nUnits = count_units(unit_list)

nUnits = 0;
for pen = 1:length(unit_list)

    pen_units = unit_list{pen,3};

    if ~isempty(pen_units)
        nUnits = nUnits + length(unique(pen_units));
    end
        
end
end