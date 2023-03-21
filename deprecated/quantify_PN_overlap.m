
overlap_ctr = 0;

for ap = 1:length(PN_units)

    allUnits_old = PN_units{ap,3};
    allUnits_new = new_PN_units{ap,3};

    in_common = intersect(allUnits_old, allUnits_new);

    if ~isempty(in_common)
        sprintf('Animal %s, Pen %s, Unit(s): %s',PN_units{ap,1},PN_units{ap,2},num2str(in_common'))
        overlap_ctr = overlap_ctr + length(in_common);
    end

end