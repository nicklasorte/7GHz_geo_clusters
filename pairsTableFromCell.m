function T = pairsTableFromCell(C)
% Robust extraction for cells like {[21.133]} into doubles.
T = table;
T.name = string(C(:,1));
T.lat1 = cellfun(@(v) v(1), C(:,2));
T.lon1 = cellfun(@(v) v(1), C(:,3));
T.lat2 = cellfun(@(v) v(1), C(:,4));
T.lon2 = cellfun(@(v) v(1), C(:,5));
end