function [assign, centers] = greedyAssignByRadius(latMid, lonMid, order, x_km)
N = numel(latMid);
tic;
assigned = false(N,1);
subclusterID = zeros(N,1);
isCenter     = false(N,1);
centerIndex  = zeros(N,1);

centers = zeros(0,1);
cid = 0;

for t = 1:N
    c = order(t);
    if assigned(c), continue, end

    dKm = deg2km(distance(latMid(c), lonMid(c), latMid, lonMid));
    members = find(~assigned & (dKm <= x_km));   % only unassigned

    cid = cid + 1;
    assigned(members) = true;

    subclusterID(members) = cid;
    centerIndex(members)  = c;
    isCenter(c) = true;

    centers(end+1,1) = c; %#ok<AGROW>
end

assign = struct("subclusterID",subclusterID, "isCenter",isCenter, "centerIndex",centerIndex);
toc;
end
