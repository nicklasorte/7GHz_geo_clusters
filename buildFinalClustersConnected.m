function out = buildFinalClustersConnected(cell_data, x_km, n_cluster, cluster_cap, link_km, maxRadius_km)
% buildFinalClustersConnected
%
% Inputs
%   cell_data    : N-by-5 cell {name, lat1, lon1, lat2, lon2}
%   x_km         : radius used for density ranking (surrounding-count)
%   n_cluster    : number of final clusters to build
%   cluster_cap  : max size per cluster
%   link_km      : edge radius for connectivity graph (controls continuity)
%   maxRadius_km : max allowed distance from the cluster center to any member (bounds area)
%
% Output
%   out.T              : table with names + finalClusterID + centerName + midpoints
%   out.centers        : indices of chosen centers
%   out.centerNames    : names of chosen centers
%   out.G              : proximity graph used for continuity
%   out.unassignedIdx  : indices not assigned (if constraints are too tight)

T0 = pairsTableFromCell(cell_data);

% --- Geodesic midpoints (Mapping Toolbox) ---
[latMid, lonMid] = geodesicMidpoints(T0.lat1, T0.lon1, T0.lat2, T0.lon2);

% --- Density ranking (surrounding count within x_km) ---
countWithinKm = countWithinRadius(latMid, lonMid, x_km);
[~, order] = sort(countWithinKm, "descend");

N = height(T0);

% --- Build proximity graph for continuity (edge if within link_km) ---
G = buildProximityGraph(latMid, lonMid, link_km);

% --- Choose centers from density order (skip duplicates) ---
centers = chooseCenters(order, n_cluster);
n_cluster = numel(centers);

% --- Grow connected, bounded clusters ---
finalClusterID = zeros(N,1);
centerIndex    = zeros(N,1);
isCenter       = false(N,1);
assigned       = false(N,1);

clusterMembers = cell(n_cluster,1);

for c = 1:n_cluster
    ctr = centers(c);
    if assigned(ctr)
        continue
    end

    % Start cluster with center
    members = ctr;
    assigned(ctr) = true;
    isCenter(ctr) = true;

    % Frontier = neighbors of current cluster that are unassigned and within maxRadius
    members = growClusterConnected(G, latMid, lonMid, members, assigned, cluster_cap, ctr, maxRadius_km);

    % Commit
    assigned(members) = true;
    finalClusterID(members) = c;
    centerIndex(members) = ctr;
    clusterMembers{c} = members;
end

% --- Attach leftovers while keeping connectivity (must connect via an edge) ---
unassigned = find(~assigned);
if ~isempty(unassigned)
    remainingCap = cluster_cap - cellfun(@numel, clusterMembers);

    for ii = unassigned(:).'
        % Candidate clusters: those with capacity and an edge to at least one member.
        bestC = 0;
        bestD = Inf;

        for c = 1:n_cluster
            if remainingCap(c) <= 0 || isempty(clusterMembers{c}), continue, end

            % Must be adjacent to current cluster to preserve connectedness
            if ~any(ismember(neighbors(G, ii), clusterMembers{c}))
                continue
            end

            ctr = centers(c);
            d = deg2km(distance(latMid(ii), lonMid(ii), latMid(ctr), lonMid(ctr)));
            if d <= maxRadius_km && d < bestD
                bestD = d;
                bestC = c;
            end
        end

        if bestC ~= 0
            finalClusterID(ii) = bestC;
            centerIndex(ii) = centers(bestC);
            clusterMembers{bestC}(end+1) = ii; %#ok<AGROW>
            remainingCap(bestC) = remainingCap(bestC) - 1;
            assigned(ii) = true;
        end
    end
end

unassignedIdx = find(~assigned);

% --- Output table with name propagated ---
T = T0;
T.latMid = latMid;
T.lonMid = lonMid;
T.countWithinKm = countWithinKm;
T.finalClusterID = finalClusterID;
T.isCenter = isCenter;
T.centerIndex = centerIndex;
T.centerName = strings(N,1);
mask = centerIndex > 0;
T.centerName(mask) = T.name(centerIndex(mask));

out = struct;
out.T = T;
out.centers = centers;
out.centerNames = T.name(centers);
out.G = G;
out.unassignedIdx = unassignedIdx;

end

% -------------------------------------------------------------------------
function T = pairsTableFromCell(C)
T = table;
T.name = string(C(:,1));
T.lat1 = cellfun(@(v) v(1), C(:,2));
T.lon1 = cellfun(@(v) v(1), C(:,3));
T.lat2 = cellfun(@(v) v(1), C(:,4));
T.lon2 = cellfun(@(v) v(1), C(:,5));
end

% -------------------------------------------------------------------------
function [latMid, lonMid] = geodesicMidpoints(lat1, lon1, lat2, lon2)
[arcDeg, az] = distance(lat1, lon1, lat2, lon2);        % Mapping Toolbox :contentReference[oaicite:2]{index=2}
[latMid, lonMid] = reckon(lat1, lon1, arcDeg/2, az);
end

% -------------------------------------------------------------------------
function countWithinKm = countWithinRadius(latMid, lonMid, x_km)
N = numel(latMid);
countWithinKm = zeros(N,1);
for i = 1:N
    dKm = deg2km(distance(latMid(i), lonMid(i), latMid, lonMid));
    countWithinKm(i) = sum(dKm <= x_km) - 1;
end
end

% -------------------------------------------------------------------------
function centers = chooseCenters(order, n_cluster)
centers = zeros(0,1);
for t = 1:numel(order)
    centers(end+1,1) = order(t); %#ok<AGROW>
    if numel(centers) == n_cluster, break, end
end
end

% -------------------------------------------------------------------------
function G = buildProximityGraph(lat, lon, link_km)
% Build an undirected graph where i--j if midpoint distance <= link_km.
% For N~few thousand, O(N^2) is acceptable; for larger N, switch to rangesearch/projected.
N = numel(lat);
S = zeros(0,1); T = zeros(0,1);

for i = 1:N
    dKm = deg2km(distance(lat(i), lon(i), lat, lon));
    nbr = find(dKm <= link_km);
    nbr(nbr <= i) = [];           % keep i<j
    if isempty(nbr), continue, end
    S = [S; repmat(i, numel(nbr), 1)]; %#ok<AGROW>
    T = [T; nbr(:)];                   %#ok<AGROW>
end

G = graph(S, T, [], N);                % MATLAB graph :contentReference[oaicite:3]{index=3}
end

% -------------------------------------------------------------------------
function members = growClusterConnected(G, latMid, lonMid, members, assigned, cap, ctr, maxRadius_km)
% Region growing: only add nodes adjacent to current cluster (keeps continuity)
% and within maxRadius_km of the center.

while numel(members) < cap
    % Candidates: neighbors of any current member
    cand = unique(cell2mat(arrayfun(@(u) neighbors(G,u), members, "UniformOutput", false)));
    cand = cand(~assigned(cand));
    if isempty(cand), break, end

    % Enforce bounded radius
    dToCenter = deg2km(distance(latMid(ctr), lonMid(ctr), latMid(cand), lonMid(cand)));
    cand = cand(dToCenter <= maxRadius_km);
    if isempty(cand), break, end

    % Pick next additions in a greedy way: closest to center first
    dToCenter = deg2km(distance(latMid(ctr), lonMid(ctr), latMid(cand), lonMid(cand)));
    [~, ord] = sort(dToCenter, "ascend");

    space = cap - numel(members);
    take = cand(ord(1:min(space, numel(ord))));

    members = [members; take(:)]; %#ok<AGROW>
    assigned(take) = true;        % reserve them immediately to prevent duplicates
end
end
