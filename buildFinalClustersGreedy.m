function out = buildFinalClustersGreedy(cell_data, x_km, n_cluster, cluster_cap)
% buildFinalClustersGreedy
%
% Inputs
%   cell_data   : N-by-5 cell {name, lat1, lon1, lat2, lon2} (cells contain 1x1 doubles)
%   x_km        : radius used to compute "surrounding count" density ranking
%   n_cluster   : desired number of final clusters
%   cluster_cap : max size of each final cluster
%
% Outputs (struct out)
%   out.T            : per-row table (name propagated) with finalClusterID, centerName, etc.
%   out.centers      : indices of chosen centers (into original rows)
%   out.centerNames  : names of chosen centers
%   out.clusterSizes : sizes of final clusters
%   out.info         : diagnostics (capacity, overflow, etc.)

T0 = pairsTableFromCell(cell_data);

% --- Midpoints (Mapping Toolbox geodesic midpoint) ---
[latMid, lonMid] = geodesicMidpoints(T0.lat1, T0.lon1, T0.lat2, T0.lon2);

% --- Density score: how many other midpoints within x_km ---
countWithinKm = countWithinRadius(latMid, lonMid, x_km);

% Sort: densest first
[~, order] = sort(countWithinKm, "descend");

N = height(T0);
capacity = n_cluster * cluster_cap;

% --- Choose cluster centers from the top of the density list ---
centers = zeros(n_cluster,1);
chosen = false(N,1);
k = 0;
for t = 1:N
    i = order(t);
    if ~chosen(i)
        k = k + 1;
        centers(k) = i;
        chosen(i) = true;
        if k == n_cluster, break, end
    end
end
if k < n_cluster
    centers = centers(1:k);
    n_cluster = k;
end

% --- Grow clusters: center + nearest unassigned until cap ---
finalClusterID = zeros(N,1);
centerIndex    = zeros(N,1);
isCenter       = false(N,1);

assigned = false(N,1);
clusterSizes = zeros(n_cluster,1);

for c = 1:n_cluster
    ctr = centers(c);
    if assigned(ctr)
        % (Shouldn’t happen with chosen logic, but keep it safe.)
        continue
    end

    % Distances from this center to all points
    dKm = deg2km(distance(latMid(ctr), lonMid(ctr), latMid, lonMid));  % Mapping Toolbox :contentReference[oaicite:1]{index=1}
    dKm(assigned) = Inf;

    % Take the nearest points up to cap
    [~, idxSort] = sort(dKm, "ascend");
    take = idxSort(1:min(cluster_cap, sum(isfinite(dKm))));

    % Assign
    assigned(take) = true;
    finalClusterID(take) = c;
    centerIndex(take) = ctr;

    isCenter(ctr) = true;
    clusterSizes(c) = numel(take);
end

% --- Assign any leftovers to nearest cluster with remaining capacity ---
left = find(~assigned);
if ~isempty(left)
    remainingCap = cluster_cap - clusterSizes;

    for ii = left(:).'
        % distances to all centers
        dToCenters = deg2km(distance(latMid(ii), lonMid(ii), latMid(centers), lonMid(centers)));
        % pick nearest center that still has capacity
        [~, ordC] = sort(dToCenters, "ascend");

        placed = false;
        for jj = ordC(:).'
            if remainingCap(jj) > 0
                finalClusterID(ii) = jj;
                centerIndex(ii) = centers(jj);
                remainingCap(jj) = remainingCap(jj) - 1;
                clusterSizes(jj) = clusterSizes(jj) + 1;
                placed = true;
                break
            end
        end

        if ~placed
            % Capacity is impossible: put into nearest cluster (overflow)
            [~, jbest] = min(dToCenters);
            finalClusterID(ii) = jbest;
            centerIndex(ii) = centers(jbest);
            clusterSizes(jbest) = clusterSizes(jbest) + 1;
        end
    end
end

% --- Output table with names propagated ---
T = T0;
T.latMid = latMid;
T.lonMid = lonMid;
T.countWithinKm = countWithinKm;
T.finalClusterID = finalClusterID;
T.isCenter = false(N,1);
T.isCenter(centers) = true;

T.centerIndex = centerIndex;
T.centerName = strings(N,1);
mask = T.centerIndex > 0;
T.centerName(mask) = T.name(T.centerIndex(mask));

% Diagnostics
overflow = max(0, N - capacity);

out = struct;
out.T = T;
out.centers = centers;
out.centerNames = T.name(centers);
out.clusterSizes = clusterSizes;
out.info = struct( ...
    "N", N, ...
    "x_km", x_km, ...
    "n_cluster", n_cluster, ...
    "cluster_cap", cluster_cap, ...
    "capacity", capacity, ...
    "overflowIfInfeasible", overflow);

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
[arcDeg, az] = distance(lat1, lon1, lat2, lon2);          % :contentReference[oaicite:2]{index=2}
[latMid, lonMid] = reckon(lat1, lon1, arcDeg/2, az);      % :contentReference[oaicite:3]{index=3}
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
