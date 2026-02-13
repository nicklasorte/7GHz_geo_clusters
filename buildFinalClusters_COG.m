function out = buildFinalClusters_COG(cell_data, x_km, n_cluster, cluster_cap, maxIter)
tic;
% buildFinalClusters_COG
%
% Goal: clusters whose "center of gravity" is central within each cluster.
% We do this by iteratively:
%   - computing each cluster centroid with meanm (Mapping Toolbox)
%   - choosing the cluster center as the member closest to that centroid (medoid)
%   - reassigning points to nearest centers with a hard-ish capacity (then spill leftovers)
%
% Inputs
%   cell_data   : N-by-5 cell {name, lat1, lon1, lat2, lon2}
%   x_km        : used only to pick initial centers (density ranking)
%   n_cluster   : desired number of clusters
%   cluster_cap : max size per cluster (if infeasible, a few points will overflow)
%   maxIter     : optional, default 5
%
% Output out:
%   out.T            : table with name, midpoints, finalClusterID, centerName
%   out.centers      : center indices (row indices into original data)
%   out.centerNames  : names of centers
%   out.clusterSizes : sizes
%   out.info         : diagnostics

if nargin < 5 || isempty(maxIter), maxIter = 5; end

T0 = pairsTableFromCell(cell_data);
N = height(T0);

% --- Midpoints (geodesic midpoint using Mapping Toolbox) ---
[latMid, lonMid] = geodesicMidpoints(T0.lat1, T0.lon1, T0.lat2, T0.lon2);

% --- Initial center selection by density ranking (surrounding count within x_km) ---
countWithinKm = countWithinRadius(latMid, lonMid, x_km);
[~, order] = sort(countWithinKm, "descend");

n_cluster = min(n_cluster, N);
centers = order(1:n_cluster);

% --- Iterate: update centers using centroid->medoid, then reassign with cap ---
finalClusterID = zeros(N,1);

for it = 1:maxIter
    % Assign to current centers with capacity
    finalClusterID = assignWithCapacity(latMid, lonMid, centers, cluster_cap);

    % Update centers: centroid (meanm) -> medoid closest to centroid
    newCenters = updateCentersByCentroidMedoid(latMid, lonMid, finalClusterID, centers);

    if isequal(newCenters, centers)
        break
    end
    centers = newCenters;
end

% Final assignment with the final centers
finalClusterID = assignWithCapacity(latMid, lonMid, centers, cluster_cap);

% Pack output
clusterSizes = accumarray(finalClusterID(finalClusterID>0), 1, [n_cluster 1], @sum, 0);

T = T0;
T.latMid = latMid;
T.lonMid = lonMid;
T.countWithinKm = countWithinKm;
T.finalClusterID = finalClusterID;

T.isCenter = false(N,1);
T.isCenter(centers) = true;

T.centerIndex = zeros(N,1);
for c = 1:n_cluster
    members = find(finalClusterID == c);
    T.centerIndex(members) = centers(c);
end

T.centerName = strings(N,1);
mask = T.centerIndex > 0;
T.centerName(mask) = T.name(T.centerIndex(mask));

capacity = n_cluster * cluster_cap;
overflow = max(0, N - capacity);

out = struct;
out.T = T;
out.centers = centers(:);
out.centerNames = T.name(centers(:));
out.clusterSizes = clusterSizes;
out.info = struct("N",N,"x_km",x_km,"n_cluster",n_cluster,"cluster_cap",cluster_cap, ...
                  "capacity",capacity,"overflowIfInfeasible",overflow,"iters",it);
toc;

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
[arcDeg, az] = distance(lat1, lon1, lat2, lon2);          % great-circle arc (deg) :contentReference[oaicite:2]{index=2}
[latMid, lonMid] = reckon(lat1, lon1, arcDeg/2, az);      % half-way along great circle
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
function idx = assignWithCapacity(latMid, lonMid, centers, cap)
% Greedy global assignment by increasing distance to centers, respecting cap.
% If infeasible, remaining points go to nearest center (overflow allowed).

N = numel(latMid);
K = numel(centers);

% Precompute distances to centers (km)
D = zeros(N,K);
for c = 1:K
    ctr = centers(c);
    D(:,c) = deg2km(distance(latMid(ctr), lonMid(ctr), latMid, lonMid));
end

% Build all (point,cluster) pairs sorted by distance
[distSorted, lin] = sort(D(:), "ascend");
[pt, cl] = ind2sub([N K], lin);

idx = zeros(N,1);
sz = zeros(K,1);

for t = 1:numel(distSorted)
    i = pt(t);
    c = cl(t);
    if idx(i) ~= 0
        continue
    end
    if sz(c) < cap
        idx(i) = c;
        sz(c) = sz(c) + 1;
    end
end

% Leftovers: nearest center (even if full)
left = find(idx == 0);
if ~isempty(left)
    [~, cNear] = min(D(left,:), [], 2);
    idx(left) = cNear;
end
end

% -------------------------------------------------------------------------
function centersOut = updateCentersByCentroidMedoid(latMid, lonMid, idx, centersIn)
% For each cluster c:
%   centroid = meanm(lat,lon) (Mapping Toolbox geographic mean) :contentReference[oaicite:3]{index=3}
%   new center = member with minimum distance to centroid (medoid)

K = numel(centersIn);
centersOut = centersIn;

for c = 1:K
    members = find(idx == c);
    if isempty(members)
        centersOut(c) = centersIn(c);
        continue
    end

    % Geographic mean (centroid)
    [latC, lonC] = meanm(latMid(members), lonMid(members));

    % Pick member closest to centroid
    dKm = deg2km(distance(latC, lonC, latMid(members), lonMid(members)));
    [~, kmin] = min(dKm);
    centersOut(c) = members(kmin);
end
end
