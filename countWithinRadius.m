function countWithinKm = countWithinRadius(latMid, lonMid, x_km)
N = numel(latMid);
countWithinKm = zeros(N,1);
tic;
for i = 1:N
    dKm = deg2km(distance(latMid(i), lonMid(i), latMid, lonMid));
    countWithinKm(i) = sum(dKm <= x_km) - 1;   % exclude self
end
toc;
end