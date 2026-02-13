function subclusters = buildSubclusterStruct(T, centers)
numClusters = max(T.subclusterID);
subclusters(numClusters,1) = struct( ...
    "id",[], "centerIndex",[], "centerName","", ...
    "memberIndex",[], "memberName", strings(0,1), ...
    "size",[]);

for k = 1:numClusters
    members = find(T.subclusterID == k);
    c = centers(k);

    subclusters(k).id = k;
    subclusters(k).centerIndex = c;
    subclusters(k).centerName  = T.name(c);
    subclusters(k).memberIndex = members;
    subclusters(k).memberName  = T.name(members);
    subclusters(k).size        = numel(members);
end
end

