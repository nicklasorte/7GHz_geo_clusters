function saveClustersToExcel(cell_p2p_data, outFinal, filename)
tic;
% saveClustersToExcel
%
% Writes one cluster per Excel sheet, preserving ALL columns from cell_p2p_data.
% Row 1 of cell_p2p_data is treated as a header and is written to every sheet.
% Converts 'missing' cell elements to "" so writecell can export.
% Writes an "Unmatched" sheet for rows whose name isn't found in outFinal.

if nargin < 3 || isempty(filename)
    filename = "clusters.xlsx";
end

% --- Extract header and body ---
header = cell_p2p_data(1,:);      % 1-by-P
body   = cell_p2p_data(2:end,:);  % (N-1)-by-P

% --- Clean 'missing' in body and header (writecell can't handle missing type) ---
header = replaceMissingInCell(header);
body   = replaceMissingInCell(body);

% --- Build lookup: name -> clusterID (normalize) ---
T = outFinal.T;

keyNames = normalizeKey(T.name);
clusterByName = containers.Map(cellstr(keyNames), num2cell(T.finalClusterID));

% --- Assign each row in body to a cluster using name in col 1 ---
N = size(body,1);
rowCluster = zeros(N,1);
missingName = false(N,1);

for i = 1:N
    nm = normalizeKey(body{i,1});
    if strlength(nm) == 0
        missingName(i) = true;
        continue
    end

    k = char(nm);  % containers.Map keys are char
    if isKey(clusterByName, k)
        rowCluster(i) = clusterByName(k);
    else
        missingName(i) = true;
    end
end

% --- Write each cluster sheet (with header) ---
clusterIDs = unique(rowCluster(rowCluster > 0));
clusterIDs = sort(clusterIDs);

for k = 1:numel(clusterIDs)
    cid = clusterIDs(k);
    rows = find(rowCluster == cid);

    sheetName = makeSafeSheetName("Cluster_" + string(cid));

    block = [header; body(rows,:)]; %#ok<AGROW>
    writecell(block, filename, "Sheet", sheetName, "Range", "A1");
end

% --- Unmatched rows sheet (with header) ---
if any(missingName)
    block = [header; body(missingName,:)];
    writecell(block, filename, "Sheet", "Unmatched", "Range", "A1");
end

% --- Summary sheet ---
summary = cell(numel(clusterIDs)+2, 2);
summary(1,:) = {"finalClusterID","numRows"};
for k = 1:numel(clusterIDs)
    cid = clusterIDs(k);
    summary{k+1,1} = cid;
    summary{k+1,2} = sum(rowCluster == cid);
end
summary{end,1} = "Unmatched";
summary{end,2} = sum(missingName);

writecell(summary, filename, "Sheet", "Summary", "Range", "A1");

fprintf("Wrote %d cluster sheets + Summary%s to %s\n", ...
    numel(clusterIDs), ternary(any(missingName), " + Unmatched", ""), filename);
toc;
end

% -------------------------------------------------------------------------
function C = replaceMissingInCell(C)
% Replace MATLAB 'missing' objects inside a cell array with "" (empty string).
for i = 1:numel(C)
    if ismissing(C{i})
        C{i} = "";   % safe for Excel
    end
end
end

% -------------------------------------------------------------------------
function k = normalizeKey(x)
% Normalize name keys for matching:
% - convert to string
% - trim whitespace
% - convert missing -> ""
% - force scalar string
x = string(x);
x = strip(x);
x(ismissing(x)) = "";
k = x;
end

% -------------------------------------------------------------------------
function s = makeSafeSheetName(s)
% Excel sheet name rules: <=31 chars, no : \ / ? * [ ]
s = string(s);
s = regexprep(s, '[:\\/\?\*\[\]]', '_');
s = strip(s);
if strlength(s) > 31
    s = extractBefore(s, 32);
end
if strlength(s) == 0
    s = "Sheet";
end
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
