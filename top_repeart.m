function [topIndices] = top_repeart(idx)
[uniqueIdx, ~, groupIdx] = unique(idx);
counts = accumarray(groupIdx, 1);

% Sort counts in descending order
[sortedCounts, sortedIdx] = sort(counts, 'descend');
topIndices = uniqueIdx(sortedIdx(1:4)); % Extract top 4

end