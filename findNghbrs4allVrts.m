function nghbrs = findNghbrs4allVrts(surf_tri)

% just a quick convenience func, neighbours of each vertex (using faces)
% Used by: orderedParcelsTransform
% Copyright Gleb Bezgin 2024

nghbrs = cell(max(max(surf_tri)),1);
for i=1:size(surf_tri,1)
    ths_a = surf_tri(i,1);
    ths_b = surf_tri(i,2);
    ths_c = surf_tri(i,3);
    nghbrs{ths_a} = [nghbrs{ths_a}; ths_b; ths_c];
    nghbrs{ths_b} = [nghbrs{ths_b}; ths_a; ths_c];
    nghbrs{ths_c} = [nghbrs{ths_c}; ths_a; ths_b];
end
% make unique
for i=1:length(nghbrs)
    nghbrs{i} = unique(nghbrs{i});
end
