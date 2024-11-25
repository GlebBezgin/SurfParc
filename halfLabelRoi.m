function [half_set,nghbr_boundary] = halfLabelRoi(ths_lbl, nghbr_lbl, hlh_wghts, labels, nghbrs)

% Part of ordered subparcels procedure; labels a given roi adjacently to a
% specified neighbour, until the centre of this given roi (thus 'half*').
% INPUT:
% ths_lbl - which roi should be half-labelled (e.g. Braak4)
% nghbr_lbl - neighbouring (e.g. Braak3)
% hlh_wghts - entire map of rois with zeros in centre, higher in perif
% labels - nVerts map to labels (e.g. Braak 1..6)
% nghbrs - to save time, have all nghbrs input (from findNghbrs*)
% OUTPUT: ones for half-roi adjacent to specified nghbr, zeros otherwise.
% Used by: orderedParcelsTransform
% Author: Gleb Bezgin
% Copyright Gleb Bezgin 2024

% For ordered parclns, ths shld b the case (change otherwise)
if abs(ths_lbl - nghbr_lbl) ~= 1
    error('Parcels are not neighbouring!');
end

half_set = zeros(size(labels));
% for dijk (GB 20210311)
nghbr_boundary = zeros(size(labels));

prcl_ths = labels==ths_lbl;
prcl_nghb = labels==nghbr_lbl;

prcl_ths_perif = labels==ths_lbl & hlh_wghts==max(hlh_wghts(prcl_ths));
prcl_ths_perif_inds = find(prcl_ths_perif);
for j=1:length(prcl_ths_perif_inds)
    ths_idx = prcl_ths_perif_inds(j);
    ths_nghbrs = nghbrs{ths_idx};
    ths_nghbrsNghb = ismember(ths_nghbrs, find(prcl_nghb));
    
    %half_set(ths_nghbrs(ths_nghbrsNghb)) = -100; % should be pos (2chng)
    
    % label actual perif next to nghbr
    if sum(ths_nghbrsNghb)
        half_set(prcl_ths_perif_inds(j)) = 1;
        % to ensure dijk mustn't loop thru all nghbr nodes
        nghbr_boundary(ths_nghbrs(ths_nghbrsNghb)) = 1;
    end
        
end

% now label the rest (until centre reached)
labels2 = labels; % changes as half_set gets deeper
labels2(find(half_set)) = nghbr_lbl;
%half_set = labels2; % testing

%iter = 1;
while numel(find( hlh_wghts(find(half_set)) == 0 ) ) == 0
    hlfst_nghbrs_cell = nghbrs(find(half_set));
    tst = [];
    for i=1:length(hlfst_nghbrs_cell)
        tst = [tst; hlfst_nghbrs_cell{i}];
    end
    hlfst_nghbrs = unique(tst);
    %hlfst_nghbrsNghb = ismember(hlfst_nghbrs, find(labels2==nghbr_lbl));
    hlfst_nghbrsThis = ismember(hlfst_nghbrs, find(labels2==ths_lbl));
    %half_set(find(hlfst_nghbrsThis)) = 1;
    half_set(hlfst_nghbrs(hlfst_nghbrsThis)) = 1;
    labels2(find(half_set)) = nghbr_lbl;
    %iter = iter + 1;
    %disp(['Iter ' num2str(iter) ', ' num2str(numel(find( hlh_wghts(find(half_set)) == 0 ) ) ) ' zero-vrts'])
end
