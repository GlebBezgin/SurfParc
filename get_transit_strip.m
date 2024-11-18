function ths_strip = get_transit_strip(ths_lbl, nghbr_lbl, prcl_ths_perif_inds, labels, nghbrs)

% Part of ordered subparcels procedure; labels a given roi adjacently to a
% specified neighbour (just at the periphery)
% INPUT:
% ths_lbl - which roi should be half-labelled (e.g. Braak4)
% nghbr_lbl - neighbouring (e.g. Braak3)
% hlh_wghts - entire map of rois with zeros in centre, higher in perif
% labels - nVerts map to labels (e.g. Braak 1..6)
% nghbrs - to save time, have all nghbrs input (from findNghbrs*)
% OUTPUT: ones for half-roi adjacent to specified nghbr, zeros otherwise.
% Used by: orderedParcelsTransform

% For ordered parclns, ths shld b the case (change otherwise)
if abs(ths_lbl - nghbr_lbl) ~= 1
    error('Parcels are not neighbouring!');
end

ths_strip = zeros(size(labels));

prcl_ths = labels==ths_lbl;
prcl_nghb = labels==nghbr_lbl;

%prcl_ths_perif = labels==ths_lbl & hlh_wghts==max(hlh_wghts(prcl_ths));
%prcl_ths_perif_inds = find(prcl_ths_perif);
for j=1:length(prcl_ths_perif_inds)
    ths_idx = prcl_ths_perif_inds(j);
    ths_nghbrs = nghbrs{ths_idx};
    ths_nghbrsNghb = ismember(ths_nghbrs, find(prcl_nghb));
    
    %half_set(ths_nghbrs(ths_nghbrsNghb)) = -100; % should be pos (2chng)
    
    % label actual perif next to nghbr
    if sum(ths_nghbrsNghb)
        ths_strip(prcl_ths_perif_inds(j)) = 1;

    end
        
end