function [roi_wghts,roi_wghtsInv] = weigh_roi_distal(roi_verts, roi_faces, strip_lohi, nghbrs)

% Same as *_transit, but only from one side (for termini, instead of centre-based)
% (more objective than weigh_roi, as it uses same apprch as trnst regions)
% isLo: 1 if lo (for last terminus), 0 if hi (for 1st terminus; should b invrtd)
% (2do: either isLo, or 2 outputs)
% Author: Gleb Bezgin
% Copyright Gleb Bezgin 2024

roi_wghts = zeros(size(roi_verts));
iter_n = 0;
roi_verts_cp = roi_verts; % undeletable copy

% starts as orig strip, then nghbrs
strip_lohi_ngh = find(strip_lohi);

while ~isempty(strip_lohi_ngh)
    iter_n = iter_n + 1;
    roi_wghts(ismember(roi_verts_cp, strip_lohi_ngh)) = iter_n;
    
    roi_verts(ismember(roi_verts, unique(strip_lohi_ngh))) = [];
    [fi,~] = ind2sub(size(roi_faces), ...
        find(ismember(roi_faces, unique(strip_lohi_ngh) )));
    roi_faces(fi,:) = [];
    
    % assign strip's neighbours
    ngh_all_lohi = unique(cat(1, nghbrs{strip_lohi_ngh}));
    % exclude inds not in roi, and reassign strip
    strip_lohi_ngh = intersect(ngh_all_lohi, roi_verts);
end

% 2bused by 1st terminus
roi_wghtsInv = -roi_wghts + max( roi_wghts );
