function roi_wghts = weigh_roi_transit(roi_verts, roi_faces, strip_lo, strip_hi, nghbrs)

% Use the 'ziplock' approach (instead of long post-hoc Dijkstra evaluation)
% to weigh transit ROI N gradually from border with (N-1) to border with
% (N+1). Note that roi_* inputs are just lists, whereas strip_* inputs
% are indexed binary arrays of size (nr_vertices,1), and nghbrs are all
% neightbours as obtained using findNghbrs4allVrts.
% Returns weights from negative near N-1 through 0 to positive near N+1.
% Used by runMcaOrdSubparc
% Author: Gleb Bezgin
% Copyright Gleb Bezgin 2024

% first weigh roi from 1 periph to ..(centre)
roi_wghts = zeros(size(roi_verts));
iter_n = 0;
roi_verts_cp = roi_verts; % undeletable copy

% starts as orig strip, then nghbrs
strip_lo_ngh = find(strip_lo);
strip_hi_ngh = find(strip_hi);

%while ~isempty(roi_verts)
%while ~isempty(roi_faces)
while ~isempty(strip_lo_ngh) && ~isempty(strip_hi_ngh)
    iter_n = iter_n + 1;
    
%{    
    % find periphery
    p_verts = findRoiPeriphery(roi_verts, roi_faces);
    roi_wghts(ismember(roi_verts_cp,p_verts)) = iter_n;
    % remove periph verts and faces
    roi_verts(ismember(roi_verts,p_verts)) = [];
    % this way it removes all triangles assoc w/vrtx; should b per edge
    [fi,~] = ind2sub(size(roi_faces),find(ismember(roi_faces,p_verts)));
    roi_faces(fi,:) = [];
%}    
    lo_and_hi = intersect(strip_lo_ngh, strip_hi_ngh);
    if isempty(lo_and_hi)
        roi_wghts(ismember(roi_verts_cp, strip_lo_ngh)) = -iter_n;
        roi_wghts(ismember(roi_verts_cp, strip_hi_ngh)) = iter_n;
    else
        % non-overlapping parts
        inLoNotHi = setdiff(strip_lo_ngh, strip_hi_ngh);
        inHiNotLo = setdiff(strip_hi_ngh, strip_lo_ngh);
        roi_wghts(ismember(roi_verts_cp, inLoNotHi)) = -iter_n;
        roi_wghts(ismember(roi_verts_cp, inHiNotLo)) = iter_n;
        % indicating it's even
        if ~rem(iter_n,2)
            roi_wghts(ismember(roi_verts_cp, lo_and_hi)) = -iter_n;
        else
            roi_wghts(ismember(roi_verts_cp, lo_and_hi)) = iter_n;
        end
    end
%{    
    % purge the current strips
    %roi_verts(ismember(roi_verts, strip_lo_ngh)) = [];
    %roi_verts(ismember(roi_verts, strip_hi_ngh)) = [];
    roi_verts(ismember(roi_verts, [strip_lo_ngh;strip_hi_ngh])) = [];
    [fi_lo,~] = ind2sub(size(roi_faces),find(ismember(roi_faces, strip_lo_ngh)));
    [fi_hi,~] = ind2sub(size(roi_faces),find(ismember(roi_faces, strip_hi_ngh)));
    %roi_faces(fi_lo,:) = []; roi_faces(fi_hi,:) = [];
    roi_faces([fi_lo;fi_hi],:) = [];
%}    
    roi_verts(ismember(roi_verts, unique([strip_lo_ngh;strip_hi_ngh]))) = [];
    [fi,~] = ind2sub(size(roi_faces), ...
        find(ismember(roi_faces, unique([strip_lo_ngh;strip_hi_ngh]) )));
    roi_faces(fi,:) = [];
    
    % assign strips neighbours
    ngh_all_lo = unique(cat(1, nghbrs{strip_lo_ngh}));
    ngh_all_hi = unique(cat(1, nghbrs{strip_hi_ngh}));
    
    % exclude inds not in roi, and reassign strips
    strip_lo_ngh = intersect(ngh_all_lo, roi_verts);
    strip_hi_ngh = intersect(ngh_all_hi, roi_verts);
end

% 1..N and -1..-N to 0..N-1 and 0..-N+1 (NO!)
%roi_wghts(roi_wghts>0) = roi_wghts(roi_wghts>0)+1;
%roi_wghts(roi_wghts<0) = roi_wghts(roi_wghts<0)-1;

% finally, convert weights to ensure sigmoid-like transition (INV)
roi_wghts(roi_wghts>0) = -roi_wghts(roi_wghts>0)+max(roi_wghts);
roi_wghts(roi_wghts<0) = -roi_wghts(roi_wghts<0)+min(roi_wghts);
