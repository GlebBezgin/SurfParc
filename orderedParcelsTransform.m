function [new_assig_nat,stages_vec] = orderedParcelsTransform(orig_assig,labels,surf)

% For a {max - max-1 - .. - 0 - .. max-1 - max} ordered parcels, where
% order is specified in labels, and using the surf.tri faces data, make
% ordered subparcels {0 - 1 - max1 - max1+1 - .. - maxLast}
% INPUT:
% orig_assig - assignments with highest val at periphery, zero at centre
% labels - ordered regions 1..N (e.g. Braak stages or hierarchy of sorts)
% surf - surface in surfstat format (2change?), where .tri is needed
% OUTPUT:
% new_assig_nat - ordered subparcels; each stage's zero is inflec pt.
% stages_vec - for each unique natural sub-parcel, gives orig stage assgmt
% =====
% orig_assig = tst_vert_wghtsINV;
% labels = BRAAK.braak_all_corr3;
% surf = surfRw;
% DEPENDENCIES:
% findNghbrs4allVrts - map each vertex to the list of its neighbours
% halfLabelRoi - label roi by criterion of its adjacency to another roi
% Copyright Gleb Bezgin 2024

% first and last ordered parcels
new_assig = zeros(size(orig_assig));
new_assig(labels==1) = orig_assig(labels==1);
new_assig(labels==max(labels)) = -orig_assig(labels==max(labels));

% get neighbour assignments
disp('Obtaining vertex neighbourhood relationships...')
nghbrs = findNghbrs4allVrts(surf.tri);

intrmd_steps = 2:max(labels)-1;
for i=1:numel(intrmd_steps)
    disp(['Half-labelling ordered parcel ' num2str( intrmd_steps(i) ) '...']);
    hlfst_lo = halfLabelRoi(intrmd_steps(i), intrmd_steps(i)-1, orig_assig, labels, nghbrs);
    hlfst_hi = halfLabelRoi(intrmd_steps(i), intrmd_steps(i)+1, orig_assig, labels, nghbrs);
    %hlfst_lo = halfLabelRoi(intrmd_steps(i), intrmd_steps(i)-1, orig_assig, labels, nghbrs);
    % find overlap between lo and hi; put zeros there
    hlfst_vrlp_inds = intersect(find(hlfst_lo), find(hlfst_hi));
    %lbl_inds = find(label==intrmd_steps(i));
    % without overlapping part
    hlfst_lo_corr_inds = setdiff(find(hlfst_lo), hlfst_vrlp_inds);
    hlfst_hi_corr_inds = setdiff(find(hlfst_hi), hlfst_vrlp_inds);
    % do assignments: '-' to lo, '+' to hi, 0 in between
    new_assig(hlfst_lo_corr_inds) = -orig_assig(hlfst_lo_corr_inds);
    new_assig(hlfst_hi_corr_inds) = orig_assig(hlfst_hi_corr_inds);
    new_assig(hlfst_vrlp_inds) = 0;
end

% convert these vals to natural order by unique occurrences
lbls_orig = [];
stages_vec = [];
for i=1:max(labels)
    % how it's weighted so far, e.g. -9..0..9 per stage
    lbls_orig = [lbls_orig; (min(new_assig(labels==i)):max(new_assig(labels==i)) )'];
    stages_vec = [stages_vec; repmat(i, length( min(new_assig(labels==i)):max(new_assig(labels==i)) ), 1)];
end
% loop through all and label w/nat order
new_assig_nat = zeros(size(new_assig));
for i=1:length(lbls_orig)
    ths_lbl = lbls_orig(i);
    ths_stage = stages_vec(i);
    new_assig_nat( labels==ths_stage & new_assig==ths_lbl ) = i;
    %ths_stage_vec = labels==ths_stage;
    %ths_stage_assig = new_assig(ths_stage_vec);
    %ths_stage_assig = ths_stage_assig(ths_stage_vec);
    %new_stage_assig(ths_stage_assig==ths_lbl) = i;
    %(ths_stage_assig==ths_lbl);
end

%{
% the rest
intrmd_steps = 2:max(labels)-1;
for i=intrmd_steps
    prcl_ths = labels==i;
    prcl_ths_perif = labels==i & orig_assig==max(orig_assig(prcl_ths));
    prcl_prev = labels==i-1;
    prcl_next = labels==i+1;
    % go through the neighbours of perif & see those bordering i-1 & i+1
    prcl_ths_perif_inds = find(prcl_ths_perif);
    for j=1:length(prcl_ths_perif_inds)
        ths_idx = prcl_ths_perif_inds(j);
        ths_nghbrs = nghbrs{ths_idx};
        ths_nghbrsPrev = ismember(ths_nghbrs, find(prcl_prev));
        ths_nghbrsNext = ismember(ths_nghbrs, find(prcl_next));
        new_assig(ths_nghbrs(ths_nghbrsPrev)) = -100; % should be pos (2chng)
        new_assig(ths_nghbrs(ths_nghbrsNext)) = 100; % should be neg
        % now repeat procedure until inflec pt (0, centre)
    end
end
%}
