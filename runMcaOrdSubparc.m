function [new_assig_nat,stages_vec] = runMcaOrdSubparc(labels,surf)

% Current version using the 'ziplock' trick (via weigh_roi_transit).
% Assumes that transitional ROIs share some neighbourhood with both
% preceding and successive ROIs.
% input: verts (braak labels), faces<-surf (obtain nghbrs from it!)
% i.e. nghbrs = findNghbrs4allVrts(surf_tri)
% Copyright Gleb Bezgin 2024

% init weights vector
vert_wghts = zeros(size(labels));

verts_uniq = nonzeros(unique(labels));

faces = surf.tri;
disp('Obtaining vertex neighbours...')
nghbrs = findNghbrs4allVrts(faces);

    % need to make initial assignment of (full) faces to rois
    faces_roi_map = zeros(size(faces,1),1);
    for i=1:size(faces,1)
        a=faces(i,1); b=faces(i,2); c=faces(i,3);
        a_roi=labels(a); b_roi=labels(b); c_roi=labels(c);
        % if the face fully belongs to an roi
        if a_roi==b_roi && b_roi==c_roi
            faces_roi_map(i) = a_roi;
        end
    end

% 1st: U, last: invU, the rest: sigmoid

% deal with first ordered roi (u-shape)
disp(['Weighing terminus ROI ' num2str(1) '...'])
this_roi = verts_uniq(1);
roi_verts = find(labels==this_roi);
roi_faces = faces(faces_roi_map==this_roi,:);
%[~,~,roi_wghts_1st] = weigh_roi(roi_verts,roi_faces,'max','linear');
%vert_wghts(labels==this_roi) = roi_wghts_1st;
perif_verts = findRoiPeriphery(roi_verts,roi_faces);
strip_lohi = get_transit_strip(this_roi, this_roi+1, perif_verts, labels, nghbrs);
[~,roi_wghts_1st_dstl] = weigh_roi_distal(roi_verts, roi_faces, strip_lohi, nghbrs);
vert_wghts(labels==this_roi) = roi_wghts_1st_dstl;

% deal with last ordered roi (inv u-shape)
disp(['Weighing terminus ROI ' num2str(verts_uniq(end)) '...'])
this_roi = verts_uniq(end);
roi_verts = find(labels==this_roi);
roi_faces = faces(faces_roi_map==this_roi,:);
%[~,roi_wghts_last,~] = weigh_roi(roi_verts,roi_faces,'max','linear');
%vert_wghts(labels==this_roi) = roi_wghts_last;
perif_verts = findRoiPeriphery(roi_verts,roi_faces);
strip_lohi = get_transit_strip(this_roi, this_roi-1, perif_verts, labels, nghbrs);
roi_wghts_last_dstl = weigh_roi_distal(roi_verts, roi_faces, strip_lohi, nghbrs);
vert_wghts(labels==this_roi) = roi_wghts_last_dstl;


% now all the transitional rois
verts_trnst = verts_uniq(2:end-1);
for i=1:length(verts_trnst)
%for i=1
%for i=2
%for i=3
    this_roi = verts_trnst(i);
    disp(['Weighing transitional ROI ' num2str(this_roi) '...'])
    roi_verts = find(labels==this_roi);
    roi_faces = faces(faces_roi_map==this_roi,:);
    perif_verts = findRoiPeriphery(roi_verts,roi_faces);
    strip_lo = get_transit_strip(this_roi, this_roi-1, perif_verts, labels, nghbrs);
    strip_hi = get_transit_strip(this_roi, this_roi+1, perif_verts, labels, nghbrs);
    roi_wghts = weigh_roi_transit(roi_verts, roi_faces, strip_lo, strip_hi, nghbrs);
    vert_wghts(labels==this_roi) = roi_wghts;
end


% convert these vals to natural order by unique occurrences
disp('Obtaining stage vector description...')
lbls_orig = [];
stages_vec = [];
for i=1:max(labels)
    % how it's weighted so far, e.g. -9..0..9 per stage
    lbls_orig = [lbls_orig; (min(vert_wghts(labels==i)):max(vert_wghts(labels==i)) )'];
    stages_vec = [stages_vec; repmat(i, length( min(vert_wghts(labels==i)):max(vert_wghts(labels==i)) ), 1)];
end
% loop through all and label w/nat order
new_assig_nat = zeros(size(vert_wghts));
disp('Converting sub-parcels to natural order...')
for i=1:length(lbls_orig)
    ths_lbl = lbls_orig(i);
    ths_stage = stages_vec(i);
    new_assig_nat( labels==ths_stage & vert_wghts==ths_lbl ) = i;
end

% checking for gaps
assig_gaps = setdiff( 0:max(new_assig_nat), unique(new_assig_nat) )';
assig_uniq = nonzeros(unique(new_assig_nat));
if ~isempty(assig_gaps)
    % because when there is more than 1 gap it doesn't print
    %{
    disp(['Warning: empty sub-parcel(s) at ' ...
        num2str(assig_gaps) '; shifting assignment by ' ...
        num2str(numel(assig_gaps)) '...'])
    %}
    disp(['Warning: empty sub-parcel(s); shifting assignment by ' ...
        num2str(numel(assig_gaps)) '...'])
    stages_vec(assig_gaps) = [];
    for i=1:length(assig_uniq)
        ths_uniq = assig_uniq(i);
        new_assig_nat(new_assig_nat==ths_uniq) = i;
    end
end

disp('Done.')
