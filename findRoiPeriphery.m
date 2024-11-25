function perif_verts = findRoiPeriphery(roi_verts,roi_faces)

% Finds periphery for provided ROI; criterion is an assignment of an edge
% to less than two triangles (in most cases, it's one triangle but for 
% smaller ROIs there is a chance that there are no triangles for a given
% edge; the rest of edges are supposed to have a double triangle 
% assignment, indicating it's not a periphery).
% INPUT:
% roi_verts - vector of vertex IDs belonging to a given ROI only
% roi_faces - triangles (fx3) with vertex IDs, from same ROI
% OUTPUT:
% perif_verts - vector of all peripheral vertices
% NOTE: verts in this method are vrtx_ids, but in getVrtxWghts signify area
% assignment.
% NOTE2: when supplying input, make sure there are only those faces which
% contain ALL existing vertices (i.e. if e.g. 1 vertex in a face is not a 
% part of the vertex list, it shouldn't be included)
% Author: Gleb Bezgin
% Copyright Gleb Bezgin 2024

% 2or1 verts, no faces, or some verts have no triangles (latter dif-t4compd)
% m/b 3rd condition is wrong (is IS, thus rmvd) - in case of e.g. a 'dumbbell'
if length(roi_verts)<3 || isempty(roi_faces) %|| ...
        %size(intersect(roi_faces,roi_verts),1) < size(roi_verts,1)
    perif_verts = roi_verts;
else
    % 0) 'faceless' vertices (if any) are periphery
    perif_verts = roi_verts(~ismember(roi_verts,roi_faces));
    % 1) init mtrx w/edges
    edge_mtrx = zeros(size(roi_verts,1));
    % more importantly, matrix with all faces adjacent w/edge
    face_mtrx = cell(size(roi_verts,1));
    % loop through faces
    for i=1:size(roi_faces,1)
        % sort as a safety measure (m/b not needed)
        abc = sort(roi_faces(i,:));
        a=abc(1); b=abc(2); c=abc(3);
        a_idx = find(roi_verts==a);
        b_idx = find(roi_verts==b);
        c_idx = find(roi_verts==c);
        % fill edge matrix with ones, 3 per face
        edge_mtrx(a_idx,b_idx) = 1;
        edge_mtrx(b_idx,c_idx) = 1;
        edge_mtrx(a_idx,c_idx) = 1;
        % fill face mtrx w/face ids
        face_mtrx{a_idx,b_idx} = [face_mtrx{a_idx,b_idx}, i];
        face_mtrx{b_idx,c_idx} = [face_mtrx{b_idx,c_idx}, i];
        face_mtrx{a_idx,c_idx} = [face_mtrx{a_idx,c_idx}, i];
    end
    % now loop through all edges and see how many have <2 assgnmnt
    edges = find(edge_mtrx);
    [edges_a,edges_b] = ind2sub(size(edge_mtrx), find(edge_mtrx));
    %perif_verts_inds = [];
    a_s = zeros(size(roi_verts)); b_s = zeros(size(roi_verts));
    for i=1:length(edges)
        this_edge = face_mtrx{edges(i)};
        % seemingly universal periphery condition
        if length(unique(this_edge)) < 2
            % add verts at the end of each edge having <2 faces
            %perif_verts_inds = [perif_verts_inds; edges_a(i); edges_b(i)];
            % mtlb suggested prealloc-n, hence this method
            a_s(edges_a(i)) = 1; b_s(edges_b(i)) = 1;
        end
    end
    % bring together sources and targets from the edges
    perif_verts_inds = [find(a_s); find(b_s)];
    % obtain result by converting back to original inds
    perif_verts = [perif_verts; unique(roi_verts(perif_verts_inds))];
end
