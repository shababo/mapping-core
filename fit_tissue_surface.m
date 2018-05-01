function [nuc_tissue_depth, plane] = fit_tissue_surface(nuc_locs, grid_spacing)

grid_edges_x1 = min(nuc_locs(:,1)):grid_spacing:max(nuc_locs(:,1));
grid_edges_x2 = min(nuc_locs(:,2)):grid_spacing:max(nuc_locs(:,2));

surface_est = zeros(length(grid_edges_x1),length(grid_edges_x2));

for i = 1:length(grid_edges_x1)
    for j = 1:length(grid_edges_x2)
        
        edge_x1 = grid_edges_x1(i);
        edge_x2 = grid_edges_x2(j);
        z_depths = nuc_locs(nuc_locs(:,1) >= edge_x1 & nuc_locs(:,1) < edge_x1 + grid_spacing & ...
            nuc_locs(:,2) >= edge_x2 & nuc_locs(:,2) < edge_x2 + grid_spacing,3);
        z_depths_sorted = sort(z_depths);
        
        surface_est(i,j) = median(z_depths_sorted(1:min(1,length(z_depths_sorted))));
        
    end
end

plane.x1 = grid_edges_x1' + grid_spacing/2;
plane.x2 = grid_edges_x2' + grid_spacing/2;
plane.z = surface_est;

[X,Y] = meshgrid(grid_edges_x1 + grid_spacing/2, grid_edges_x2 + grid_spacing/2);

plane.fit = [];
A = [Y(:) X(:) ones(length(X(:)),1)];
b = surface_est(:);
A(isnan(b),:) = [];
b(isnan(b)) = [];

plane.fit = inv(A'*A)*A'*b;

res = b - A*plane.fit;
A(abs(res) > 20,:) = [];
b(abs(res) > 20) = [];

plane.fit = inv(A'*A)*A'*b;
Z = plane.fit(1)*Y + plane.fit(2)*X + plane.fit(3);

nuc_tissue_depth = nuc_locs(:,3) - (plane.fit(1)*nuc_locs(:,1) + plane.fit(2)*nuc_locs(:,2) + plane.fit(3));

figure; 
scatter3(nuc_locs(:,1),nuc_locs(:,2),-nuc_locs(:,3),30,'.'); hold on
mesh(Y,X,-Z)