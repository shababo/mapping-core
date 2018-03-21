load('12_13_17_0/12_13_17_0_data.mat')
cell_locs = get_rowmat_from_structarray(experiment_setup.neurons,'location');
dists = pdist(cell_locs);
figure; histogram(dists)
title('Pairwise Nuclear Distances, 532 cells, 320 x 320 x 90 um')

xlabel('distance (um)')
ylabel('count')

figure; histogram(dists)
title('Pairwise Nuclear Distances, 532 cells, 320 x 320 x 90 um')
xlabel('distance (um)')
ylabel('count')
xlim([0 50])
ylim([0 500])

%%

distsmat = squareform(dists);
distsmat(distsmat == 0) = Inf;
for i = 1:size(distsmat,1)
    this_min = min(distmat
    for j = 1:i
        distsmat(i,j) = Inf;
    end
    
mindists = min(distsmat);
figure; 
subplot(121)
histogram(mindists)
title('Pairwise Nuclear Distances, 532 cells, 320 x 320 x 90 um')
xlabel('distance (um)')
ylabel('count')


%%


z_dists = [];
for i = 1:size(cell_locs,1)
    for j = i+1:size(cell_locs,1)
        
        if abs(cell_locs(i,1) - cell_locs(j,1)) <= 10 && ...
            abs(cell_locs(i,2) - cell_locs(j,2)) <= 10
        
            z_dists = [z_dists abs(cell_locs(i,3) - cell_locs(j,3))];
            
        end
    end
end


%%

figure; histogram(z_dists)
title('Pairwise Axial Nuclear Distances, 532 cells, 320 x 320 x 90 um')

xlabel('distance (um)')
ylabel('count')

%%


z_dists = [];
for i = 1:size(cell_locs,1)
    these_dists = [];
    these_inds = [];
    for j = 1:size(cell_locs,1)
        
        if i ~= j && abs(cell_locs(i,1) - cell_locs(j,1)) <= 10 && ...
            abs(cell_locs(i,2) - cell_locs(j,2)) <= 10
        
            these_dists = [these_dists abs(cell_locs(i,3) - cell_locs(j,3))];
            these_inds = [these_inds j];
            
        end
    end
    [z_dist_this, ind] = min(these_dists);
    if these_inds(ind) > i % remember we skip i == j
        z_dists = [z_dists z_dist_this];
    end
end


%%

subplot(122); histogram(z_dists)
title('Pairwise Axial Nuclear Distances, 532 cells, 320 x 320 x 90 um')

xlabel('distance (um)')
ylabel('count')

