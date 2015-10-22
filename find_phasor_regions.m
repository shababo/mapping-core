function [N, R] = find_phasor_regions(I)
% [N, R] = FIND_PHASOR_REGIONS(IMAGE) 
% Sample function that generates ten randomly placed circles to
% be replaced by function that identifies N regions from image I
% and populates R matrix where R(:,:,n) is the nth region
	
    
    rng(1);
    
	rows = size(I, 1);
	columns = size(I, 2);
    
	


    radius = 5;
    spacing = 10;
    
    pad = 100;
    center = [radius+1 radius+1];

    circle = zeros(rows,columns);
    for i = 1:2*radius
        for j = 1:2*radius
            if norm([i j] - center) <= (radius - 1)
                circle(i,j) = 1;
            end
        end
    end

    circle_ind = find(circle(:));

    %N = (floor((columns - 2*radius + 1)/spacing) + 1)*(floor((rows - 2*radius + 1)/spacing) + 1) + 1;
    N = 25;
    R = zeros(rows, columns, N);
    rois_avail = 1:N;

    i = 1;
    center = [pad+1 pad+1];
    while center(1) < rows - radius - pad
        while center(2) < columns - radius - pad
            
            center_ind = sub2ind([rows columns],center(1),center(2));
            ind_vec = circle_ind + center_ind;
            ind_vec = ind_vec(ind_vec <= columns*rows);

            this_roi = randsample(rois_avail,1);
            %this_roi = i;
            
            tmp_matrix = zeros(rows,columns);
            tmp_matrix(ind_vec) = this_roi;
            R(:,:,this_roi) = tmp_matrix;
            rois_avail = setdiff(rois_avail,this_roi);

            center(2) = center(2) + spacing;
            
            if i >= N
                break
            end
            i = i + 1;
            
        end
        
        if i >= N
                break
        end
            
        
        center(2) = pad + 1;
        center(1) = center(1) + spacing;
    end
    
    
    %size(R);
    %N
    figure; imagesc(sum(R,3));   
end