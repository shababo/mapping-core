function draw_regions( regions )
%DRAW_REGIONS Summary of this function goes here
%   Detailed explanation goes here

n = size(regions,3);

reg_image = sum(bsxfun(@times,regions,reshape(1:25,[1 1 25])),3);

figure;
imagesc(reg_image);
hold on;



offset = 0;

for i = 1:n
    
    
    location = find(regions(:,:,i));
    if isempty(location)
        offset = 1;
        continue
    end
    
    [x,y] = ind2sub([256 256],location(1));
    text(y+4,x,num2str(i-offset));


end



