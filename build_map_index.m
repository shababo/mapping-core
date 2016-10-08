function map_index = build_map_index(hologram_info)

num_stims = length(hologram_info);

x_positions = sort(unique([hologram_info.x]));
y_positions = sort(unique([hologram_info.y]));

if length(unique(diff(x_positions))) > 1 || ...
   length(unique(diff(y_positions))) > 1

    disp('Not a uniformly spaced grid... gudbye.');
    return
    
end

x_diff = diff(x_positions(1:2));
y_diff = diff(y_positions(1:2));
x_offset = 1 - x_positions(1)/x_diff;
y_offset = 1 - y_positions(1)/y_diff;

map_index = zeros(num_stims,2);

for i = 1:num_stims
    
    map_index(i,1) = hologram_info(i).x/x_diff + x_offset;
    map_index(i,2) = hologram_info(i).y/y_diff + y_offset;
    
end
    