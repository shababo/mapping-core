function events = find_square_events(trace, threshold)

when_on = trace > threshold;


change_points = (find(diff(when_on) ~= 0)) + 1;

length(change_points)
events = [];

event_c = 1;

for i = 1:2:length(change_points)
    
    events(event_c,:) = change_points(i:i+1);
    event_c = event_c+1;
end