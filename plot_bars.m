function plot_bars(bar_limits, trial_length)

%% THIS FUNCTION IS NOT FUNCTIONAL RIGHT NOW
bar_corner_time = trial_length/10/20000;
bar_corner_y = -offset + vert_offset;

plot([bar_corner_time; bar_corner_time], bar_corner_y + [0; bar_limits(2)], '-k',  bar_corner_time + [0; bar_limits(1)], [bar_corner_y; bar_corner_y], '-k', 'LineWidth', 2)
text(bar_corner_time - bar_limits(1)/2,bar_corner_y + bar_limits(2)/2, [num2str(bar_limits(2)) ' pA'], 'HorizontalAlignment','right')
text(bar_corner_time + bar_limits(1)/2,bar_corner_y - bar_limits(2)/2, [num2str(bar_limits(1)*1000) ' ms'], 'HorizontalAlignment','center')
