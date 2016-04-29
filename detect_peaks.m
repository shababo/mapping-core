function detection_results = detect_peaks(traces,threshold,min_window,return_bit_vec,start_ind,threshold_min)

if return_bit_vec
    detection_results = zeros(size(traces));
else
    detection_results = cell(size(traces,1),1);
end

for i = 1:size(traces,1)
    
    [~, event_times] = findpeaks(traces(i,:),'MinPeakHeight',max(threshold_min,threshold*std(traces(i,start_ind:end))),'MinPeakDistance',min_window); %,'MinPeakProminence',30
    event_times(event_times < start_ind) = [];
    
    if return_bit_vec
        for j = 1:length(event_times)
            detection_results(i,event_times(j)) = 1;
        end
    else
        detection_results{i} = event_times;
    end
end 