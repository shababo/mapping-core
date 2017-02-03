function detection_results = detect_peaks(traces,threshold,min_window,return_bit_vec,start_ind,threshold_min,do_z,do_hp)

if return_bit_vec
    detection_results = zeros(size(traces));
else
    detection_results = cell(size(traces,1),1);
end

for i = 1:size(traces,1)
    

%     [~, event_times, w, p] = findpeaks(traces(i,:),...
%         'MinPeakHeight', max(threshold_min,threshold*std(traces(i,start_ind:end))),...
%         'MinPeakDistance', min_window,'MinPeakProminence',40,'MaxPeakWidth',200,...
%         'MinPeakWidth',1,'WidthReference','halfheight');
    if do_hp
        filtered_trace = highpass_filter(traces(i,:),20000);
    else
        filtered_trace = traces(i,:);
    end
    if do_z
        scale = std(filtered_trace(1:100));
    else
        scale = 1;
    end
    [~, event_times, w, p] = findpeaks(filtered_trace,...
        'MinPeakHeight',max(threshold_min,threshold*scale),...
        'MinPeakDistance',min_window);%,'MinPeakProminence',5);%,'MaxPeakWidth',10);
%     if ~isempty(w)
%         event_times
%         w
%         p
%     end
    

    event_times(event_times < start_ind) = [];
    
    if return_bit_vec
%         for j = 1:length(event_times)
            detection_results(i,event_times) = 1;
%         end
    else
        detection_results{i} = event_times;
    end
end 
