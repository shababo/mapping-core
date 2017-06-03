function [target_locs, stim_inds, stim_pow] = build_stim_locs(trial_metadata,z_inds)

stim_locs = [];
count = 1;
stim_locs_mat = [];
stim_pow = [];
for i = 1:length(trial_metadata)
    
    this_seq = trial_metadata(i).sequence;
    this_stim_key = trial_metadata(i).stim_key;
    for j = 1:length(this_seq)
        stim_pow = [stim_pow; this_seq(j).target_power];
        for k = 1:size(this_stim_key,3)
            stim_locs(count,:,k) = this_stim_key(this_seq(j).precomputed_target_index,:,k);
            if ~any(isnan(stim_locs(count,:,k)))
                stim_locs(count,3,k) = z_inds(i);
                stim_locs_mat = [stim_locs_mat; stim_locs(count,:,k)];
            else
                stim_locs(count,3,k) = NaN;
            end
        end
        count = count + 1;
    end
end
size(stim_locs);
size(stim_locs_mat)
target_locs = unique(stim_locs_mat,'rows');
% size(target_locs)
% assignin('base','target_locs',target_locs)
% assignin('base','stim_locs',stim_locs)
stim_inds = zeros(size(stim_locs,1),size(stim_locs,3));
for i = 1:size(stim_locs,1)
    i
    for k = 1:size(stim_locs,3)
%         k
%         if i == 1 && k == 1
%             find(any(bsxfun(@minus,target_locs,stim_locs(i,:,k)),2) == 0)
%             size(find(any(bsxfun(@minus,target_locs,stim_locs(i,:,k)),2) == 0))
%         end
        if ~any(isnan(stim_locs(i,:,k)))
            target_ind = find(all(bsxfun(@minus,target_locs,stim_locs(i,:,k)) == 0,2),2);
%             size(target_ind)
            if size(target_ind,1) > 1 || size(target_ind,2) > 1
                disp('fail')
                return
            end
            stim_inds(i,k) = target_ind(1);
        else
            stim_inds(i,k) = NaN;
        end
    end
end