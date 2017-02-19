function stacked_stacks = stack_stacks(stacks)

num_stacks = length(stacks);

length_trial = size(stacks{1},2);
num_trials = sum(cellfun(@(x) size(x,1),stacks));
stacked_stacks = zeros(num_trials,length_trial);

count = 1;
for i = 1:num_stacks
    stacked_stacks(count:count+size(stacks{i},1)-1,:) = stacks{i};
    count = count+size(stacks{i},1);
end

            