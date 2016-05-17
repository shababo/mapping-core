function hack(slice_num, cell_num, run_num, date, min_or_max, tag, varargin)

if ~isstr(slice_num)
    slice_num = num2str(slice_num);
end
    
mapping_eval = [''...
'[traces] = '...
    'get_mapping_data(''' date '_slice' slice_num '_cell' num2str(cell_num) '' tag '.mat'',' num2str(run_num) ',''' min_or_max ''',1,1'];

for i = 1:2:length(varargin)
    mapping_eval = [mapping_eval ',''' varargin{i} ''',' num2str(varargin{i+1})];
end
mapping_eval = [mapping_eval ');']

eval(mapping_eval);

assignin('base',...
    ['traces_by_location_' date '_s' slice_num 'c' num2str(cell_num) '_r' num2str(run_num)],...
    traces);
% assignin('base',...
%     ['current_image_' date '_s' slice_num 'c' num2str(cell_num) '_r' num2str(run_num)],...
%     current_image);