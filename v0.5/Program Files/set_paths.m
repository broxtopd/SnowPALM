function set_paths(base_path,opt)
% Set up model paths
if strcmpi(opt,'add')
    addpath([base_path filesep 'Model'])
    addpath([base_path filesep 'Misc'])
elseif strcmpi(opt,'remove')
    rmpath([base_path filesep 'Model'])
    rmpath([base_path filesep 'Misc'])
end
