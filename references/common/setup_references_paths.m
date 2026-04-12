function refs_root = setup_references_paths()
%SETUP_REFERENCES_PATHS Add the repo reference folders needed by scripts.

this_dir = fileparts(mfilename('fullpath'));
refs_root = find_references_root(this_dir);

addpath(fullfile(refs_root, 'common'));
addpath(fullfile(refs_root, 'components'));
addpath(fullfile(refs_root, 'simulink', 'builders'));
addpath(fullfile(refs_root, 'simulink', 'runs'));
addpath(fullfile(refs_root, 'simulink', 'tests'));
addpath(fullfile(refs_root, 'analysis'));
addpath(fullfile(refs_root, 'analysis', 'jitter'));
addpath(fullfile(refs_root, 'third_party', 'dstoolbox'));
end

function refs_root = find_references_root(start_dir)
current_dir = start_dir;
while true
    if exist(fullfile(current_dir, 'third_party', 'dstoolbox'), 'dir')
        refs_root = current_dir;
        return;
    end

    parent_dir = fileparts(current_dir);
    if strcmp(parent_dir, current_dir)
        error('setup_references_paths:RootNotFound', ...
            'Could not locate the references root from %s.', start_dir);
    end
    current_dir = parent_dir;
end
end
