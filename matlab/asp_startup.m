function asp_startup()
%ASP_STARTUP Put every ASP directory on the MATLAB/Octave path.
%
%   Run this once per session from anywhere:
%       run('<repo>/matlab/asp_startup.m')
%
%   The reference model deliberately uses no toolboxes.  It runs unmodified
%   on MATLAB R2018b+ and GNU Octave 7+.

here = fileparts(mfilename('fullpath'));

subdirs = {'config','fx','model','core','analysis','verify','studies','export','golden'};
for k = 1:numel(subdirs)
    p = fullfile(here, subdirs{k});
    if exist(p, 'dir')
        addpath(p);
    end
end
addpath(here);

end
