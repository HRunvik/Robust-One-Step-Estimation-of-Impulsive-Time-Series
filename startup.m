link = mfilename('fullpath');
link = link(1:end-7); % s-t-a-r-t-u-p is 7 chars

% path to folders
addpath(genpath([link 'Auxil/']));
addpath(genpath([link 'RobustRisk/']));