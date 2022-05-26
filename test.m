[filepath, ~, ~] = fileparts(mfilename('fullpath')); % get the directory of the current m-file
cd(filepath);
X = randn(101, 101, 101);
S = rpfft3(X)

