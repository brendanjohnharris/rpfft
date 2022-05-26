[filepath, ~, ~] = fileparts(mfilename('fullpath')); % get the directory of the current m-file
cd(filepath);
X = 1
S = rpfft3(X, "subsample", false)

