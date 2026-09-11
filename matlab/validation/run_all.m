% run_all.m - Run all MSL MATLAB validation scripts and print a summary.
%
% Usage:
%   matlab -batch "run('matlab/validation/run_all.m')"
%
% Each module script is expected to read comparison data written by the
% corresponding C++ test under test_result/<module>/matlab_compare.
%
% Validation scripts execute in the caller workspace and may overwrite
% variables. The runner therefore stores its state in prefixed variables and
% restores them after every script.

msl_validation_dir = fileparts(mfilename('fullpath'));
msl_validation_modules = {'matrix', 'signal', 'integral', 'difference', ...
    'interp', 'polynomial', 'equation', 'ode'};

msl_validation_results = strings(size(msl_validation_modules));
for msl_validation_index = 1:numel(msl_validation_modules)
    msl_validation_module = msl_validation_modules{msl_validation_index};
    fprintf('\n===== %s =====\n', upper(msl_validation_module));

    script_dir = fullfile(msl_validation_dir, msl_validation_module);
    script_path = fullfile(script_dir, ...
        ['validate_' msl_validation_module '.m']);

    saved_results = msl_validation_results;
    saved_index = msl_validation_index;
    try
        run(script_path);
        msl_validation_results = saved_results;
        msl_validation_index = saved_index;
        msl_validation_results(msl_validation_index) = "PASS";
    catch error_info
        msl_validation_results = saved_results;
        msl_validation_index = saved_index;
        fprintf(2, 'FAILED: %s\n', error_info.message);
        msl_validation_results(msl_validation_index) = "FAIL";
    end
end

fprintf('\nMSL MATLAB validation summary\n');
for msl_validation_index = 1:numel(msl_validation_modules)
    fprintf('%-12s %s\n', msl_validation_modules{msl_validation_index}, ...
        msl_validation_results(msl_validation_index));
end
if license('test', 'Signal_Toolbox') ~= 1
    fprintf(['Note: Signal Processing Toolbox unavailable; ' ...
        'Butterworth/Welch/Savitzky-Golay checks were SKIPPED.\n']);
end

msl_validation_passed = sum(msl_validation_results == "PASS");
msl_validation_failed = sum(msl_validation_results == "FAIL");
fprintf('Total: %d passed, %d failed\n', msl_validation_passed, ...
    msl_validation_failed);

if msl_validation_failed > 0
    error('%d module validation(s) failed.', msl_validation_failed);
end
