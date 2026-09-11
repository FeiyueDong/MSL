if ~exist('script_dir', 'var') || isempty(script_dir)
    script_dir = fileparts(mfilename('fullpath'));
end
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'integral', 'matlab_compare');

summary_file = fullfile(data_dir, 'integral_summary.txt');
cum_file = fullfile(data_dir, 'integral_cumulative.txt');
if ~isfile(summary_file) || ~isfile(cum_file)
    error('Missing comparison data. Run "xmake run test_integral" first.');
end

summary = readmatrix(summary_file);
cum_data = readmatrix(cum_file);

x = cum_data(:, 1);
y = cum_data(:, 2);
cum_uniform = cum_data(:, 3);
cum_nonuniform = cum_data(:, 4);

expected = [
    trapz(y);
    trapz(x, y);
    integral(@(xx) xx.^2, 0.0, 2.0);
    integral(@(xx) xx.^2, 0.0, 2.0);
    integral(@(xx) xx.^2, 0.0, 4.0);
    integral(@sin, 0.0, pi)
];

assert(max(abs(summary - expected)) < 1e-10, ...
    'Integral summary comparison failed.');
assert(max(abs(cum_uniform - cumtrapz(y))) < 1e-12, ...
    'Uniform cumtrapz comparison failed.');
assert(max(abs(cum_nonuniform - cumtrapz(x, y))) < 1e-12, ...
    'Non-uniform cumtrapz comparison failed.');

% ---- Cumulative Simpson against analytic primitives ----
cs_uniform = readmatrix(fullfile(data_dir, 'integral_cumsimpson_uniform.txt'));
cs_x = cs_uniform(:, 1);
cs_uniform_error = max(abs(cs_uniform(:, 3) - cs_x.^3 / 3));
assert(cs_uniform_error < 1e-12, ...
    sprintf('Uniform cumulative Simpson comparison failed (error %.3e).', ...
    cs_uniform_error));

cs_nonuniform_data = readmatrix(fullfile(data_dir, ...
    'integral_cumsimpson_nonuniform.txt'));
cs_nonuniform_error = max(abs(cs_nonuniform_data(:, 3) - ...
    cs_nonuniform_data(:, 1).^3 / 3));
assert(cs_nonuniform_error < 1e-12, ...
    sprintf(['Non-uniform cumulative Simpson comparison failed ' ...
    '(error %.3e).'], cs_nonuniform_error));

% For an even number of points the last interval uses the documented
% trapezoidal fallback
cs_even = readmatrix(fullfile(data_dir, 'integral_cumsimpson_even.txt'));
cs_even_expected = cs_even(:, 1).^3 / 3;
cs_even_expected(end) = cs_even(end - 1, 1)^3 / 3 + ...
    0.5 * (cs_even(end - 1, 2) + cs_even(end, 2));
cs_even_error = max(abs(cs_even(:, 3) - cs_even_expected));
assert(cs_even_error < 1e-12, ...
    sprintf(['Even-length cumulative Simpson comparison failed ' ...
    '(error %.3e).'], cs_even_error));

% ---- Matrix integration: Simpson vs analytic, trapz vs MATLAB ----
integrand = readmatrix(fullfile(data_dir, 'integral_matrix.txt'));
matrix_results = readmatrix(fullfile(data_dir, ...
    'integral_matrix_results.txt'));
simpson_expected = [integral(@(t) t.^2, 0, 4), integral(@(t) t.^3, 0, 4)];
assert(max(abs(matrix_results(1, :) - simpson_expected)) < 1e-10, ...
    'Matrix Simpson comparison against analytic integrals failed.');
assert(max(abs(matrix_results(2, :) - trapz(integrand))) < 1e-12, ...
    'Matrix trapz comparison against MATLAB trapz failed.');

fprintf(['Integral validation passed. summary error = %.3e, cumtrapz ' ...
    'error = %.3e, cumulative Simpson error = %.3e\n'], ...
    max(abs(summary - expected)), ...
    max(abs(cum_uniform - cumtrapz(y))), ...
    max([cs_uniform_error, cs_nonuniform_error, cs_even_error]));
