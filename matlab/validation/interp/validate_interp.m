if ~exist('script_dir', 'var') || isempty(script_dir)
    script_dir = fileparts(mfilename('fullpath'));
end
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'interp', 'matlab_compare');

data_file = fullfile(data_dir, 'interp_data.txt');
interp_file = fullfile(data_dir, 'interp_results.txt');
nearest_file = fullfile(data_dir, 'interp_nearest.txt');
if ~isfile(data_file) || ~isfile(interp_file) || ~isfile(nearest_file)
    error('Missing comparison data. Run "xmake run test_interp" first.');
end

data_points = readmatrix(data_file);
x = data_points(:, 1);
y = data_points(:, 2);

data = readmatrix(interp_file);
xq = data(:, 1);

linear_ref = interp1(x, y, xq, 'linear');
spline_ref = interp1(x, y, xq, 'spline');
pchip_ref = interp1(x, y, xq, 'pchip');
makima_ref = interp1(x, y, xq, 'makima');
poly_coeffs = polyfit(x, y, numel(x) - 1);
poly_ref = polyval(poly_coeffs, xq);

assert(max(abs(data(:, 2) - linear_ref)) < 1e-12, ...
    'Linear interpolation comparison failed.');
assert(max(abs(data(:, 3) - spline_ref)) < 1e-10, ...
    'Spline interpolation comparison failed.');
assert(max(abs(data(:, 4) - pchip_ref)) < 1e-12, ...
    'PCHIP interpolation comparison failed.');
assert(max(abs(data(:, 5) - makima_ref)) < 1e-12, ...
    'Akima/MAKIMA interpolation comparison failed.');
assert(max(abs(data(:, 6) - poly_ref)) < 1e-8, ...
    'Polynomial interpolation comparison failed.');

nearest = readmatrix(nearest_file);
nearest_ref = interp1(x, y, nearest(:, 1), 'nearest');
assert(max(abs(nearest(:, 2) - nearest_ref)) < 1e-12, ...
    'Nearest interpolation comparison failed.');

fprintf(['Interp validation passed. linear error = %.3e, spline error = ' ...
    '%.3e, pchip error = %.3e, makima error = %.3e, polynomial error = ' ...
    '%.3e\n'], max(abs(data(:, 2) - linear_ref)), ...
    max(abs(data(:, 3) - spline_ref)), ...
    max(abs(data(:, 4) - pchip_ref)), ...
    max(abs(data(:, 5) - makima_ref)), ...
    max(abs(data(:, 6) - poly_ref)));
