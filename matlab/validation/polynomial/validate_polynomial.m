if ~exist('script_dir', 'var') || isempty(script_dir)
    script_dir = fileparts(mfilename('fullpath'));
end
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'polynomial', 'matlab_compare');

coeff_file = fullfile(data_dir, 'polynomial_coefficients.txt');
values_file = fullfile(data_dir, 'polynomial_values.txt');
if ~isfile(coeff_file) || ~isfile(values_file)
    error('Missing comparison data. Run "xmake run test_polynomial" first.');
end

coeffs_msl = readmatrix(coeff_file);
values = readmatrix(values_file);

x_fit = [-2.0; -1.0; 0.0; 1.0; 2.0];
y_fit = 1.0 - 2.0 * x_fit + 0.5 * x_fit.^2;

coeffs_matlab_desc = polyfit(x_fit, y_fit, 2);
coeffs_matlab_asc = flip(coeffs_matlab_desc(:));

assert(max(abs(coeffs_msl - coeffs_matlab_asc)) < 1e-10, ...
    'Polynomial fit coefficient comparison failed.');

query = values(:, 1);
expected = 1.0 + 2.0 * query + 3.0 * query.^2;
assert(max(abs(values(:, 2) - expected)) < 1e-12, ...
    'Polynomial object evaluation comparison failed.');
assert(max(abs(values(:, 3) - expected)) < 1e-12, ...
    'polyval comparison failed.');

fprintf('Polynomial validation passed. coeff error = %.3e\n', ...
    max(abs(coeffs_msl - coeffs_matlab_asc)));
