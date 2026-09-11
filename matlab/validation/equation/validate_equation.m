if ~exist('script_dir', 'var') || isempty(script_dir)
    script_dir = fileparts(mfilename('fullpath'));
end
repo_root = fullfile(script_dir, '..', '..', '..');
data_file = fullfile(repo_root, 'test_result', 'equation', ...
    'matlab_compare', 'equation_roots.txt');

if ~isfile(data_file)
    error('Missing comparison data. Run "xmake run test_equation" first.');
end

data = readmatrix(data_file);
roots_msl = data(:, 1);
residuals = abs(data(:, 2));

sqrt2_ref = fzero(@(x) x.^2 - 2.0, [0.0, 2.0]);
fixed_point_ref = fzero(@(x) cos(x) - x, [0.0, 1.0]);

expected = [
    sqrt2_ref;
    sqrt2_ref;
    sqrt2_ref;
    sqrt2_ref;
    sqrt2_ref;
    fixed_point_ref
];

root_error = abs(roots_msl - expected);

assert(all(root_error(1:5) < 1e-10), ...
    'sqrt(2) root comparison failed.');
assert(root_error(6) < 1e-10, ...
    'cos(x)-x root comparison failed.');
assert(all(residuals < 1e-9), ...
    'Root residual comparison failed.');

fprintf('Equation validation passed. max root error = %.3e\n', ...
    max(root_error));
