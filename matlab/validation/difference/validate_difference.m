script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'difference', 'matlab_compare');

required = {'difference_vector.txt', 'difference_matrix.txt', ...
    'difference_row_diff.txt', 'difference_col_diff.txt', ...
    'difference_nonuniform_input.txt', ...
    'difference_forward_nonuniform.txt', ...
    'difference_central_nonuniform.txt', ...
    'difference_second_derivative.txt', ...
    'difference_field.txt', 'difference_grad_x.txt', ...
    'difference_grad_y.txt', 'difference_laplacian.txt', ...
    'difference_vector_x.txt', 'difference_vector_y.txt', ...
    'difference_divergence.txt', 'difference_curl.txt', ...
    'difference_savgol.txt'};
for i = 1:numel(required)
    if ~isfile(fullfile(data_dir, required{i}))
        error('Missing comparison data. Run "xmake run test_difference" first.');
    end
end

read = @(name) readmatrix(fullfile(data_dir, name));

% ---- Vector diff and central gradient ----
vec_data = read('difference_vector.txt');
y = vec_data(:, 1);
diff_msl = vec_data(1:end-1, 2);
central_msl = vec_data(:, 3);

diff_ref = diff(y);
central_ref = [
    y(2) - y(1);
    (y(3:end) - y(1:end-2)) / 2.0;
    y(end) - y(end-1)
];

assert(max(abs(diff_msl - diff_ref)) < 1e-12, ...
    'Vector diff comparison failed.');
assert(max(abs(central_msl - central_ref)) < 1e-12, ...
    'Central gradient comparison failed.');

% ---- Matrix diff ----
M = read('difference_matrix.txt');
row_diff = read('difference_row_diff.txt');
col_diff = read('difference_col_diff.txt');

assert(max(abs(row_diff - diff(M, 1, 1)), [], 'all') < 1e-12, ...
    'Matrix row diff comparison failed.');
assert(max(abs(col_diff - diff(M, 1, 2)), [], 'all') < 1e-12, ...
    'Matrix col diff comparison failed.');

% ---- Non-uniform gradients on f(x) = x^2 ----
nonuniform = read('difference_nonuniform_input.txt');
nu_x = nonuniform(:, 1);
nu_y = nonuniform(:, 2);
nu_forward = read('difference_forward_nonuniform.txt');
nu_central = read('difference_central_nonuniform.txt');

forward_ref = diff(nu_y) ./ diff(nu_x);
assert(max(abs(nu_forward - forward_ref)) < 1e-12, ...
    'Non-uniform forward gradient comparison failed.');

% Interior weighted central differences are exact for x^2; edges are
% first-order one-sided differences
central_expected = 2 * nu_x;
central_expected(1) = (nu_y(2) - nu_y(1)) / (nu_x(2) - nu_x(1));
central_expected(end) = (nu_y(end) - nu_y(end-1)) / (nu_x(end) - nu_x(end-1));
assert(max(abs(nu_central - central_expected)) < 1e-12, ...
    'Non-uniform central gradient comparison failed.');

% ---- Second derivative of a quadratic ----
second = read('difference_second_derivative.txt');
second_error = max(abs(second(:, 3) - 2.0));
assert(second_error < 1e-12, ...
    sprintf('Second derivative comparison failed (error %.3e).', ...
    second_error));

% ---- 2D gradient and Laplacian on f(x, y) = x^2 + y^2 ----
field = read('difference_field.txt');
grad_x = read('difference_grad_x.txt');
grad_y = read('difference_grad_y.txt');
laplacian_field = read('difference_laplacian.txt');

[n_rows, n_cols] = size(field);
x_coord = (0:n_cols-1);
y_coord = (0:n_rows-1)';
expected_grad_x = repmat(2 * x_coord, n_rows, 1);
expected_grad_y = repmat(2 * y_coord, 1, n_cols);

interior = 2:(n_rows-1);
interior_cols = 2:(n_cols-1);
grad_x_error = max(abs(grad_x(interior, interior_cols) - ...
    expected_grad_x(interior, interior_cols)), [], 'all');
grad_y_error = max(abs(grad_y(interior, interior_cols) - ...
    expected_grad_y(interior, interior_cols)), [], 'all');
assert(grad_x_error < 1e-12 && grad_y_error < 1e-12, ...
    '2D gradient comparison failed.');

laplacian_error = max(abs(laplacian_field(interior, interior_cols) - 4.0), ...
    [], 'all');
assert(laplacian_error < 1e-12, ...
    sprintf('Laplacian comparison failed (error %.3e).', laplacian_error));

% ---- Divergence and curl on linear fields ----
divergence_field = read('difference_divergence.txt');
curl_field = read('difference_curl.txt');
divergence_error = max(abs(divergence_field - 2.0), [], 'all');
curl_error = max(abs(curl_field - 2.0), [], 'all');
assert(divergence_error < 1e-12, ...
    sprintf('Divergence comparison failed (error %.3e).', divergence_error));
assert(curl_error < 1e-12, ...
    sprintf('Curl comparison failed (error %.3e).', curl_error));

% ---- Savitzky-Golay derivative (requires Signal Processing Toolbox) ----
has_sgolay = exist('sgolay', 'file') == 2 && license('test', 'Signal_Toolbox') == 1;
if has_sgolay
    savgol = read('difference_savgol.txt');
    savgol_ref = sgolay_derivative(savgol(:, 2), 2, 5, 1.0);
    savgol_error = max(abs(savgol(:, 3) - savgol_ref));
    assert(savgol_error < 1e-10, ...
        sprintf('Savitzky-Golay derivative comparison failed (error %.3e).', ...
        savgol_error));
    fprintf(['Difference validation passed. vector diff error = %.3e, ' ...
        'SG derivative error = %.3e\n'], ...
        max(abs(diff_msl - diff_ref)), savgol_error);
else
    fprintf(['Difference validation passed. vector diff error = %.3e\n' ...
        'Difference validation: Savitzky-Golay check SKIPPED ' ...
        '(Signal Processing Toolbox not available).\n'], ...
        max(abs(diff_msl - diff_ref)));
end

function deriv = sgolay_derivative(y, order, framelen, dx)
    [~, g] = sgolay(order, framelen);
    n = numel(y);
    half = (framelen - 1) / 2;
    deriv = zeros(n, 1);
    for i = 1:n
        if i <= half
            window = y(1:framelen);
            x0 = (i - 1) - half;
        elseif i > n - half
            window = y(n-framelen+1:n);
            x0 = (i - 1) - (n - framelen) - half;
        else
            window = y(i-half:i+half);
            x0 = 0;
        end
        value = 0;
        for k = 1:order
            value = value + k * x0^(k-1) * (g(:, k+1)' * window);
        end
        deriv(i) = value / dx;
    end
end
