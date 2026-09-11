if ~exist('script_dir', 'var') || isempty(script_dir)
    script_dir = fileparts(mfilename('fullpath'));
end
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'matrix', 'matlab_compare');

required = {'matrix_A.txt', 'matrix_B.txt', 'matrix_product.txt', ...
    'matrix_svd_input.txt', 'matrix_svd_reconstruction.txt', ...
    'matrix_svd_U.txt', 'matrix_svd_S.txt', 'matrix_svd_V.txt', ...
    'matrix_lu_A.txt', 'matrix_lu_L.txt', 'matrix_lu_U.txt', ...
    'matrix_lu_permutation.txt', ...
    'matrix_qr_tall_A.txt', 'matrix_qr_tall_Q.txt', 'matrix_qr_tall_R.txt', ...
    'matrix_qr_wide_A.txt', 'matrix_qr_wide_Q.txt', 'matrix_qr_wide_R.txt', ...
    'matrix_eig_A.txt', 'matrix_eig_V.txt', 'matrix_eig_D.txt', ...
    'matrix_geig_A.txt', 'matrix_geig_B.txt', 'matrix_geig_V.txt', ...
    'matrix_geig_D.txt', ...
    'matrix_pinv_A.txt', 'matrix_pinv.txt', ...
    'matrix_inverse_A.txt', 'matrix_inverse.txt', ...
    'matrix_determinant.txt', ...
    'matrix_tsvd_A.txt', 'matrix_tsvd_U.txt', 'matrix_tsvd_V.txt', ...
    'matrix_tsvd_singular_values.txt'};
for i = 1:numel(required)
    if ~isfile(fullfile(data_dir, required{i}))
        error('Missing comparison data. Run "xmake run test_matrix" first.');
    end
end

read = @(name) readmatrix(fullfile(data_dir, name));
readc = @(name) read_complex(fullfile(data_dir, name));

A = read('matrix_A.txt');
B = read('matrix_B.txt');
C = read('matrix_product.txt');
assert(max(abs(C - A * B), [], 'all') < 1e-12, ...
    'Matrix product comparison failed.');

% ---- SVD ----
SVD_input = read('matrix_svd_input.txt');
SVD_reconstructed = read('matrix_svd_reconstruction.txt');
SVD_U = read('matrix_svd_U.txt');
SVD_S = read('matrix_svd_S.txt');
SVD_V = read('matrix_svd_V.txt');
assert(max(abs(SVD_reconstructed - SVD_input), [], 'all') < 1e-10, ...
    'Matrix SVD reconstruction comparison failed.');
[~, S_ref, ~] = svd(SVD_input);
assert(max(abs(diag(SVD_S) - diag(S_ref))) < 1e-10, ...
    'SVD singular value comparison failed.');
assert(max(abs(SVD_U' * SVD_U - eye(size(SVD_U, 2))), [], 'all') < 1e-10, ...
    'SVD left singular vectors are not orthonormal.');
assert(max(abs(SVD_V' * SVD_V - eye(size(SVD_V, 2))), [], 'all') < 1e-10, ...
    'SVD right singular vectors are not orthonormal.');

% ---- LU: P * A = L * U ----
LU_A = read('matrix_lu_A.txt');
LU_L = read('matrix_lu_L.txt');
LU_U = read('matrix_lu_U.txt');
LU_p = read('matrix_lu_permutation.txt');
assert(max(abs(LU_A(LU_p, :) - LU_L * LU_U), [], 'all') < 1e-12, ...
    'LU comparison P * A = L * U failed.');
assert(max(abs(LU_L - tril(LU_L)), [], 'all') < 1e-12 ...
    && max(abs(diag(LU_L) - 1)) < 1e-12, ...
    'LU L is not unit lower triangular.');
assert(max(abs(LU_U - triu(LU_U)), [], 'all') < 1e-12, ...
    'LU U is not upper triangular.');

% ---- QR (economy) ----
qr_cases = {'tall', 'wide'};
for k = 1:numel(qr_cases)
    name = qr_cases{k};
    QR_A = read(sprintf('matrix_qr_%s_A.txt', name));
    QR_Q = read(sprintf('matrix_qr_%s_Q.txt', name));
    QR_R = read(sprintf('matrix_qr_%s_R.txt', name));
    QR_rank = min(size(QR_A));
    assert(max(abs(QR_A - QR_Q * QR_R), [], 'all') < 1e-10, ...
        sprintf('QR reconstruction failed for %s matrices.', name));
    assert(max(abs(QR_Q' * QR_Q - eye(QR_rank)), [], 'all') < 1e-10, ...
        sprintf('QR Q is not orthonormal for %s matrices.', name));
    assert(max(abs(QR_R - triu(QR_R)), [], 'all') < 1e-10, ...
        sprintf('QR R is not upper triangular for %s matrices.', name));
end

% ---- Eigenvalues: A * V = V * D ----
EIG_A = read('matrix_eig_A.txt');
EIG_V = readc('matrix_eig_V.txt');
EIG_D = readc('matrix_eig_D.txt');
assert(max(abs(EIG_A * EIG_V - EIG_V * EIG_D), [], 'all') < 1e-10, ...
    'Eigenvalue residual A * V - V * D is too large.');
eig_ref = sort(eig(EIG_A), 'ComparisonMethod', 'real');
eig_msl = sort(diag(EIG_D), 'ComparisonMethod', 'real');
assert(max(abs(eig_ref - eig_msl)) < 1e-10, ...
    'Eigenvalue comparison against MATLAB eig failed.');

% ---- Generalized eigenvalues: A * V = B * V * D ----
GEIG_A = read('matrix_geig_A.txt');
GEIG_B = read('matrix_geig_B.txt');
GEIG_V = readc('matrix_geig_V.txt');
GEIG_D = readc('matrix_geig_D.txt');
assert(max(abs(GEIG_A * GEIG_V - GEIG_B * GEIG_V * GEIG_D), [], 'all') < 1e-10, ...
    'Generalized eigenvalue residual is too large.');
geig_ref = sort(eig(GEIG_A, GEIG_B), 'ComparisonMethod', 'real');
geig_msl = sort(diag(GEIG_D), 'ComparisonMethod', 'real');
assert(max(abs(geig_ref - geig_msl)) < 1e-10, ...
    'Generalized eigenvalue comparison against MATLAB eig failed.');

% ---- Pseudoinverse and Moore-Penrose identities ----
PINV_A = read('matrix_pinv_A.txt');
PINV = read('matrix_pinv.txt');
assert(max(abs(PINV - pinv(PINV_A)), [], 'all') < 1e-10, ...
    'Pseudoinverse comparison against MATLAB pinv failed.');
assert(max(abs(PINV_A * PINV * PINV_A - PINV_A), [], 'all') < 1e-10 ...
    && max(abs(PINV * PINV_A * PINV - PINV), [], 'all') < 1e-10, ...
    'Moore-Penrose identities failed.');

% ---- Inverse and determinant ----
INV_A = read('matrix_inverse_A.txt');
INV = read('matrix_inverse.txt');
DET = read('matrix_determinant.txt');
assert(max(abs(INV - inv(INV_A)), [], 'all') < 1e-12, ...
    'Inverse comparison against MATLAB inv failed.');
assert(abs(DET - det(INV_A)) < 1e-12, ...
    'Determinant comparison against MATLAB det failed.');

% ---- Truncated SVD ----
TSVD_A = read('matrix_tsvd_A.txt');
TSVD_U = read('matrix_tsvd_U.txt');
TSVD_V = read('matrix_tsvd_V.txt');
TSVD_s = read('matrix_tsvd_singular_values.txt');
tsvd_ref = svd(TSVD_A);
assert(max(abs(TSVD_s - tsvd_ref(1:numel(TSVD_s)))) < 1e-8, ...
    'Truncated SVD singular value comparison failed.');
assert(max(abs(TSVD_U' * TSVD_U - eye(size(TSVD_U, 2))), [], 'all') < 1e-10, ...
    'Truncated SVD U is not orthonormal.');
assert(max(abs(TSVD_V' * TSVD_V - eye(size(TSVD_V, 2))), [], 'all') < 1e-10, ...
    'Truncated SVD V is not orthonormal.');

fprintf(['Matrix validation passed. product error = %.3e, LU residual = ' ...
    '%.3e, SVD singular value error = %.3e\n'], ...
    max(abs(C - A * B), [], 'all'), ...
    max(abs(LU_A(LU_p, :) - LU_L * LU_U), [], 'all'), ...
    max(abs(diag(SVD_S) - diag(S_ref))));

function M = read_complex(path)
    raw = readmatrix(path);
    M = raw(:, 1:2:end) + 1i * raw(:, 2:2:end);
end
