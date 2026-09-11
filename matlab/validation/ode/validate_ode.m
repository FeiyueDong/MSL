if ~exist('script_dir', 'var') || isempty(script_dir)
    script_dir = fileparts(mfilename('fullpath'));
end
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'ode', 'matlab_compare');

exp_file = fullfile(data_dir, 'ode_exp_eval.txt');
osc_file = fullfile(data_dir, 'ode_oscillator_rk4.txt');

if ~isfile(exp_file) || ~isfile(osc_file)
    error('Missing comparison data. Run "xmake run test_ode" first.');
end

exp_data = readmatrix(exp_file);
t = exp_data(:, 1);
euler_y = exp_data(:, 2);
heun_y = exp_data(:, 3);
rk4_y = exp_data(:, 4);
ode45_y = exp_data(:, 5);

exact_exp = exp(t);
assert(max(abs(euler_y - exact_exp)) < 2e-1, ...
    'Euler exponential comparison failed.');
assert(max(abs(heun_y - exact_exp)) < 5e-3, ...
    'Heun exponential comparison failed.');
assert(max(abs(rk4_y - exact_exp)) < 3e-6, ...
    'RK4 exponential comparison failed.');

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[~, y_matlab] = ode45(@(tt, yy) yy, t, 1.0, opts);
assert(max(abs(ode45_y - y_matlab)) < 1e-7, ...
    'ODE45 MATLAB comparison failed.');

osc_data = readmatrix(osc_file);
t_osc = osc_data(:, 1);
x_msl = osc_data(:, 2);
v_msl = osc_data(:, 3);

x_exact = cos(t_osc);
v_exact = -sin(t_osc);

assert(max(abs(x_msl - x_exact)) < 1e-8, ...
    'RK4 oscillator position comparison failed.');
assert(max(abs(v_msl - v_exact)) < 1e-8, ...
    'RK4 oscillator velocity comparison failed.');

fprintf('ODE validation passed. max ode45 error = %.3e\n', ...
    max(abs(ode45_y - y_matlab)));
