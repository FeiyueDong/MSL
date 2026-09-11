script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'signal', 'matlab_compare');

fft_file = fullfile(data_dir, 'signal_fft.txt');
window_file = fullfile(data_dir, 'signal_windows.txt');
filter_file = fullfile(data_dir, 'signal_filter.txt');
if ~isfile(fft_file) || ~isfile(window_file) || ~isfile(filter_file)
    error('Missing comparison data. Run "xmake run test_signal" first.');
end

fft_data = readmatrix(fft_file);
x = fft_data(:, 1);
X_msl = fft_data(:, 2) + 1i * fft_data(:, 3);
restored = fft_data(:, 4);
X_ref = fft(x);

assert(max(abs(X_msl - X_ref)) < 1e-12, ...
    'FFT comparison failed.');
assert(max(abs(restored - real(ifft(X_ref)))) < 1e-12, ...
    'IFFT comparison failed.');

windows = readmatrix(window_file);
n = size(windows, 1);
k = (0:n-1)';
hann_ref = 0.5 - 0.5 * cos(2 * pi * k / (n - 1));
hamming_ref = 0.54 - 0.46 * cos(2 * pi * k / (n - 1));

assert(max(abs(windows(:, 1) - hann_ref)) < 1e-12, ...
    'Hann window comparison failed.');
assert(max(abs(windows(:, 2) - hamming_ref)) < 1e-12, ...
    'Hamming window comparison failed.');

filter_data = readmatrix(filter_file);
filter_ref = filter([0.5, 0.5], 1.0, filter_data(:, 1));
assert(max(abs(filter_data(:, 2) - filter_ref)) < 1e-12, ...
    'Filter comparison failed.');

% ---- Butterworth / filter / filtfilt / freqz cross validation ----
butter_cases_file = fullfile(data_dir, 'butterworth_cases.txt');
butter_response_file = fullfile(data_dir, 'butterworth_response.txt');
butter_filter_file = fullfile(data_dir, 'butterworth_filter.txt');
butter_filtfilt_file = fullfile(data_dir, 'butterworth_filtfilt.txt');

has_signal_toolbox = license('test', 'Signal_Toolbox') == 1 ...
    && exist('butter', 'file') == 2 ...
    && exist('freqz', 'file') == 2 ...
    && exist('filtfilt', 'file') == 2;

if ~has_signal_toolbox
    fprintf('Signal validation passed. FFT error = %.3e\n', ...
        max(abs(X_msl - X_ref)));
    fprintf(['Signal validation: Butterworth/filter/filtfilt/freqz ', ...
        'checks SKIPPED (Signal Processing Toolbox not available).\n']);
elseif ~all(cellfun(@isfile, {butter_cases_file, butter_response_file, ...
        butter_filter_file, butter_filtfilt_file}))
    fprintf('Signal validation passed. FFT error = %.3e\n', ...
        max(abs(X_msl - X_ref)));
    fprintf(['Signal validation: Butterworth comparison data missing; ', ...
        'run "xmake run test_signal".\n']);
else
    cases = readmatrix(butter_cases_file);
    response = readmatrix(butter_response_file);
    filter_cmp = readmatrix(butter_filter_file);
    filtfilt_cmp = readmatrix(butter_filtfilt_file);

    f = response(:, 1);
    x_butter = filter_cmp(:, 1);
    n_cases = size(cases, 1);
    max_response_error = 0;
    max_filter_error = 0;
    max_filtfilt_error = 0;

    for case_index = 1:n_cases
        [b, a] = design_butter_case(cases(case_index, :));

        h_ref = freqz(b, a, pi * f');
        response_error = max(abs(abs(h_ref(:)) - response(:, 1 + case_index)));
        max_response_error = max(max_response_error, response_error);
        assert(response_error < 1e-8, ...
            ['Butterworth response comparison failed for case %d ' ...
             '(error %.3e).'], case_index, response_error);

        y_ref = filter(b, a, x_butter);
        filter_error = max(abs(y_ref - filter_cmp(:, 1 + case_index)));
        max_filter_error = max(max_filter_error, filter_error);
        assert(filter_error < 1e-8, ...
            ['Butterworth filter comparison failed for case %d ' ...
             '(error %.3e).'], case_index, filter_error);

        y_ff = filtfilt(b, a, x_butter);
        filtfilt_error = max(abs(y_ff - filtfilt_cmp(:, 1 + case_index)));
        max_filtfilt_error = max(max_filtfilt_error, filtfilt_error);
        assert(filtfilt_error < 1e-6, ...
            ['Butterworth filtfilt comparison failed for case %d ' ...
             '(error %.3e).'], case_index, filtfilt_error);
    end

    fprintf(['Signal validation passed. FFT error = %.3e, Butterworth ' ...
        'response error = %.3e, filter error = %.3e, filtfilt error = ' ...
        '%.3e\n'], max(abs(X_msl - X_ref)), max_response_error, ...
        max_filter_error, max_filtfilt_error);
end

function [b, a] = design_butter_case(case_row)
    type = case_row(1);
    order = case_row(2);
    low = case_row(3);
    high = case_row(4);
    switch type
        case 1
            [b, a] = butter(order, low);
        case 2
            [b, a] = butter(order, low, 'high');
        case 3
            [b, a] = butter(order, [low, high]);
        case 4
            [b, a] = butter(order, [low, high], 'stop');
        otherwise
            error('Unknown Butterworth case type %d.', type);
    end
end
