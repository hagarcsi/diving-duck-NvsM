% PowerPhaseAccUnN  –  acceleration during the POWER PHASE: video vs logger
%
% For every digitized trial:
%   VIDEO : acceleration of Pt9 (the logger's position on the duck) along
%           the body axis (tail -> neck), from the 3D digitization.
%   LOGGER: forward acceleration from analyze_peaks_and_calibration
%           (ForwardRaw, low-passed here at cutHz like the video - NOT the
%           2 Hz version, which would smear a ~120 ms window), at the same
%           moments, using tauX from
%           MagSyncUnN:  tLogger = tVideo + tauX.
%
% Compared per power phase (not point by point, because the sync is only
% known to a few tens of ms):
%   aMean = mean acceleration over the power phase       (m/s^2)
%   dv    = change of forward speed over the power phase (m/s) = aMean * T
% plus a sensitivity check: the logger mean recomputed with the sync
% shifted by -shiftMax..+shiftMax frames. If that range is small compared
% with the video-logger difference, the sync uncertainty does not matter.
%
% Inputs (all made by the earlier steps):
%   MagDistUnN_results.mat  – Pt9 position per frame (P, incl. estimated)
%   MagSyncUnN_results.mat  – tauX per trial
%   <logger>_acce_results.mat (ResOnA) – forward acceleration of the logger
%   <trial>xyzpts.csv in resFolder – neck (Pt2) and tail base (Pt3), for
%                                    the body axis
% Power-phase frames: ORIGINAL video frames (= row in the xyzpts file),
% from the summary spreadsheet, columns "power phase start/end (frm)".
% (HemiduckUnN uses the "MATLAB" columns, which are shifted where xlsread
% dropped leading empty rows – NOT the same frame numbers.)
%
% NOTE: HemiduckUnN's ACCE is the acceleration of the BODY CENTRE (mid
% neck-tail), not Pt9, and is not smoothed – do not compare it with the
% logger. This code computes the Pt9 acceleration itself.

function [Cmp] = PowerPhaseAccUnN()

close all

FR       = 250;     % video, fps
FRacc    = 100;     % logger, Hz
resFolder = 'C:\Users\hagar\OneDrive\Desktop\for Phd!\question3\diving';
accFolder = 'C:\Users\hagar\OneDrive\Desktop\for Phd!\נטה אדומת ראש\שחייה מתחת למים\20260706\אחרי כיול';   % xyzpts files

cutHz      = 10;    % low-pass cutoff (Hz), 4th-order Butterworth, zero-phase.
                    % Applied to Pt9 position BEFORE the 2 derivatives (the
                    % digitizing errors become huge in acceleration otherwise)
                    % and, for a fair comparison, to the accelerometer too.
                    % Smaller = smoother, but real fast changes are lost.
axisMode   = 'mean';  % 'mean'     = ONE forward direction per trial: the mean
                     %              tail->neck direction (as BoL in HemiduckUnN2).
                     %              Stable: no wobble that drags lateral (Y)
                     %              noise into X.
                     % 'perFrame' = tail->neck recomputed every frame (old way)
axisSmooth = 15;    % gaussian window (frames) on the body axis, 'perFrame' only
shiftMax   = 15;    % frames (= 60 ms), range of the sync sensitivity check
ctx        = 0.25;  % s of context shown before/after the power phase
doPlot     = true;

% power phase, ORIGINAL frames (start, end)
PP = {
'duck1round01xyzpts.csv'    84  118
'duck1round02xyzpts.csv'    99  135
'duck1round04xyzpts.csv'   142  177
'duck1round07xyzpts.csv'   127  160
'duck1round08xyzpts.csv'   132  164
'duck1round09xyzpts.csv'   147  181
'duck1round12xyzpts.csv'    95  130
'duck3round01xyzpts.csv'   137  167
'duck3round03xyzpts.csv'   157  195
'duck3round06xyzpts.csv'   169  201
'duck3round07xyzpts.csv'    96  129
'duck3round08xyzpts.csv'   121  159
'duck3round11xyzpts.csv'   155  194
'duck4round01xyzpts.csv'   131  159
'duck4round02xyzpts.csv'    93  120
'duck4round03xyzpts.csv'    98  126
'duck4round05xyzpts.csv'    85  119
'duck4round10xyzpts.csv'   120  147
'duck5round01xyzpts.csv'   101  130
'duck5round02xyzpts.csv'    97  127
'duck5round06xyzpts.csv'   111  146
'duck5round17xyzpts.csv'   109  134
'duck5round18xyzpts.csv'    96  129
};

load(fullfile(resFolder, 'MagDistUnN_results.mat'), 'MagUnN');
load(fullfile(resFolder, 'MagSyncUnN_results.mat'), 'Sync');

Cmp = struct([]);
cacheName = {};  cacheR = {};          % logger results, loaded once each
shifts = -shiftMax : shiftMax;
PL = cell(0, 15);                      % plot data, drawn per duck after the loop
warned3 = false;                       % note about missing logger axes, printed once

for i = 1:size(PP, 1)
    fn = PP{i,1};  f1 = PP{i,2};  f2 = PP{i,3};
    V = MagUnN(strcmp({MagUnN.file}, fn));
    S = Sync(strcmp({Sync.file}, fn));
    if isempty(V) || isempty(S)
        fprintf('%-24s SKIPPED: not in MagDistUnN / MagSyncUnN results\n', fn);
        continue
    end

    %% ---- VIDEO: Pt9 acceleration along the body axis ----
    P  = V.P;                              % Pt9 (incl. neck/tail estimate)
    nF = size(P, 1);
    if f2 > nF
        fprintf('%-24s SKIPPED: power phase ends at frame %d, data has %d\n', fn, f2, nF);
        continue
    end
    okP = ~isnan(P(:,1));
    if nnz(okP) < 20
        fprintf('%-24s SKIPPED: too few Pt9 frames\n', fn);
        continue
    end
    [bV, aV] = butter(4, cutHz / (FR/2));
    % filter only from the first to the last frame that has data. (Filling
    % the ends with a constant made the duck "stop" there, which showed up
    % as a huge false deceleration at the end of each recording.)
    seg = find(okP, 1) : find(okP, 1, 'last');
    P(seg, :) = filtfilt(bV, aV, fillmissing(P(seg, :), 'linear'));
    P(~okP, :) = NaN;
    Vv = zeros(nF, 3);  Av = zeros(nF, 3);
    for k = 1:3
        Vv(:,k) = gradient(P(:,k))  .* FR;     % m/s
        Av(:,k) = gradient(Vv(:,k)) .* FR;     % m/s^2
    end

    % body axis = tail base -> neck base (unit vector per frame)
    A = readmatrix(fullfile(accFolder, fn));
    E = A(:, 4:6) - A(:, 7:9);
    if nnz(~isnan(E(:,1))) < 2
        fprintf('%-24s SKIPPED: no neck/tail for the body axis\n', fn);
        continue
    end
    if strcmp(axisMode, 'mean')
        Eu = E ./ vecnorm(E, 2, 2);            % unit vector per frame (NaN where missing)
        Em = mean(Eu, 1, 'omitnan');
        E  = repmat(Em ./ norm(Em), nF, 1);    % the same direction in every frame
    else
        E = fillmissing(E, 'linear', 'EndValues', 'nearest');
        E = smoothdata(E, 'gaussian', axisSmooth);
        E = E ./ vecnorm(E, 2, 2);
    end

    % body frame per frame: forward = E, up = lab Z made perpendicular to E,
    % lateral = up x forward. No roll: only 2 body points are digitized.
    Up = [0 0 1] - sum(E .* [0 0 1], 2) .* E;
    Up = Up ./ vecnorm(Up, 2, 2);
    Lt = cross(Up, E, 2);

    vVid    = sum(Vv .* E,  2);            % forward speed, video (m/s)
    aVid    = sum(Av .* E,  2);            % forward acceleration, video
    aVidLat = sum(Av .* Lt, 2);            % lateral
    aVidUp  = sum(Av .* Up, 2);            % up (body frame)
    w    = (f1 : f2)';                     % power-phase frames
    T    = (f2 - f1) / FR;                 % s
    aVidMean = mean(aVid(w), 'omitnan');
    dvVid    = trapz(aVid(w)) / FR;
    % robust version: straight line through the forward SPEED over the power
    % phase (one derivative only, all frames used). slope = mean acceleration
    okw = ~isnan(vVid(w));
    if nnz(okw) > 3
        pv = polyfit((w(okw) - f1) ./ FR, vVid(w(okw)), 1);
        aVidFit = pv(1);
    else
        aVidFit = NaN;
    end
    dvVidFit = aVidFit * T;
    fracPt9  = mean(V.src(w) == 1);        % share of frames on real Pt9
    nanVid   = nnz(isnan(aVid(w)));

    %% ---- LOGGER: forward acceleration at the same moments ----
    k = find(strcmp(cacheName, S.accFile), 1);
    if isempty(k)
        [~, nm] = fileparts(S.accFile);
        L = load(fullfile(resFolder, [nm '_acce_results.mat']), 'ResOnA');
        cacheName{end+1} = S.accFile;  cacheR{end+1} = L.ResOnA; %#ok<AGROW>
        k = numel(cacheName);
    end
    R = cacheR{k};
    % row-based time, the same clock MagSyncUnN used for tauX
    tA   = (0 : numel(R.Time_s) - 1)' ./ FRacc;
    [bA, aA] = butter(4, cutHz / (FRacc/2));
    aRaw = filtfilt(bA, aA, fillmissing(R.Acceleration(:, 1), 'linear'));  % ForwardRaw, same low-pass as the video
    aFlt = R.Acceleration(:, 2);           % ForwardFiltered (2 Hz), info only

    % lateral and up of the logger, in the same body frame as the video:
    % forward = ForwardAxisSensor, up = gravity at rest (the sensor reads
    % +g upward), lateral = up x forward. Needs LinearAccSensor and
    % GravityStatic in ResOnA (saved by analyze_peaks_and_calibration).
    has3 = isfield(R, 'LinearAccSensor') && isfield(R, 'GravityStatic');
    if has3
        fS = R.ForwardAxisSensor(:)';  fS = fS ./ norm(fS);
        uS = R.GravityStatic(:)';
        uS = uS - dot(uS, fS) .* fS;  uS = uS ./ norm(uS);
        lS = cross(uS, fS);
        aLatRaw = filtfilt(bA, aA, fillmissing(R.LinearAccSensor * lS', 'linear'));
        aUpRaw  = filtfilt(bA, aA, fillmissing(R.LinearAccSensor * uS', 'linear'));
    else
        aLatRaw = nan(size(aRaw));  aUpRaw = aLatRaw;
        if ~warned3
            fprintf(['NOTE: ResOnA has no LinearAccSensor / GravityStatic - accelerometer shows ' ...
                     'X (forward) only. Add them in analyze_peaks_and_calibration.\n']);
            warned3 = true;
        end
    end

    tLog = @(frames) (frames - 1) ./ FR + S.tauX;   % video frame -> logger time
    aImuW  = interp1(tA, aRaw, tLog(w), 'linear');
    aImuMean = mean(aImuW, 'omitnan');
    dvImu    = trapz(aImuW) / FR;
    aImuFltMean = mean(interp1(tA, aFlt, tLog(w), 'linear'), 'omitnan');
    aImuLatMean = mean(interp1(tA, aLatRaw, tLog(w), 'linear'), 'omitnan');
    aImuUpMean  = mean(interp1(tA, aUpRaw,  tLog(w), 'linear'), 'omitnan');

    % sensitivity: logger mean with the sync shifted by +-shiftMax frames
    mShift = nan(size(shifts));
    for q = 1:numel(shifts)
        mShift(q) = mean(interp1(tA, aRaw, tLog(w + shifts(q)), 'linear'), 'omitnan');
    end

    % shape similarity inside the window (information only: ~12 logger
    % samples, and very sensitive to the sync)
    ok = ~isnan(aVid(w)) & ~isnan(aImuW);
    rShape = NaN;
    if nnz(ok) > 5, c = corrcoef(aVid(w(ok)), aImuW(ok)); rShape = c(1,2); end

    %% ---- store ----
    tok = regexp(fn, 'duck(\d+)round(\d+)', 'tokens', 'once');
    C.file       = fn;
    C.duck       = str2double(tok{1});
    C.round      = str2double(tok{2});
    C.accFile    = S.accFile;
    C.frameStart = f1;
    C.frameEnd   = f2;
    C.T_s        = T;
    C.tauX       = S.tauX;
    C.syncCorr   = S.corr;
    C.aVidMean   = aVidMean;     % m/s^2
    C.aImuMean   = aImuMean;     % m/s^2 (raw)
    C.aImuFltMean = aImuFltMean; % m/s^2 (2 Hz filtered) – expect smeared
    C.aVidLat    = mean(aVidLat(w), 'omitnan');   % m/s^2, lateral
    C.aVidUp     = mean(aVidUp(w),  'omitnan');   % m/s^2, up
    C.aImuLat    = aImuLatMean;                   % NaN if not saved by the logger code
    C.aImuUp     = aImuUpMean;
    C.dvVid      = dvVid;        % m/s
    C.dvImu      = dvImu;        % m/s
    C.aVidFit    = aVidFit;      % m/s^2, slope of the video speed (robust)
    C.dvVidFit   = dvVidFit;     % m/s
    C.diff       = aImuMean - aVidMean;
    C.imuShiftMin = min(mShift); % range of the logger mean over the shifts
    C.imuShiftMax = max(mShift);
    C.rShape     = rShape;
    C.fracPt9    = fracPt9;      % 1 = whole window on digitized Pt9
    C.nanVid     = nanVid;
    if isempty(Cmp), Cmp = C; else, Cmp(end+1) = C; end %#ok<AGROW>

    fprintf(['%-24s frames %3d-%3d | video %7.2f | accelerometer %7.2f m/s^2 ' ...
             '(shift +-%d fr: %6.2f..%6.2f) | dv video %5.2f / accelerometer %5.2f m/s | Pt9 %3.0f%%\n'], ...
        fn, f1, f2, aVidMean, aImuMean, shiftMax, C.imuShiftMin, C.imuShiftMax, ...
        dvVidFit, dvImu, 100*fracPt9);
    if fracPt9 < 1
        fprintf('    <-- %.0f%% of the window uses the neck/tail estimate of Pt9\n', 100*(1-fracPt9));
    end
    if nanVid > 0
        fprintf('    <-- %d frames without a video estimate\n', nanVid);
    end

    %% ---- keep what the per-duck plot needs ----
    if doPlot
        fc = (max(1, f1 - round(ctx*FR)) : min(nF, f2 + round(ctx*FR)))';
        iL = tA >= tLog(fc(1)) & tA <= tLog(fc(end));
        % accelerometer speed = running sum (integral) of its acceleration,
        % shifted so that it starts at the video speed at power-phase start
        vAcc  = cumtrapz(tA, aRaw);
        v0 = vVid(f1);
        if isnan(v0), v0 = mean(vVid(w), 'omitnan'); end
        vAccW = vAcc(iL) - interp1(tA, vAcc, tLog(f1)) + v0;
        PL(end+1, :) = {C.duck, C.round, 1000 .* (fc - f1) ./ FR, ...
            [aVid(fc), aVidLat(fc), aVidUp(fc)], 1000 .* (tA(iL) - tLog(f1)), ...
            [aRaw(iL), aLatRaw(iL), aUpRaw(iL)], 1000*T, ...
            [aVidMean, C.aVidLat, C.aVidUp], [aImuMean, C.aImuLat, C.aImuUp], ...
            S.corr, fracPt9, vVid(fc), vAccW, dvVidFit, dvImu}; %#ok<AGROW>
    end
end

if isempty(Cmp), return, end

%% ---- summary ----
x = [Cmp.aVidMean];  y = [Cmp.aImuMean];
ok = ~isnan(x) & ~isnan(y);
fprintf('\nPower phase, %d trials:\n', nnz(ok));
fprintf('   mean acceleration  video %.2f +- %.2f | accelerometer %.2f +- %.2f m/s^2 (mean +- SD)\n', ...
    mean(x(ok)), std(x(ok)), mean(y(ok)), std(y(ok)));
fprintf('   accelerometer - video     %.2f +- %.2f m/s^2\n', mean(y(ok) - x(ok)), std(y(ok) - x(ok)));
if nnz(ok) > 2
    c = corrcoef(x(ok), y(ok));
    pf = polyfit(x(ok), y(ok), 1);
    fprintf('   correlation r = %.2f, accelerometer = %.2f * video %+.2f\n', c(1,2), pf(1), pf(2));
    if c(1,2) < -0.3
        fprintf(['   <-- NEGATIVE correlation: the accelerometer X (forward) axis may be reversed ' ...
                 '(forwardSign in analyze_peaks_and_calibration)\n']);
    end
end
sens = [Cmp.imuShiftMax] - [Cmp.imuShiftMin];
fprintf('   sync sensitivity (+-%d frames): accelerometer mean changes by %.2f m/s^2 on average (max %.2f)\n', ...
    shiftMax, mean(sens), max(sens));
fprintf('   (compare with |accelerometer - video| = %.2f on average: if much smaller, the sync uncertainty does not matter)\n', ...
    mean(abs(y(ok) - x(ok))));

if doPlot
    figure('Name', 'Power phase: video vs accelerometer', 'Color', 'w'); hold on
    cl = [0 0.45 0.74; 0.30 0.65 0.15; 0.90 0.45 0.05];   % forward, lateral, up
    X3 = {x, [Cmp.aVidLat], [Cmp.aVidUp]};
    Y3 = {y, [Cmp.aImuLat], [Cmp.aImuUp]};
    hs = gobjects(1, 3);
    errorbar(x, y, y - [Cmp.imuShiftMin], [Cmp.imuShiftMax] - y, 'o', ...
        'Color', cl(1,:), 'CapSize', 3);          % sync sensitivity, forward only
    for ax = 1:3
        hs(ax) = plot(X3{ax}, Y3{ax}, 'o', 'Color', cl(ax,:), ...
            'MarkerFaceColor', cl(ax,:), 'MarkerSize', 5);
    end
    allv = [X3{:}, Y3{:}, [Cmp.imuShiftMin], [Cmp.imuShiftMax]];
    lim = [min(allv, [], 'omitnan'), max(allv, [], 'omitnan')];
    lim = lim + [-0.1 0.1] .* max(diff(lim), 1);
    plot(lim, lim, 'k--');
    xlim(lim); ylim(lim); axis square; box on
    text(x, y, "  " + erase(string({Cmp.file}), "xyzpts.csv"), 'FontSize', 7);
    xlabel('video: mean acceleration over the power phase, Pt9 (m/s^2)');
    ylabel('accelerometer: mean acceleration over the power phase (m/s^2)');
    legend(hs, {'X (forward)', 'Y (lateral)', 'Z (up)'}, 'Location', 'best');
    title(sprintf('Power phase, mean per axis  (bars = X, sync shifted \\pm%d frames, dashed = 1:1)', shiftMax));
end

%% ---- Plots: ONE FIGURE PER DUCK AND PER AXIS, one panel per trial ----
%   3 separate figures per duck: X (forward), Y (lateral), Z (up), body frame
%   blue o = video Pt9,  red x = accelerometer (raw)   gray = power phase
%   panel title: mean over the power phase on this axis, video | accelerometer
%   red panel title = less than 90% of the window on the digitized Pt9
if doPlot && ~isempty(PL)
    cVid = [0 0.45 0.74];       % video
    cAcc = [0.85 0.10 0.10];    % accelerometer
    axName = {'X (forward)', 'Y (lateral)', 'Z (up)'};
    axShort = {'X', 'Y', 'Z'};
    duckOf = [PL{:, 1}];
    for dk = unique(duckOf)
        q = find(duckOf == dk);
        for ax = 1:3
            figure('Name', sprintf('Power phase duck%d %s', dk, axShort{ax}), ...
                   'Color', 'w', 'WindowState', 'maximized');
            tl = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
            title(tl, sprintf(['duck %d  –  %s acceleration around the power phase (gray):  ' ...
                'o = video Pt9,  x = accelerometer'], dk, axName{ax}));
            for n = q
                [~, rd, tr, av, tg, ag, Tms, mv, ml, sc, fp] = PL{n, 1:11};
                nexttile; hold on
                yv = [av(:, ax); ag(:, ax)];
                yl = [min(yv, [], 'omitnan'), max(yv, [], 'omitnan')];
                if any(isnan(yl)), yl = [-1 1]; end
                yl = yl + [-0.1 0.1] .* max(diff(yl), 1);
                patch([0 Tms Tms 0], yl([1 1 2 2]), [0.92 0.92 0.92], 'EdgeColor', 'none');
                hv = plot(tr, av(:, ax), '-o', 'Color', cVid, 'LineWidth', 1, ...
                    'MarkerSize', 3, 'MarkerIndices', 1:3:numel(tr));
                hl = plot(tg, ag(:, ax), ':x', 'Color', cAcc, 'LineWidth', 1, ...
                    'MarkerSize', 6);
                yline(0, 'k:');
                ylim(yl); box on
                title(sprintf('round %02d: %s video %.2f | accelerometer %.2f (sync %.2f, Pt9 %.0f%%)', ...
                    rd, axShort{ax}, mv(ax), ml(ax), sc, 100*fp), 'FontWeight', 'normal', ...
                    'FontSize', 8, 'Color', [fp < 0.9, 0, 0]);
                if n == q(1)
                    legend([hv, hl], {'video', 'accelerometer'}, 'Location', 'best', 'FontSize', 7);
                end
            end
            xlabel(tl, 'time from power-phase start (ms)');
            ylabel(tl, sprintf('%s acceleration (m/s^2)', axName{ax}));
        end
    end
end

%% ---- Plots: FORWARD SPEED, one figure per duck, one panel per trial ----
%   blue = video speed along X (one derivative only - much less noise)
%   red  = accelerometer speed (integral of its acceleration), started at
%          the video speed at power-phase start (an accelerometer cannot
%          know the speed itself, only how much it changes)
%   if the accelerometer is right, both lines rise and fall together
if doPlot && ~isempty(PL)
    duckOf = [PL{:, 1}];
    for dk = unique(duckOf)
        q = find(duckOf == dk);
        figure('Name', sprintf('Power phase duck%d speed', dk), 'Color', 'w', ...
               'WindowState', 'maximized');
        tl = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, sprintf(['duck %d  –  X (forward) SPEED around the power phase (gray):  ' ...
            'o = video,  x = accelerometer (integrated)'], dk));
        for n = q
            [rd, tr, tg, Tms, sc, fp] = PL{n, [2 3 5 7 10 11]};
            [vv, va, dvv, dva] = PL{n, 12:15};
            nexttile; hold on
            yv = [vv(:); va(:)];
            yl = [min(yv, [], 'omitnan'), max(yv, [], 'omitnan')];
            if any(isnan(yl)), yl = [-1 1]; end
            yl = yl + [-0.1 0.1] .* max(diff(yl), 0.2);
            patch([0 Tms Tms 0], yl([1 1 2 2]), [0.92 0.92 0.92], 'EdgeColor', 'none');
            hv = plot(tr, vv, '-o', 'Color', [0 0.45 0.74], 'LineWidth', 1, ...
                'MarkerSize', 3, 'MarkerIndices', 1:3:numel(tr));
            hl = plot(tg, va, ':x', 'Color', [0.85 0.10 0.10], 'LineWidth', 1, 'MarkerSize', 6);
            ylim(yl); box on
            title(sprintf('round %02d: speed change video %.2f | accelerometer %.2f m/s (sync %.2f, Pt9 %.0f%%)', ...
                rd, dvv, dva, sc, 100*fp), 'FontWeight', 'normal', 'FontSize', 8, ...
                'Color', [fp < 0.9, 0, 0]);
            if n == q(1)
                legend([hv, hl], {'video', 'accelerometer'}, 'Location', 'best', 'FontSize', 7);
            end
        end
        xlabel(tl, 'time from power-phase start (ms)');
        ylabel(tl, 'X (forward) speed (m/s)');
    end

    % summary: speed change over the power phase, one point per trial
    figure('Name', 'Power phase: speed change video vs accelerometer', 'Color', 'w'); hold on
    xd = [Cmp.dvVidFit];  yd = [Cmp.dvImu];
    plot(xd, yd, 'o', 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
    lim = [min([xd yd], [], 'omitnan'), max([xd yd], [], 'omitnan')];
    lim = lim + [-0.1 0.1] .* max(diff(lim), 0.2);
    plot(lim, lim, 'k--');
    xlim(lim); ylim(lim); axis square; box on
    text(xd, yd, "  " + erase(string({Cmp.file}), "xyzpts.csv"), 'FontSize', 7);
    xlabel('video: X speed change over the power phase (m/s)');
    ylabel('accelerometer: X speed change over the power phase (m/s)');
    okd = ~isnan(xd) & ~isnan(yd);
    if nnz(okd) > 2
        c = corrcoef(xd(okd), yd(okd));
        title(sprintf('speed change per power phase (dashed = 1:1),  r = %.2f', c(1,2)));
    else
        title('speed change per power phase (dashed = 1:1)');
    end
end

save(fullfile(resFolder, 'PowerPhaseAccUnN_results.mat'), 'Cmp');
writetable(struct2table(Cmp), fullfile(resFolder, 'PowerPhaseAccUnN_results.csv'));
fprintf('\nSaved: %s (+ .csv)\n', fullfile(resFolder, 'PowerPhaseAccUnN_results.mat'));
end