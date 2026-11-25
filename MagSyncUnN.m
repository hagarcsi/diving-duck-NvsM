% MagSyncUnN  –  video <-> accelerometer sync (UNDERWATER), part 2
%
% Uses the results of MagDistUnN (part 1). For every trial:
%   1. finds the magnetometer peak in the accelerometer data, with
%      SUB-SAMPLE precision (parabola through the samples around the peak;
%      the accelerometer samples only every 10 ms)
%   2. finds the peak of the expected signal (Bgrid) in the video, also
%      sub-frame
%   3. offset between the two clocks:  tAcc = tVid + tau
%        tauPeak = peak to peak (simple)
%        tauX    = whole peak SHAPE (measured vs Bgrid), cross-correlation
%                  around tauPeak  <- main result
%      and, for comparison, what the offset would be if the peak were
%      anchored to the minimum of r (tauMin) or to the 90 deg crossing
%      (tauCross). The differences are printed in ms.
%
% Accelerometer (logger) files, e.g. Yosef1duck3_S1.csv:
%   one long recording per duck (all rounds in one file), 100 Hz,
%   magnetometer compX/Y/Z = cols 10:12. Separator ';' or ',' (a file
%   opened and saved in Excel becomes ',' and its time column loses the
%   milliseconds). Time is therefore ALWAYS built from the row number at
%   FRacc – all files checked so far have no missing samples.
%
% Magnetometer processing:
%   - single-sample spikes removed (sensor glitches)
%   - hard-iron calibration: sphere fitted to the quiet data (no magnet),
%     its centre subtracted, so |B| does not change with the duck's
%     orientation. Kept only if it makes the quiet |B| MORE constant
%     (with few orientations, e.g. duck1, the fit is poorly determined).
%   - |B| of the calibrated vector; baseline = median of the first/last
%     baseEdge seconds of the window, subtracted.
%
% WHICH PASS BELONGS TO WHICH VIDEO TRIAL:
%   Every pass by the magnet is one round, so roundNN = the NN-th pass of
%   that duck, in time order. If a duck has several logger files (duck3:
%   Yosef2 first, then Yosef1), the rounds continue from file to file in
%   the order given in sessionOrder.
%   Weak passes (duck far from the magnet) can fall below passThr: add
%   them by hand in extraPasses. False detections: remove in dropPasses.
%   CHECK: for every trial the code also tries all the other passes of
%   that duck. If another pass matches the video's expected peak shape
%   clearly better, it warns – a sign that the counting is off by one
%   somewhere (a missed or an extra pass).
%
%   AccMap = optional override: logger file + tPass for a trial
%   ('' / [] = by round number).
%
% VALIDATION with the older pipeline (analyze_peaks_and_calibration.m),
% if its saved ResOnA files are listed in oldFiles:
%   1. its peak list vs the passes found here (old-only peaks of ~8200
%      are sensor glitches, which the old findpeaks counted)
%   2. swim direction in each video (+X / -X) vs the notebook direction
%      (East / West) of the pass given to it. One direction should always
%      be the same sign; a trial that breaks this = wrong round/pass.

function [Sync] = MagSyncUnN()

close all

FRacc = 100;        % accelerometer sampling rate, Hz
resFile   = 'C:\Users\hagar\OneDrive\Desktop\for Phd!\question3\diving\MagDistUnN_results.mat';
accFolder = 'C:\Users\hagar\OneDrive\Desktop\for Phd!\נטה אדומת ראש\שחייה מתחת למים\20260706\נתוני תאוצה';
outFolder = fileparts(resFile);
magCols   = 10:12;  % compX, compY, compZ

% order of the logger files for a duck with more than one file
% (rounds continue from file to file). Ducks not listed: one file.
sessionOrder = {
    3, {'Yosef2duck3_S1.csv', 'Yosef1duck3_S1.csv'}    % started with Yosef2
};

% rounds in the notebook with NO pass by the magnet (the duck did not
% swim past it), per duck. Later rounds then move up by one pass.
%   e.g.  noPassRounds = {3, [4]};   -> duck3 round 4 had no pass
% duck3 round 5 = first round of the 2nd session (Yosef1duck3): it has a
% video (not digitized) but no pass was detected in the logger - probably
% a weak pass below passThr. Session 2 = 7 rounds (5-11) but 6 passes.
% Without this, rounds 6-11 got the pass of the NEXT round (video swim
% direction did not match the notebook in rounds 6, 7, 8, 11).
noPassRounds = {3, [5]};

% passes too weak to be detected (logger file, approximate time in s)
% NOTE: the weak, wide events at ~260-340 s (duck3 Yosef2 335 s, duck4
% 260 s, duck5 284 s) fall right after the calibration and look like
% handling, not a swim – not added.
extraPasses = cell(0, 2);
% detections that are NOT a pass (logger file, approximate time in s),
% e.g.  dropPasses = {'Yosef2duck4_S1.csv', 1003.6};
dropPasses = cell(0, 2);

% VALIDATION with the older pipeline (analyze_peaks_and_calibration.m):
% its saved results (ResOnA, one .mat per logger file) hold the peak
% times AND the swim direction from the notebook (East/West). Used to
%   1. compare its peak list with the passes found here
%   2. check that the swim direction in every video matches the notebook
%      direction of the pass it was given (a shifted round shows up here)
% Found automatically: every *_acce_results.mat in accResFolder
% (analyze_peaks_and_calibration saves one per logger file there).
% The older duck1 file (duck_acce_calibration_results.mat) is no longer
% used: Yosef1duck1_S1_acce_results.mat replaces it (same file, newer run,
% sensor glitches removed). Both together would set duck1's directions twice.
accResFolder = 'C:\Users\hagar\OneDrive\Desktop\for Phd!\question3\diving';
dOld = dir(fullfile(accResFolder, '*_acce_results.mat'));
dOld = dOld(~strcmp({dOld.name}, 'all_acce_results.mat'));   % combined file (ResAll)
oldFiles = fullfile(accResFolder, {dOld.name})';
if isempty(oldFiles)
    fprintf('No *_acce_results.mat in %s - run analyze_peaks_and_calibration first\n', accResFolder);
end

% Pass for every video trial, from round_to_imu_peak.csv (notebook film
% numbers -> logger peak). Needed because the logger also has passes with
% NO video (the "X" rows in the notebook: e.g. duck1's first pass, duck4
% round 2's first pass), so "roundNN = NN-th pass" shifts the rounds.
% Checked: the swim direction in every video (+X/-X) matches this list
% for all 23 trials; the round-number method gave a different pass in 14.
% duck3, Yosef2 session: 4 films (confirmed) = 4 passes, rounds 01-04.
% The csv assumed 5 (incl. the weak event at 355.8 s = handling, not a
% pass), so for round01/03 it is one pass early - corrected here.
% Rounds 06-11 = Yosef1duck3 passes 1-6 (as in the csv); round05 = a
% video with no pass.
% Leave '' / [] to use the round number instead.
AccMap = {
% video file                  logger file              tPass (s)   notebook
'duck1round01xyzpts.csv'      'Yosef1duck1_S1.csv'      428.37   % R1F1, E->W
'duck1round02xyzpts.csv'      'Yosef1duck1_S1.csv'      477.98   % R1F2, W->E
'duck1round04xyzpts.csv'      'Yosef1duck1_S1.csv'      556.57   % R1F4, W->E
'duck1round07xyzpts.csv'      'Yosef1duck1_S1.csv'     1108.52   % R2F1, W->E
'duck1round08xyzpts.csv'      'Yosef1duck1_S1.csv'     1149.11   % R2F2, E->W
'duck1round09xyzpts.csv'      'Yosef1duck1_S1.csv'     1174.36   % R2F3, W->E
'duck1round12xyzpts.csv'      'Yosef1duck1_S1.csv'     1258.06   % R2F6, E->W
'duck3round01xyzpts.csv'      'Yosef2duck3_S1.csv'      512.66   % R1F1, W->E
'duck3round03xyzpts.csv'      'Yosef2duck3_S1.csv'     1056.40   % R1F3, W->E
'duck3round06xyzpts.csv'      'Yosef1duck3_S1.csv'      327.44   % R2F1, W->E
'duck3round07xyzpts.csv'      'Yosef1duck3_S1.csv'      485.43   % R2F2, E->W
'duck3round08xyzpts.csv'      'Yosef1duck3_S1.csv'      566.37   % R2F3, W->E
'duck3round11xyzpts.csv'      'Yosef1duck3_S1.csv'      851.39   % R2F6, E->W
'duck4round01xyzpts.csv'      'Yosef2duck4_S1.csv'      299.92   % R1F1, W->E
'duck4round02xyzpts.csv'      'Yosef2duck4_S1.csv'      376.26   % R1F2, E->W
'duck4round03xyzpts.csv'      'Yosef2duck4_S1.csv'      406.90   % R1F3, W->E
'duck4round05xyzpts.csv'      'Yosef2duck4_S1.csv'      482.66   % R1F5, W->E
'duck4round10xyzpts.csv'      'Yosef2duck4_S1.csv'     1706.62   % R2F4, W->E
'duck5round01xyzpts.csv'      'Yosef1duck5_S1.csv'      340.88   % R1F1, W->E
'duck5round02xyzpts.csv'      'Yosef1duck5_S1.csv'      836.84   % R1F2, W->E
'duck5round06xyzpts.csv'      'Yosef1duck5_S1.csv'     1071.78   % R1F6, W->E
'duck5round17xyzpts.csv'      'Yosef1duck5_S1.csv'     2653.06   % R3F5, E->W
'duck5round18xyzpts.csv'      'Yosef1duck5_S1.csv'     2691.09   % R3F6, W->E
};

prm.winHalf  = 1.5;     % s, window = tPass +- winHalf (a pass is < 1 s wide)
prm.baseEdge = 0.3;     % s at each end of the window used for the baseline
prm.nPk      = 4;       % samples on each side of the peak for the parabola
prm.maxLag   = 0.15;    % s, cross-correlation search range around tauPeak
prm.lagStep  = 0.0005;  % s, cross-correlation step (0.5 ms)
minCorr  = 0.90;        % warn if the shape match is worse than this
betterBy = 0.05;        % warn if another pass matches better by this much

% pass detection
passThr  = 100;     % |B| above the local baseline to count as a pass
maxBase  = 300;     % ignore events on a baseline above this (not the
                    % magnet: e.g. the logger taken out at the end)
calMax   = 200;     % |B| below this = quiet data for the calibration
doPlot   = true;

load(resFile, 'MagUnN');

%% ---- Read + calibrate every logger file once, list passes ----
d = dir(fullfile(accFolder, '*.csv'));
d = d(~contains({d.name}, 'xyzpts'));        % skip the video digitizations
names = unique([{d.name}'; AccMap(~cellfun(@isempty, AccMap(:,2)), 2)]);
Acc = struct('file', {}, 'duck', {}, 't', {}, 'B', {}, 'passes', {}, 'dir', {});
for k = 1:numel(names)
    fA = fullfile(accFolder, names{k});
    if ~isfile(fA)
        fprintf('Logger file not found: %s\n', names{k});
        continue
    end
    try
        Braw = readAcc(fA, magCols);
    catch
        Braw = [];
    end
    if size(Braw, 1) < 1000
        fprintf('Not a logger file (skipped): %s\n', names{k});
        continue
    end
    tA = (0:size(Braw,1)-1)' ./ FRacc;
    Braw = despike(Braw);
    [off, R] = hardIron(Braw, calMax);
    mRaw = vecnorm(Braw, 2, 2);
    Bc   = vecnorm(Braw - off, 2, 2);
    q = mRaw < calMax;
    sRaw = diff(prctile(mRaw(q), [5 95]));
    sCal = diff(prctile(Bc(q),   [5 95]));
    calTxt = sprintf('calibrated (sphere radius %.0f)', R);
    if sCal >= sRaw                           % calibration does not help
        Bc = mRaw;  off = [0 0 0];
        calTxt = 'NOT calibrated (fit did not help)';
    end
    P = findPasses(tA, Bc, FRacc, passThr, maxBase);

    % hand-added weak passes / removed false ones
    for e = find(strcmp(extraPasses(:,1)', names{k}))
        P(end+1, :) = [extraPasses{e,2}, NaN, NaN]; %#ok<AGROW>
    end
    for e = find(strcmp(dropPasses(:,1)', names{k}))
        P(abs(P(:,1) - dropPasses{e,2}) < 1, :) = [];
    end
    P = sortrows(P, 1);

    tok = regexp(names{k}, 'duck(\d+)', 'tokens', 'once');
    duck = NaN;  if ~isempty(tok), duck = str2double(tok{1}); end
    Acc(end+1) = struct('file', names{k}, 'duck', duck, 't', tA, 'B', Bc, ...
                        'passes', P, 'dir', strings(size(P, 1), 1)); %#ok<AGROW>
    fprintf('%s (duck %d): %.1f min, %d passes, %s, quiet |B| spread %.0f -> %.0f\n', ...
        names{k}, duck, tA(end)/60, size(P, 1), calTxt, sRaw, min(sRaw, sCal));
end

%% ---- Validation 1: compare with the older pipeline's peak lists ----
for f = 1:numel(oldFiles)
    if ~isfile(oldFiles{f})
        fprintf('\nOld results not found (skipped): %s\n', oldFiles{f});
        continue
    end
    O = load(oldFiles{f});
    if ~isfield(O, 'ResOnA'), continue, end
    for r = 1:numel(O.ResOnA)
        R = O.ResOnA(r);
        [~, nm] = fileparts(char(string(R.FileName)));
        a = find(arrayfun(@(x) startsWith(erase(x.file, '.csv'), nm) || ...
                               startsWith(nm, erase(x.file, '.csv')), Acc), 1);
        if isempty(a)
            fprintf('\nOld results for %s: no matching logger file read here\n', nm);
            continue
        end
        T = R.Peaks;
        tOld = T.PeakTime_s(:);  vOld = T.MaxVecMag(:);  dOld = string(T.Direction(:));
        P = Acc(a).passes(:,1);
        used = false(size(P));
        fprintf('\nValidation vs older pipeline – %s: %d old peaks, %d passes here\n', ...
            Acc(a).file, numel(tOld), numel(P));
        for j = 1:numel(tOld)
            [dd, ip] = min(abs(P - tOld(j)));
            if ~isempty(dd) && dd < 0.5
                used(ip) = true;
                Acc(a).dir(ip) = dOld(j);
                fprintf('   old #%2d %7.1f s %-5s = pass at %7.1f s\n', j, tOld(j), dOld(j), P(ip));
            else
                why = '';  if vOld(j) > 5000, why = '  (sensor glitch, ~8200)'; end
                fprintf('   old #%2d %7.1f s %-5s   <-- NOT a pass here%s\n', j, tOld(j), dOld(j), why);
            end
        end
        for ip = find(~used)'
            fprintf('   pass at %7.1f s   <-- not in the old list\n', P(ip));
        end
    end
end

%% ---- Number the passes of every duck = round numbers ----
% Rounds{duck} = [logger index, tPass, height, width] , row r = round r
Rounds = {};
dk = unique([Acc.duck]);
for duck = dk(~isnan(dk))
    aIdx = find([Acc.duck] == duck);
    s = find(cellfun(@(x) isequal(x, duck), sessionOrder(:,1)), 1);
    if ~isempty(s)                              % files in session order
        [ok, pos] = ismember(sessionOrder{s,2}, {Acc.file});
        aIdx = pos(ok);
    end
    L = zeros(0, 4);
    for a = aIdx
        P = Acc(a).passes;
        L = [L; [repmat(a, size(P, 1), 1), P]]; %#ok<AGROW>
    end
    Rounds{duck} = L; %#ok<AGROW>
    fprintf('\nduck %d: %d passes\n', duck, size(L, 1));
    for r = 1:size(L, 1)
        hTxt = sprintf('height %5.0f  width %3.0f ms', L(r,3), 1000*L(r,4));
        if isnan(L(r,3)), hTxt = 'added by hand (extraPasses)'; end
        fprintf('   pass %2d: %-22s %7.1f s (%02d:%05.2f)  %s\n', r, Acc(L(r,1)).file, ...
            L(r,2), floor(L(r,2)/60), mod(L(r,2), 60), hTxt);
    end
end
fprintf('\n');

%% ---- Loop over trials: pass = round number, then sync ----
problems = cell(0, 2);
Sync = [];
for i = 1:numel(MagUnN)
    V = MagUnN(i);
    tok = regexp(V.file, 'duck(\d+)round(\d+)', 'tokens', 'once');
    duck = str2double(tok{1});  round_ = str2double(tok{2});
    k = find(strcmp(AccMap(:,1), V.file), 1);
    fOver = '';  tOver = [];
    if ~isempty(k), fOver = AccMap{k,2};  tOver = AccMap{k,3}; end
    L = [];
    if duck <= numel(Rounds), L = Rounds{duck}; end

    % the pass for this trial
    if ~isempty(fOver) && ~isempty(tOver)               % override
        a = find(strcmp({Acc.file}, fOver), 1);
        if isempty(a)
            problems(end+1, :) = {V.file, ['logger file not read: ' fOver]}; %#ok<AGROW>
            fprintf('%-24s SKIPPED: logger file not read: %s\n', V.file, fOver);
            continue
        end
        tP = tOver;  how = 'AccMap';
    else
        s = find(cellfun(@(x) isequal(x, duck), noPassRounds(:,1)), 1);
        skip = [];  if ~isempty(s), skip = noPassRounds{s,2}; end
        if ismember(round_, skip)
            problems(end+1, :) = {V.file, 'round listed in noPassRounds'}; %#ok<AGROW>
            fprintf('%-24s SKIPPED: round %d is in noPassRounds\n', V.file, round_);
            continue
        end
        iPass = round_ - nnz(skip < round_);            % pass number
    end
    if isempty(fOver) || isempty(tOver)
        if iPass <= size(L, 1)                          % by round number
            a = L(iPass, 1);  tP = L(iPass, 2);
            how = sprintf('round %d = pass %d', round_, iPass);
        else
            problems(end+1, :) = {V.file, sprintf('pass %d, but duck %d has %d passes', ...
                iPass, duck, size(L, 1))}; %#ok<AGROW>
            fprintf('%-24s SKIPPED: needs pass %d, but duck %d has only %d passes\n', ...
                V.file, iPass, duck, size(L, 1));
            continue
        end
    end
    S = syncOne(Acc(a).t, Acc(a).B, tP, V, prm);
    if isempty(S)
        problems(end+1, :) = {V.file, 'no usable window'}; %#ok<AGROW>
        fprintf('%-24s SKIPPED: no usable window\n', V.file);
        continue
    end

    % CHECK: do other passes of this duck match the video better?
    S.betterRound = NaN;  S.betterCorr = NaN;
    for r = 1:size(L, 1)
        if L(r,1) == a && abs(L(r,2) - tP) < 1, continue, end
        Rr = syncOne(Acc(L(r,1)).t, Acc(L(r,1)).B, L(r,2), V, prm);
        if ~isempty(Rr) && Rr.corr > S.corr + betterBy && ...
                (isnan(S.betterCorr) || Rr.corr > S.betterCorr)
            S.betterRound = r;  S.betterCorr = Rr.corr;
        end
    end

    S.file    = V.file;
    S.duck    = duck;
    S.round   = round_;
    S.accFile = Acc(a).file;
    S.tPass   = tP;
    % swim direction in the video (+1 = toward +X) and in the notebook
    S.vidDir  = sign(V.P(V.valid(2), 1) - V.P(V.valid(1), 1));
    [dIp, ip] = min(abs(Acc(a).passes(:,1) - tP));   % nearest detected pass
    if isempty(dIp) || dIp > 0.5, ip = []; end
    S.noteDir = "";  if ~isempty(ip), S.noteDir = Acc(a).dir(ip); end
    if isempty(Sync), Sync = S; else, Sync(end+1) = orderfields(S, Sync(1)); end %#ok<AGROW>

    % --- report ---
    fprintf('%-24s %s, %s @ %7.1f s | tau = %9.4f s | corr %.3f\n', V.file, ...
        how, S.accFile, S.tPass, S.tauX, S.corr);
    if ~isnan(S.betterRound)
        fprintf(['    <-- pass %d (%.1f s) matches this video better (corr %.3f): ' ...
                 'a round without a pass, or an extra pass, before it?\n'], S.betterRound, ...
                 L(S.betterRound, 2), S.betterCorr);
    end
    dTxt = {'-X', '?', '+X'};
    if S.noteDir ~= "" && S.noteDir ~= "?"
        fprintf('    video swims toward %s | notebook: %s\n', dTxt{S.vidDir + 2}, S.noteDir);
    end
    fprintf('    vs peak-to-peak %+5.1f ms | min r %+5.1f ms | crossing %+5.1f ms\n', ...
        1000*(S.tauPeak - S.tauX), 1000*(S.tauMin - S.tauX), 1000*(S.tauCross - S.tauX));
    if S.pkSign < 0,  fprintf('    <-- magnetometer shows a DIP, not a peak\n'); end
    if S.atEdge,      fprintf('    <-- peak at the edge of the window\n'); end
    if S.lagEdge,     fprintf('    <-- best lag at the edge of the search range, check!\n'); end
    if S.corr < minCorr
        fprintf('    <-- poor shape match (%.2f): check the pass\n', S.corr);
    end
end

%% ---- Validation 2: video swim direction vs notebook direction ----
if ~isempty(Sync) && any([Sync.noteDir] ~= "")
    nd = [Sync.noteDir];  vd = [Sync.vidDir];
    has = nd ~= "" & nd ~= "?";             % "?" = direction not filled yet
    fprintf('\nSwim direction, video vs notebook (%d trials with a notebook direction):\n', nnz(has));
    dTxt = {'-X', '?', '+X'};
    for d = unique(nd(has))
        sgn = vd(has & nd == d);
        maj = sign(sum(sgn));  if maj == 0, maj = sgn(1); end
        fprintf('   %-5s = video %s in %d of %d trials\n', d, dTxt{maj + 2}, nnz(sgn == maj), numel(sgn));
        for q = find(has & nd == d & vd ~= maj)
            fprintf('      <-- %s: video %s but notebook %s – wrong round/pass?\n', ...
                Sync(q).file, dTxt{vd(q) + 2}, d);
        end
    end
end

%% ---- Checks: same pass twice, order vs round number ----
if ~isempty(Sync)
    key = strcat({Sync.accFile}, '@', arrayfun(@(s) sprintf('%.1f', s.tPass), Sync, 'Uni', 0));
    [u, ~, j] = unique(key);
    for q = find(accumarray(j(:), 1) > 1)'
        fprintf('\n<-- SAME PASS for %s: %s\n', strjoin({Sync(j == q).file}, ', '), u{q});
    end
    for a = unique({Sync.accFile})
        s = Sync(strcmp({Sync.accFile}, a{1}));
        [~, o] = sort([s.round]);
        if any(diff([s(o).tPass]) <= 0)
            fprintf('\n<-- ORDER: in %s the passes do not follow the round numbers:\n', a{1});
            for q = o
                fprintf('      round %2d -> %7.1f s\n', s(q).round, s(q).tPass);
            end
        end
    end
end

if ~isempty(problems)
    fprintf('\n%d trials SKIPPED:\n', size(problems, 1));
    for k = 1:size(problems, 1)
        fprintf('   %-24s %s\n', problems{k, :});
    end
end

%% ---- Plots ----
if doPlot && ~isempty(Sync)
    for i = 1:numel(Sync)
        S = Sync(i);
        V = MagUnN(strcmp({MagUnN.file}, S.file));
        figure('Name', S.file, 'Color', 'w');

        % A: measured vs expected, on the VIDEO clock after the shift
        subplot(2,1,1); hold on
        plot(S.ta - S.tauX, S.dev ./ max(S.dev), 'k.-', 'MarkerSize', 10);
        plot(V.t, V.Bgrid ./ max(V.Bgrid), '-', 'Color', [0 0.45 0.75], 'LineWidth', 1.5);
        xline(S.tPkVid, 'b:');
        xline(V.tMin,   'r--');
        xline(V.tCross, 'm--');
        xlim(S.tPkVid + [-1 1]); box on
        xlabel('video time (s)'); ylabel('normalized');
        title({strrep(S.file, '_', '\_'), ...
            sprintf('A. %s @ %.1f s (black, shifted by %.4f s) vs expected (blue), corr %.3f', ...
            strrep(S.accFile, '_', '\_'), S.tPass, S.tauX, S.corr)});
        legend({'magnetometer', 'Bgrid (video)', 'Bgrid peak', 'min r', '90 deg crossing'}, ...
            'Location', 'best');

        % B: correlation vs lag
        subplot(2,1,2); hold on
        plot(1000*(S.lags - S.tauX), S.c, 'k-', 'LineWidth', 1.2);
        xline(0, 'g-');
        xline(1000*(S.tauPeak - S.tauX), 'b:');
        xlabel('lag relative to best (ms)'); ylabel('correlation'); box on
        title('B. Shape match vs shift  (green = best, blue = peak-to-peak)');
    end

    % summary: differences between the anchoring methods
    figure('Name', 'Sync summary', 'Color', 'w'); hold on
    n = numel(Sync);
    plot(1:n, 1000*([Sync.tauPeak]  - [Sync.tauX]), 'bo', 'MarkerFaceColor', 'b');
    plot(1:n, 1000*([Sync.tauMin]   - [Sync.tauX]), 'r^', 'MarkerFaceColor', 'r');
    plot(1:n, 1000*([Sync.tauCross] - [Sync.tauX]), 'ms', 'MarkerFaceColor', 'm');
    yline(0, 'k-');
    yline([-10 10], 'k:');                 % one accelerometer sample
    set(gca, 'XTick', 1:n, 'XTickLabel', strrep({Sync.file}, 'xyzpts.csv', ''));
    xtickangle(60); box on
    ylabel('difference from shape match (ms)');
    title('Anchoring methods vs shape match  (dotted = \pm1 accelerometer sample)');
    legend({'peak to peak', 'min r', '90 deg crossing'}, 'Location', 'best');
end

save(fullfile(outFolder, 'MagSyncUnN_results.mat'), 'Sync', 'problems');
end

%% ---- helper: sync one video trial to one pass ----
function R = syncOne(tA, Bc, tPass, V, prm)
R = [];
in = tA >= tPass - prm.winHalf & tA <= tPass + prm.winHalf & ~isnan(Bc);
ta = tA(in);  bm = Bc(in);
if numel(ta) < 20, return, end

% baseline and peak
edge = ta < ta(1) + prm.baseEdge | ta > ta(end) - prm.baseEdge;
base = median(bm(edge));
dev  = bm - base;
[~, ip] = max(abs(dev));
pkSign = sign(dev(ip));                 % +1 = peak, -1 = dip
dev = pkSign .* dev;                    % peak always upward from here
tPkAcc = subPeak(ta, dev, ip, prm.nPk);

% expected signal (video), sub-frame peak
tv  = V.t;
mdl = V.Bgrid;
okv = ~isnan(mdl);
[~, iv] = max(mdl);
tPkVid = subPeak(tv, mdl, iv, prm.nPk);
tauPeak = tPkAcc - tPkVid;

% shape match: Bgrid shifted by tau vs measured, correlation
lags = tauPeak + (-prm.maxLag : prm.lagStep : prm.maxLag);
c = nan(size(lags));
for j = 1:numel(lags)
    m = interp1(tv(okv), mdl(okv), ta - lags(j));   % model at acc times
    ok = ~isnan(m);
    if nnz(ok) < 10, continue, end
    Rc = corrcoef(m(ok), dev(ok));
    c(j) = Rc(1,2);
end
if all(isnan(c)), return, end
[cMax, jm] = max(c);

R.tauX     = lags(jm);        % s, main result: tAcc = tVid + tauX
R.corr     = cMax;            % shape match (1 = perfect)
R.tauPeak  = tauPeak;         % s, peak to peak
R.tauMin   = tPkAcc - V.tMin;     % s, if the peak were at min r
R.tauCross = tPkAcc - V.tCross;   % s, if the peak were at the crossing
R.tPkAcc   = tPkAcc;          % s, magnetometer peak (logger clock)
R.tPkVid   = tPkVid;          % s, Bgrid peak (video clock)
R.pkSign   = pkSign;          % +1 peak, -1 dip
R.base     = base;            % baseline |B|
R.height   = max(dev);        % peak height above baseline
R.atEdge   = ip <= prm.nPk || ip > numel(ta) - prm.nPk;
R.lagEdge  = (jm == 1 || jm == numel(lags));
R.ta       = ta;              % logger time in window
R.dev      = dev;             % |B| - baseline (sign-corrected)
R.lags     = lags;
R.c        = c;
end

%% ---- helper: read the magnetometer columns (';' or ',' separated) ----
function B = readAcc(f, magCols)
fid = fopen(f, 'r', 'n', 'UTF-8');
hdr = fgetl(fid);
if contains(hdr, ';'), dl = ';'; else, dl = ','; end
nCol = max(magCols);
fmt = [repmat('%s', 1, nCol) '%*[^\n]'];
C = textscan(fid, fmt, 'Delimiter', dl);
fclose(fid);
B = zeros(numel(C{1}), numel(magCols));
for k = 1:numel(magCols)
    B(:,k) = str2double(C{magCols(k)});
end
bad = all(B == 0, 2) | any(isnan(B), 2);   % e.g. last row of zeros
B(bad,:) = [];
end

%% ---- helper: remove single-sample spikes (sensor glitches) ----
function B = despike(B)
m = vecnorm(B, 2, 2);
dm = abs(m - movmedian(m, 5));
sp = dm > 5 * max(median(dm), 1) & dm > 200;
for k = 1:size(B, 2)
    mk = movmedian(B(:,k), 5);
    B(sp, k) = mk(sp);
end
end

%% ---- helper: hard-iron offset = centre of a sphere fitted to quiet data ----
function [off, R] = hardIron(B, calMax)
X = B(vecnorm(B, 2, 2) < calMax, :);
A = [2 * X, ones(size(X, 1), 1)];
p = A \ sum(X.^2, 2);
off = p(1:3)';
R = sqrt(p(4) + off * off');
end

%% ---- helper: list magnet passes [tPeak height width] ----
function P = findPasses(t, m, FR, thr, maxBase)
base = movmedian(m, round(3 * FR));
dev = m - base;
hi = dev > thr & base < maxBase;
P = zeros(0, 3);
d = diff([0; hi; 0]);
s = find(d == 1);  e = find(d == -1) - 1;
k = 1;
while k <= numel(s)                                % merge events < 1 s apart
    s1 = s(k);  e1 = e(k);
    while k < numel(s) && t(s(k+1)) - t(e1) < 1
        k = k + 1;  e1 = e(k);
    end
    seg = s1 : e1;
    [h, j] = max(dev(seg));
    half = seg(dev(seg) > h / 2);
    P(end+1, :) = [t(seg(j)), h, t(half(end)) - t(half(1))]; %#ok<AGROW>
    k = k + 1;
end
end

%% ---- helper: sub-sample peak (parabola through 2*n+1 points) ----
function tp = subPeak(t, y, i, n)
tp = t(i);
idx = max(1, i-n) : min(numel(y), i+n);
idx = idx(~isnan(y(idx)));
if numel(idx) < 3 || i == idx(1) || i == idx(end), return, end
p = polyfit(t(idx) - t(i), y(idx), 2);
if p(1) >= 0, return, end                  % not a peak shape
tp = t(i) - p(2) / (2*p(1));
tp = min(max(tp, t(idx(1))), t(idx(end))); % stay inside the points used
end