% MagDistUnN  –  distance between the accelerometer and the magnet (UNDERWATER)
%
% Part 1 of the video <-> accelerometer sync:
%   For every digitized trial, computes the distance r(t) between the
%   accelerometer (Pt9) and the magnet centre, and the expected
%   magnetometer signal shape 1/r^3.
%
%   Pt2 cols  4: 6 = Neck base
%   Pt3 cols  7: 9 = Tail base
%   Pt9 cols 25:27 = estimated accelerometer position on the duck
%   Pt10 cols 28:30, Pt11 cols 31:33 = OPTIONAL, two sharp points on the
%                    logger itself (e.g. its front and back corner).
%                    Where BOTH are digitized, their midpoint replaces Pt9
%                    (less digitizing noise; same acceleration, the logger
%                    is rigid). It is shifted by the mean offset Pt9 -
%                    midpoint over the frames that have both, so the
%                    position (and the distance to the magnet) stays on Pt9.
%                    The distance Pt10-Pt11 must be constant: its SD is
%                    printed as the digitizing precision (devInfo).
%   Magnet file     = 2 corners of the magnet (pt1, pt2), same for all
%                     underwater trials (cameras did not move)
%
% Extending Pt9 outside the frames where it was digitized:
%   The accelerometer is on the body, so its offset from the neck / tail
%   base is roughly constant over a short time. For each trial
%   separately, the 3D offset (Pt9 - ref) is averaged over the nOff
%   frames where both are digitized and that are CLOSEST to the gap
%   (first frames for the gap before Pt9, last frames for the gap after).
%   The offset is a 3D vector, so swim direction is already in it.
%   Neck is used first (closer to the accelerometer = smaller error),
%   tail base only where the neck is missing.
%
%   src (per frame): 0 = no data, 1 = Pt9 (or Pt10/Pt11 midpoint, see
%                    fromDev) digitized (or short gap filled),
%                    2 = estimated from neck, 3 = estimated from tail base
%
% Magnet:
%   A plate standing perpendicular to the floor (W x H x T), magnetized
%   through its thickness: axis = horizontal, normal to the broad face.
%   The two digitized corners = the TOP-BACK edge (along the width).
%   r and s are measured from the midpoint of that edge (magC).
%
% Expected signal (all |B| up to a constant):
%   invr3 = 1/r^3                       (distance only, from magC)
%   Bdip  = point dipole at the magnet's volume centre
%   Bgrid = the whole plate as nW x nH x nT small dipoles, summed.
%           Best model at 10-30 cm, where the magnet is not a "point":
%           the top of the magnet (closer to the duck) weighs more.
%   (assumes Z = vertical)
%
% Anchor point – the "90 deg crossing":
%   s(t) = position along the swim path relative to the magnet
%          (projection of magnet->accelerometer on the path direction).
%   s = 0 when the vector to the magnet is at 90 deg to the path.
%   s is ~linear in time, so the crossing frame is found by a linear fit
%   of nFit frames, inside the data or EXTRAPOLATED outside it (can be a
%   negative frame, or after the last frame).
%   Tested by truncating 5 trials 0-20 frames after the true crossing
%   (pathMode 'X', nFit 35): mean error 1.1, max 2.6 frames.
%
% NOTE: coordinates are NOT flipped by swim direction (flg) here.
%       The magnet is in the raw camera frame, so the duck must stay in
%       the same frame.

function [MagUnN] = MagDistUnN()

close all

FR = 250;           % filming rate, fps
maxGap = 5;         % fill NaN gaps up to this many frames (linear)
smoothWin = 7;      % gaussian smoothing window, frames (~28 ms) ; 0 = off
nOff = 7;           % frames (nearest the gap) used for the neck/tail offset
                    % small = the extension joins Pt9 without a jump.
                    % (15 gave jumps of up to 4 cm at the seam, because
                    %  the offset slowly changes as the body moves)
maxSeam = 3;        % warn if the jump at a seam > maxSeam x normal step
useRefs = true;     % extend Pt9 using neck / tail base
useDev  = true;     % use Pt10 + Pt11 (2 sharp points on the logger) where digitized
devCols = {28:30, 31:33};   % Pt10, Pt11 columns in the xyzpts file
edgeFrames = 3;     % min this close to the start/end of the data = suspicious
nFit = 35;          % frames for the linear fit of s(t) (90 deg crossing)
                    % (tested: 15 -> mean err 2.9 / max 7.2 frames,
                    %          35 -> mean err 1.1 / max 2.6 frames)
pathMode = 'wall';  % direction for the 90 deg crossing:
                    %  'wall'   = along the magnet's face (parallel to the
                    %             wall) -> crossing = when the accelerometer
                    %             passes the magnet's centre line, normal to
                    %             its face. Follows the magnet even if the
                    %             wall is tilted vs the lab X axis
                    %  'X'      = lab X axis -> crossing = when the
                    %             accelerometer passes the magnet's X.
                    %             Same definition for all trials
                    %  'global' = main direction of the whole trajectory
                    %  'local'  = direction of motion near the crossing
                    %  (the crossing frame is sensitive to this choice,
                    %   ~1 frame per degree at r = 0.2 m)
dirWin = 30;        % frames nearest the crossing ('local' mode only)
maxExtrap = 40;     % warn if the crossing is extrapolated further than this
doPlot = true;

% magnet size (m) and grid for the field model
magW = 0.100;       % width (measured; digitized corners give 9.3 cm)
                    % [] = use the distance between the digitized corners
magH = 0.145;       % height (measured)
magT = 0.025;       % thickness (measured)
nGrid = [10 10 3];  % dipoles along width, height, thickness

% direction of the magnet's broad face (horizontal), used for its axis:
%   'corners' = from the 2 digitized corners. NOTE: only ~10 cm apart and
%               Y (depth from the cameras) is the noisiest axis (~3x X),
%               so the angle can be off by 10-25 deg
%   'wall'    = parallel to wallDir (the magnet was placed parallel to the
%               pool wall). Best: digitize 2 points FAR apart (>50 cm)
%               along a straight edge parallel to the wall, and put
%               their difference here
magFaceMode = 'wall';
wallFile = 'magnetptsxyzpts.csv';  % DLTdv file (in magFolder) with 2 points
                    %     along a line parallel to the wall. '' = use wallDir
wallCols = 7:12;    % columns of those 2 points (here pt3, pt4 = the window)
frontCols = 13:15;  % top FRONT corner of the magnet (pt5), in magFile.
                    % [] = unknown (thickness centred on the back edge)
wallDir = [1 0 0];  % used only if wallFile is ''

folder    = 'C:\Users\hagar\OneDrive\Desktop\for Phd!\נטה אדומת ראש\שחייה מתחת למים\20260706\נתוני תאוצה';   % duck digitizations (with Pt9)
magFolder = 'C:\Users\hagar\OneDrive\Desktop\for Phd!\question3\diving';   % magnet digitization
outFolder = magFolder;                                                     % where results are saved
magFile   = 'magnetptsxyzpts.csv';

FN = {'duck1round01xyzpts.csv'
'duck1round02xyzpts.csv'
'duck1round04xyzpts.csv'
'duck1round07xyzpts.csv'
'duck1round08xyzpts.csv'
'duck1round09xyzpts.csv'
'duck1round12xyzpts.csv'
'duck3round01xyzpts.csv'
'duck3round03xyzpts.csv'
'duck3round06xyzpts.csv'
'duck3round07xyzpts.csv'
'duck3round08xyzpts.csv'
'duck3round11xyzpts.csv'
'duck4round01xyzpts.csv'
'duck4round02xyzpts.csv'
'duck4round03xyzpts.csv'
'duck4round05xyzpts.csv'
'duck4round10xyzpts.csv'
'duck5round01xyzpts.csv'
'duck5round02xyzpts.csv'
'duck5round06xyzpts.csv'
'duck5round17xyzpts.csv'
'duck5round18xyzpts.csv'};

refCols  = {4:6, 7:9};          % neck, tail base (priority order)
refName  = {'neck', 'tail'};

%% ---- Magnet centre ----
M = readmatrix(fullfile(magFolder, magFile));
M = M(:, 1:6);                                % pt1, pt2 = magnet corners
M = mean(M(all(~isnan(M), 2), :), 1);         % average all digitized rows
magC = (M(1:3) + M(4:6)) ./ 2;                % midpoint of the two corners
fprintf('Magnet centre: [%.4f %.4f %.4f], corner distance %.3f m\n', ...
        magC, norm(M(4:6) - M(1:3)));

% magnet frame: eW along the width (edge), eZ up, magAx = normal to face
eW = M(4:6) - M(1:3);
if isempty(magW), magW = norm(eW); end
fprintf('Magnet face from corners: %.1f deg from X (in XY)\n', atan2d(eW(2), eW(1)));
if strcmp(magFaceMode, 'wall')
    if ~isempty(wallFile)
        Wp = readmatrix(fullfile(magFolder, wallFile));
        Wp = Wp(:, wallCols);
        Wp = Wp(all(~isnan(Wp), 2), :);              % rows with both points
        dirs = Wp(:, 4:6) - Wp(:, 1:3);
        angs = atan2d(dirs(:,2), dirs(:,1));
        Wp = mean(Wp, 1);                            % average if several
        wallDir = Wp(4:6) - Wp(1:3);
        fprintf('Wall line: %d rows, angle spread between rows %.1f deg (std)\n', ...
            numel(angs), std(angs));
        fprintf('Wall line: %.1f cm long, %.1f deg from X (in XY)%s\n', ...
            100*norm(wallDir), atan2d(wallDir(2), wallDir(1)), ...
            repmat('  <-- short, angle less accurate', 1, double(norm(wallDir) < 0.3)));
    end
    if dot(wallDir, eW) < 0, wallDir = -wallDir; end
    eW = wallDir;
    eW(3) = 0;                                % face edge is horizontal
end
eW = eW ./ norm(eW);
eZ = [0 0 1] - dot([0 0 1], eW) .* eW;       % vertical, made exactly
eZ = eZ ./ norm(eZ);                          % perpendicular to the edge
magAx = cross(eW, eZ);
magAx = magAx ./ norm(magAx);
fprintf('Magnet axis (assumed): [%.3f %.3f %.3f]\n', magAx);

% magnet volume: from the top-back edge down by H, and from the back face
% toward the front by T. The front side comes from the front corner
% (frontCols); without it the thickness is centred on the back edge.
Mall = readmatrix(fullfile(magFolder, magFile));
tOff = -0.5;                                  % thickness fraction offset
if ~isempty(frontCols) && size(Mall, 2) >= frontCols(end)
    Fp = Mall(:, frontCols);
    Fp = Fp(all(~isnan(Fp), 2), :);
    if ~isempty(Fp)
        Fp = mean(Fp, 1);
        [dBack, kB] = min([norm(Fp - M(1:3)), norm(Fp - M(4:6))]);
        backC = M(3*kB-2 : 3*kB);
        frontSign = sign(dot(Fp - backC, magAx));
        magAx = frontSign .* magAx;           % axis now points to the front
        tOff = 0;                             % volume = back face -> front
        fprintf(['Front corner check: %.1f cm from the nearest back corner ' ...
                 '(magnet thickness %.1f cm)\n'], 100*dBack, 100*magT);
    end
end
magCentre = magC - (magH/2) .* eZ + magT .* (tOff + 0.5) .* magAx;
[gw, gh, gt] = ndgrid(((1:nGrid(1)) - 0.5) ./ nGrid(1) - 0.5, ...
                      -((1:nGrid(2)) - 0.5) ./ nGrid(2), ...
                      ((1:nGrid(3)) - 0.5) ./ nGrid(3) + tOff);
magSrc = magC + magW .* gw(:) .* eW + magH .* gh(:) .* eZ + ...
         magT .* gt(:) .* magAx;              % [nSrc x 3] dipole positions
fprintf('Magnet %.1f x %.1f x %.1f cm, volume centre [%.4f %.4f %.4f]\n', ...
        100*[magW magH magT], magCentre);

%% ---- Loop over trials ----
problems = cell(0, 2);          % {file, reason} – skipped trials
MagUnN = struct('file', {});
for i = 1:numel(FN)
    f = fullfile(folder, FN{i});
    if ~isfile(f)
        problems(end+1, :) = {FN{i}, 'file not found'}; %#ok<AGROW>
        fprintf('%-24s SKIPPED: file not found\n', FN{i});
        continue
    end
    try
        A = readmatrix(f);
    catch err
        problems(end+1, :) = {FN{i}, ['could not read: ' err.message]}; %#ok<AGROW>
        fprintf('%-24s SKIPPED: could not read\n', FN{i});
        continue
    end
    if size(A, 2) < 27
        problems(end+1, :) = {FN{i}, sprintf('only %d columns (need 27, Pt9 = 25:27)', size(A, 2))}; %#ok<AGROW>
        fprintf('%-24s SKIPPED: only %d columns\n', FN{i}, size(A, 2));
        continue
    end
    if ~any(all(~isnan(A(:, 25:27)), 2))
        problems(end+1, :) = {FN{i}, 'Pt9 (cols 25:27) not digitized'}; %#ok<AGROW>
        fprintf('%-24s SKIPPED: Pt9 not digitized\n', FN{i});
        continue
    end
    nF = size(A, 1);
    frame = (1:nF)';
    t = (frame - 1) ./ FR;                    % s, frame 1 = 0

    % --- Pt9 as digitized (gaps filled, smoothed) ---
    Praw = A(:, 25:27);
    [P9, a9, b9] = cleanPoint(Praw, maxGap, smoothWin);

    % --- Pt10 + Pt11 on the logger: midpoint replaces Pt9 where both exist ---
    fromDev = false(nF, 1);
    devInfo = struct('nFrames', 0, 'offset', [NaN NaN NaN], 'nOverlap', 0, ...
                     'distMean_mm', NaN, 'distSD_mm', NaN, 'seamMax_mm', NaN, ...
                     'nOnlyPt10', 0, 'nOnlyPt11', 0);
    if useDev && size(A, 2) >= devCols{2}(end)
        D10 = cleanPoint(A(:, devCols{1}), maxGap, smoothWin);
        D11 = cleanPoint(A(:, devCols{2}), maxGap, smoothWin);
        both = ~isnan(D10(:,1)) & ~isnan(D11(:,1));
        if nnz(both) >= 3
            Pd = (D10 + D11) ./ 2;
            % precision: raw distance between the 2 points (rigid logger)
            dRaw = vecnorm(A(:, devCols{1}) - A(:, devCols{2}), 2, 2);
            devInfo.distMean_mm = 1000 * mean(dRaw, 'omitnan');
            devInfo.distSD_mm   = 1000 * std(dRaw, 'omitnan');
            ov = both & ~isnan(P9(:,1));
            if nnz(ov) >= 3
                dOff = mean(P9(ov,:) - Pd(ov,:), 1);
            else
                dOff = [0 0 0];
                fprintf('%-24s NOTE: Pt10/Pt11 do not overlap Pt9 - midpoint used without offset\n', FN{i});
            end
            P9old = P9;
            P9(both,:) = Pd(both,:) + dOff;
            fromDev = both;
            % only ONE of the two visible: that point + its mean offset to the
            % (Pt9-aligned) midpoint, measured on the frames that have both
            nOne = [0 0];
            Dq = {D10, D11};
            for q = 1:2
                only = ~isnan(Dq{q}(:,1)) & ~both;
                if any(only)
                    qOff = mean(Pd(both,:) + dOff - Dq{q}(both,:), 1);
                    P9(only,:) = Dq{q}(only,:) + qOff;
                    fromDev = fromDev | only;
                    nOne(q) = nnz(only);
                end
            end
            devInfo.nOnlyPt10 = nOne(1);
            devInfo.nOnlyPt11 = nOne(2);
            % jump where the logger points meet the old Pt9 (should be ~noise)
            k  = find(fromDev(1:end-1) ~= fromDev(2:end) & ...
                      ~isnan(P9old(1:end-1,1)) & ~isnan(P9old(2:end,1)));   % seam between k and k+1
            jm = vecnorm(P9(k+1,:) - P9(k,:), 2, 2);
            devInfo.nFrames = nnz(both);
            devInfo.offset = dOff;
            devInfo.nOverlap = nnz(ov);
            if ~isempty(jm), devInfo.seamMax_mm = 1000 * max(jm); end
            a9 = find(~isnan(P9(:,1)), 1, 'first');
            b9 = find(~isnan(P9(:,1)), 1, 'last');
            stepOld = 1000 * median(vecnorm(diff(P9old), 2, 2), 'omitnan');   % normal step per frame
            if any(nOne)
                fprintf('%-24s single logger point used on %d frames (only Pt10) + %d (only Pt11)\n', FN{i}, nOne(1), nOne(2));
            end
            fprintf(['%-24s Pt10/Pt11 midpoint on %d frames | Pt10-Pt11 = %.1f +- %.1f mm (SD = digitizing precision) | ' ...
                     'offset to Pt9 %.1f mm | max jump at seam %.1f mm (normal step %.1f mm)\n'], ...
                FN{i}, devInfo.nFrames, devInfo.distMean_mm, devInfo.distSD_mm, 1000*norm(dOff), devInfo.seamMax_mm, stepOld);
        end
    end
    P = P9;
    src = zeros(nF, 1);
    src(~isnan(P(:,1))) = 1;

    % --- extend with neck / tail base ---
    offInfo = struct('ref', {}, 'side', {}, 'offset', {}, 'err', {}, ...
                     'nUsed', {}, 'frames', {});
    if useRefs && ~isempty(a9)
        for q = 1:numel(refCols)
            R = cleanPoint(A(:, refCols{q}), maxGap, smoothWin);
            both = find(~isnan(R(:,1)) & ~isnan(P9(:,1)));
            if numel(both) < 3, continue, end

            for side = {'before', 'after'}
                if strcmp(side{1}, 'before')
                    use  = both(1 : min(nOff, numel(both)));
                    gapF = (1 : a9-1)';
                else
                    use  = both(max(1, end-nOff+1) : end);
                    gapF = (b9+1 : nF)';
                end
                off = P9(use,:) - R(use,:);
                mOff = mean(off, 1);
                err = sqrt(sum((off - mOff).^2, 2));   % fit error, m

                % fill only frames still empty and where the ref exists
                fill = gapF(src(gapF) == 0 & ~isnan(R(gapF,1)));
                P(fill,:) = R(fill,:) + mOff;
                src(fill) = 1 + q;

                if ~isempty(fill)
                    offInfo(end+1) = struct('ref', refName{q}, ...
                        'side', side{1}, 'offset', mOff, ...
                        'err', [mean(err) max(err)], ...
                        'nUsed', numel(use), 'frames', fill([1 end])'); %#ok<AGROW>
                end
            end
        end
    end

    % --- jump where the extension meets the real Pt9 (should be ~1 step) ---
    seamJ = [NaN NaN];
    if a9 > 1 && ~isnan(P(a9-1, 1))
        seamJ(1) = norm(P(a9,:) - P(a9-1,:));
    end
    if b9 < nF && ~isnan(P(b9+1, 1))
        seamJ(2) = norm(P(b9+1,:) - P(b9,:));
    end
    st = sqrt(sum(diff(P9(a9:b9,:)).^2, 2));
    stepMed = median(st(~isnan(st)));

    % --- distance and expected signal ---
    D = P - magC;
    r = sqrt(sum(D.^2, 2));
    invr3 = 1 ./ r.^3;
    Bdip  = vecnorm(dipoleField(P, magCentre, magAx), 2, 2);
    Bgrid = vecnorm(dipoleField(P, magSrc, magAx), 2, 2);
    [~, iG] = max(Bgrid);

    valid = find(src > 0);
    f1 = valid(1);  f2 = valid(end);
    [rMin, iMin] = min(r);
    atEdge = (iMin - f1 < edgeFrames) || (f2 - iMin < edgeFrames);
    [~, iB] = max(Bdip);

    % --- 90 deg crossing: s(t) along the swim path ---
    % path direction: first from the whole trajectory, then refined from
    % the dirWin frames nearest the crossing (the LOCAL direction of motion
    % is what defines "perpendicular" at the closest approach)
    switch pathMode
        case 'wall',  u = eW .* sign(eW(1) + (eW(1) == 0)); nIter = 1;  % points to +X
        case 'X',     u = [1 0 0];             nIter = 1;
        case 'global',u = pathDir(D(valid, :)); nIter = 1;
        case 'local', u = pathDir(D(valid, :)); nIter = 3;
    end
    for iter = 1:nIter
        s = D * u';                            % m, along the path
        sv = s(valid);
        k = find(sign(sv(1:end-1)) ~= sign(sv(2:end)), 1);
        if ~isempty(k)                         % crossing inside the data
            c0 = max(1, min(k - floor(nFit/2), numel(valid) - nFit + 1));
            use = valid(c0 : min(c0 + nFit - 1, numel(valid)));
            crossMode = 'data';
        elseif abs(sv(1)) <= abs(sv(end))      % magnet is before the data
            use = valid(1 : min(nFit, numel(valid)));
            crossMode = 'extrap-before';
        else                                   % magnet is after the data
            use = valid(max(1, end-nFit+1) : end);
            crossMode = 'extrap-after';
        end
        pf = polyfit(frame(use), s(use), 1);
        fCross = -pf(2) / pf(1);               % fractional frame (can be < 1)
        if strcmp(pathMode, 'local')           % refine from nearby frames
            [~, ord] = sort(abs(frame(valid) - fCross));
            near = valid(sort(ord(1 : min(dirWin, numel(ord)))));
            u = pathDir(D(near, :));
        end
    end
    % --- comparison: crossing with BOTH fixed directions, every time ---
    uWall = eW .* sign(eW(1) + (eW(1) == 0));
    fCrossWall = crossFrame(D, frame, valid, uWall, nFit);
    fCrossX    = crossFrame(D, frame, valid, [1 0 0], nFit);

    extrapBy = max([0, f1 - fCross, fCross - f2]);
    % distance from the magnet at the crossing (data or extrapolated)
    Pc = interp1(frame(use), P(use,:), fCross, 'linear', 'extrap');
    rCross = norm(Pc - magC);

    % --- store ---
    MagUnN(i).file    = FN{i};
    MagUnN(i).frame   = frame;
    MagUnN(i).t       = t;
    MagUnN(i).P       = P;          % accelerometer position (incl. estimated)
    MagUnN(i).Praw    = Praw;       % Pt9 as digitized
    MagUnN(i).src     = src;        % 0 none,1 Pt9,2 neck,3 tail
    MagUnN(i).fromDev = fromDev;    % true = frame uses the Pt10/Pt11 midpoint (src is 1 there)
    MagUnN(i).devInfo = devInfo;    % Pt10/Pt11 use + precision (mm)
    MagUnN(i).offInfo = offInfo;    % offsets used + fit error (m)
    MagUnN(i).seamJump = seamJ;     % jump at [start end] seam, m (NaN = none)
    MagUnN(i).stepMed  = stepMed;   % typical Pt9 step per frame, m
    MagUnN(i).D       = D;          % vector magnet -> accelerometer
    MagUnN(i).r       = r;          % distance, m
    MagUnN(i).invr3   = invr3;      % expected signal shape, 1/m^3
    MagUnN(i).Bdip    = Bdip;       % point dipole at volume centre (a.u.)
    MagUnN(i).fBdip   = iB;         % frame of max Bdip
    MagUnN(i).Bgrid   = Bgrid;      % whole-plate model (a.u.)  <- main
    MagUnN(i).fBgrid  = iG;         % frame of max Bgrid
    MagUnN(i).s       = s;          % position along the path vs magnet, m
    MagUnN(i).pathDir = u;          % path direction used for s
    MagUnN(i).fCross  = fCross;     % 90 deg crossing frame (fractional)
    MagUnN(i).tCross  = (fCross - 1) / FR;  % s, same time base as t
    MagUnN(i).crossMode = crossMode;% 'data' / 'extrap-before' / 'extrap-after'
    MagUnN(i).extrapBy = extrapBy;  % frames outside the data (0 = inside)
    MagUnN(i).fCrossWall = fCrossWall;  % crossing along the magnet face
    MagUnN(i).fCrossX    = fCrossX;     % crossing along lab X
    MagUnN(i).rCross  = rCross;     % distance at the crossing, m
    MagUnN(i).speed   = abs(pf(1)) * FR;    % m/s along the path
    MagUnN(i).valid   = [f1 f2];    % first/last frame with any estimate
    MagUnN(i).validPt9 = [a9 b9];   % first/last frame with real Pt9
    MagUnN(i).rMin    = rMin;
    MagUnN(i).fMin    = iMin;       % frame of closest approach
    MagUnN(i).tMin    = t(iMin);
    MagUnN(i).srcMin  = src(iMin);  % was the min on real Pt9 (1) or estimated?
    MagUnN(i).atEdge  = atEdge;     % true = closest point not inside the data

    % --- report ---
    flag = '';
    if atEdge, flag = '  <-- min at edge, peak may be cut'; end
    if src(iMin) > 1, flag = [flag '  (min from ' refName{src(iMin)-1} ')']; end
    fprintf('%-24s Pt9 %3d-%3d, all %3d-%3d | r_min = %.3f m at frame %3d%s\n', ...
            FN{i}, a9, b9, f1, f2, rMin, iMin, flag);
    xflag = '';
    if extrapBy > maxExtrap, xflag = '  <-- far extrapolation, check!'; end
    fprintf('    90deg crossing at frame %6.1f (%s, %4.1f frames outside data) | r = %.3f m%s\n', ...
            fCross, crossMode, extrapBy, rCross, xflag);
    fprintf('    max: Bdip at %3d | Bgrid at %3d\n', iB, iG);
    fprintf('    crossing minus Bgrid peak:  wall %+5.1f frames | X %+5.1f frames\n', ...
            fCrossWall - iG, fCrossX - iG);
    seamTxt = {'before Pt9 starts', 'after Pt9 ends'};
    for k = 1:2
        if seamJ(k) > maxSeam * stepMed
            fprintf('    <-- jump of %.1f cm at the seam (normal step %.1f cm): check the extension %s\n', ...
                100*seamJ(k), 100*stepMed, seamTxt{k});
        end
    end
    for k = 1:numel(offInfo)
        fprintf('    %-4s %-6s frames %3d-%3d | offset err mean %.1f / max %.1f cm\n', ...
            offInfo(k).ref, offInfo(k).side, offInfo(k).frames, ...
            100*offInfo(k).err);
    end
end

MagUnN(arrayfun(@(s) isempty(s.file), MagUnN)) = [];   % drop skipped files

%% ---- Summary: which crossing definition sits closer to the Bgrid peak ----
if ~isempty(MagUnN)
    dW = [MagUnN.fCrossWall] - [MagUnN.fBgrid];
    dX = [MagUnN.fCrossX]    - [MagUnN.fBgrid];
    fprintf('\nCrossing minus Bgrid peak, all %d trials (frames):\n', numel(MagUnN));
    fprintf('   wall : mean %+5.1f | mean |diff| %4.1f | max |diff| %4.1f\n', ...
            mean(dW), mean(abs(dW)), max(abs(dW)));
    fprintf('   X    : mean %+5.1f | mean |diff| %4.1f | max |diff| %4.1f\n', ...
            mean(dX), mean(abs(dX)), max(abs(dX)));
    fprintf('   (smaller |diff| = that definition follows the magnetic peak better)\n');
end

%% ---- Report skipped files (screen + text file to fix later) ----
probFile = fullfile(outFolder, 'MagDistUnN_problems.txt');
if isempty(problems)
    fprintf('\nAll %d files processed.\n', numel(FN));
    if isfile(probFile), delete(probFile); end      % no stale list
else
    fprintf('\n%d of %d files SKIPPED (list saved to %s):\n', ...
            size(problems, 1), numel(FN), probFile);
    fid = fopen(probFile, 'w', 'n', 'UTF-8');
    fprintf(fid, 'MagDistUnN – skipped files, %s\n', datestr(now));
    fprintf(fid, 'Folder: %s\n\n', folder);
    for k = 1:size(problems, 1)
        fprintf('   %-24s %s\n', problems{k, :});
        fprintf(fid, '%-24s %s\n', problems{k, :});
    end
    fclose(fid);
end

%% ---- Plots: r(t) and 1/r^3, colored by source ----
if doPlot
    cols = [0 0 0; 0.85 0.45 0; 0 0.45 0.75];   % Pt9, neck, tail
    for i = 1:numel(MagUnN)
        S = MagUnN(i);
        % A: distance | B: expected signal | C: position along X
        %   line color = where the accelerometer position came from
        %   (black = Pt9 digitized, orange = from neck, blue = from tail)
        %   gray dotted lines = first/last frame of real Pt9
        %   magenta dashed line = 90 deg crossing (passing the magnet's X)
        figure('Name', S.file, 'Color', 'w');
        xl = [min(S.valid(1), floor(S.fCross)) - 5, ...
              max(S.valid(2), ceil(S.fCross)) + 5];
        for p = 1:3
            subplot(3,1,p); hold on
            switch p
                case 1, y = S.r;
                case 2, y = S.invr3 ./ max(S.invr3);
                case 3, y = S.s;
            end
            hs = gobjects(1, 3);
            for q = 1:3
                yq = y; yq(S.src ~= q) = NaN;
                hs(q) = plot(S.frame, yq, '-', 'Color', cols(q,:), 'LineWidth', 1.6);
            end
            xline(S.validPt9, ':', 'Color', [0.5 0.5 0.5]);
            xline(S.fCross, 'm--', 'LineWidth', 1.2);
            xlim(xl); box on
            switch p
                case 1
                    plot(S.fMin, S.r(S.fMin), 'ro', 'MarkerFaceColor', 'r');
                    ylabel('r (m)');
                    title({strrep(S.file, '_', '\_'), ...
                        'A. Distance accelerometer - magnet  (red dot = minimum)'});
                    legend(hs, {'Pt9 digitized', 'from neck', 'from tail'}, ...
                        'Location', 'best');
                case 2
                    h1 = plot(S.frame, S.Bdip ./ max(S.Bdip), ':', ...
                        'Color', [0.3 0.3 0.3], 'LineWidth', 1.4);
                    h2 = plot(S.frame, S.Bgrid ./ max(S.Bgrid), '-', ...
                        'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
                    ylabel('normalized (max = 1)');
                    title('B. Expected magnetometer signal');
                    legend([hs(1) h1 h2], {'1/r^3 (color = source)', ...
                        'B point dipole', 'B whole plate'}, 'Location', 'best');
                case 3
                    yline(0, 'k:');
                    plot(S.fCross, 0, 'mo', 'MarkerFaceColor', 'm');
                    ylabel('s (m)'); xlabel('frame (250 fps)');
                    title(sprintf('C. Position along the path vs magnet (pathMode ''%s'', 0 = passing the magnet)', pathMode));
            end
        end
    end
end

save(fullfile(outFolder, 'MagDistUnN_results.mat'), 'MagUnN', 'magC', ...
     'magAx', 'magCentre', 'magSrc', 'problems');
end

%% ---- helper: field of equal dipoles (moment along ax), summed ----
% P   : [nF x 3] field points,  src : [nS x 3] dipole positions
% B   : [nF x 3], arbitrary units (only the shape matters)
function B = dipoleField(P, src, ax)
B = zeros(size(P));
for k = 1:size(src, 1)
    d = P - src(k,:);
    r = sqrt(sum(d.^2, 2));
    md = (d * ax') ./ r;                       % m . d_hat
    B = B + (3 .* md .* d ./ r - ax) ./ r.^3;
end
B = B ./ size(src, 1);
end

%% ---- helper: 90 deg crossing frame for a fixed direction u ----
% same method as in the main loop: s = position along u, linear fit of
% nFit frames around (or nearest to) s = 0
function fCross = crossFrame(D, frame, valid, u, nFit)
s = D * u';
sv = s(valid);
k = find(sign(sv(1:end-1)) ~= sign(sv(2:end)), 1);
if ~isempty(k)
    c0 = max(1, min(k - floor(nFit/2), numel(valid) - nFit + 1));
    use = valid(c0 : min(c0 + nFit - 1, numel(valid)));
elseif abs(sv(1)) <= abs(sv(end))
    use = valid(1 : min(nFit, numel(valid)));
else
    use = valid(max(1, end-nFit+1) : end);
end
pf = polyfit(frame(use), s(use), 1);
fCross = -pf(2) / pf(1);
end

%% ---- helper: main direction of a set of points (PCA), pointing to +X ----
function u = pathDir(X)
[~, ~, V] = svd(X - mean(X, 1), 'econ');
u = V(:, 1)';
if u(1) < 0, u = -u; end
end

%% ---- helper: fill short gaps + smooth, inside the digitized range ----
function [Pc, a, b] = cleanPoint(P, maxGap, win)
Pc = P;
ok = ~isnan(P(:,1));
a = find(ok, 1, 'first');
b = find(ok, 1, 'last');
if isempty(a), return, end
Pc(a:b,:) = fillmissing(P(a:b,:), 'linear', 'MaxGap', maxGap);
if win > 0
    gap = isnan(Pc(:,1));                        % long gaps stay NaN
    Pc(a:b,:) = smoothdata(Pc(a:b,:), 'gaussian', win, 'omitnan');
    Pc(gap,:) = NaN;
end
end