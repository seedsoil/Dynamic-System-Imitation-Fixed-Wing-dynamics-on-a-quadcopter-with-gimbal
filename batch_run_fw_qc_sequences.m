%% batch_run_fw_qc_sequences.m
% Batch driver for the functionized compare:
%   S = run_fw_qc_compare(seq, Name=Value, ...)
%
% Features
%   • Option A: explicit seq list (default below)
%   • Option B: auto-generate N random seq
%   • Serial or parallel (parfor) execution
%   • Forces exports under <projectRoot>/Plots/
%   • Writes a summary CSV at the end

%% --------------------- SETUP ---------------------
funcName = 'run_fw_qc_compare.m';
fp = which(funcName);
if isempty(fp)
    error('Function file not found on path: %s', funcName);
end
projectRoot  = fileparts(fp);
basePlotsDir = fullfile(projectRoot, 'Plots');      % force common output root
if ~isfolder(basePlotsDir), mkdir(basePlotsDir); end   % ok if exists (mkdir/isfolder)  % (docs) :contentReference[oaicite:0]{index=0}

% Execution mode
useParallel = false;    % set true to use a parallel pool (auto-starts if available) :contentReference[oaicite:1]{index=1}

% Common options passed to run_fw_qc_compare
commonOpts = struct( ...
    'PlotsRoot',        string(basePlotsDir), ...
    'ExportVariants',   ["both","on"], ...
    'SavePNG',          true, ...
    'SavePDF',          false, ...
    'Verbose',          true, ...
    'MarkImpulses',     false);

%% --------------------- SEQUENCES ---------------------
% ===== Option A: explicit list of sequences =====
seqList = {
    [1 +0; 1 +0; 1 +0; 1 +0; 1 +0; 1 +0; 1 +0; 1 +0]                                     % seq_1p0_1p0_1p0_1p0_1p0_1p0_1p0_1p0
    [2 +5; 1 +0; 1 +0; 1 +0; 2 -5; 1 +0; 1 +0; 1 +0]                                     % seq_2p5_1p0_1p0_1p0_2n5_1p0_1p0_1p0
    [3 +0; 1 +0; 3 +5; 1 +0; 3 +0; 1 +0; 3 -5; 1 +0]                                     % seq_3p0_1p0_3p5_1p0_3p0_1p0_3n5_1p0
    [3 +0; 1 +0; 4 +5; 1 +0; 3 +0; 1 +0; 4 -5; 1 +0]                                     % seq_3p0_1p0_4p5_1p0_3p0_1p0_4n5_1p0
    [3 +2; 2 +5; 4 -4; 3 +4; 2 -2; 4 +2; 4 -5; 2 -1]                                     % seq_3p2_2p5_4n4_3p4_2n2_4p2_4n5_2n1
};

% ===== Option B: ONLY used if seqList is empty/missing =====
if ~exist('seqList','var') || isempty(seqList)
    numRuns      = 1;                      % total random runs
    numPulses    = 8;                      % rows per seq
    channelSet   = [1 2 3 4];              % 1=THR 2=ELEV 3=AIL 4=RUD
    scaleChoices = [2 3 4 5 6];
    rng(1);
    seqList = cell(numRuns,1);
    for k = 1:numRuns
        ch  = channelSet(randi(numel(channelSet),numPulses,1));
        mag = scaleChoices(randi(numel(scaleChoices),numPulses,1));
        sgn = sign(randn(numPulses,1)); sgn(sgn==0)=1;
        seqList{k} = [ch, mag.*sgn];
    end
end

%% --------------------- RUN ALL ---------------------
nRuns  = numel(seqList);
status = strings(nRuns,1);
elapsed= zeros(nRuns,1);
tok    = strings(nRuns,1);
metrics= cell(nRuns,1);

tAll = tic;
fprintf('\n[Batch] Starting %d runs (%s)…\n', nRuns, char(tern(useParallel,'parallel','serial')));

if useParallel
    % A pool starts automatically when needed (Parallel Computing Toolbox). :contentReference[oaicite:2]{index=2}
    parfor i = 1:nRuns
        seq = seqList{i}; %#ok<PFBNS>
        tRun = tic;
        try
            Si = run_fw_qc_compare(seq, ...
                 PlotsRoot=commonOpts.PlotsRoot, ExportVariants=commonOpts.ExportVariants, ...
                 SavePNG=commonOpts.SavePNG, SavePDF=commonOpts.SavePDF, ...
                 Verbose=commonOpts.Verbose, MarkImpulses=commonOpts.MarkImpulses);  % name=value syntax :contentReference[oaicite:3]{index=3}
            status(i)  = "ok";
            elapsed(i) = toc(tRun);
            tok(i)     = string(Si.seqToken);
            metrics{i} = Si.metrics;
        catch ME
            status(i)  = "fail";
            elapsed(i) = toc(tRun);
            tok(i)     = string(seq_to_token(seq));
            metrics{i} = struct();
            fprintf(2,'[Batch] #%d %s FAILED: %s\n', i, char(tok(i)), ME.message);
        end
    end
else
    for i = 1:nRuns
        seq = seqList{i};
        tRun = tic;
        tki  = seq_to_token(seq);
        fprintf('[Batch] Run %d/%d: %s\n', i, nRuns, char(tki));
        try
            Si = run_fw_qc_compare(seq, ...
                 PlotsRoot=commonOpts.PlotsRoot, ExportVariants=commonOpts.ExportVariants, ...
                 SavePNG=commonOpts.SavePNG, SavePDF=commonOpts.SavePDF, ...
                 Verbose=commonOpts.Verbose, MarkImpulses=commonOpts.MarkImpulses);  % name=value syntax :contentReference[oaicite:4]{index=4}
            status(i)  = "ok";
            elapsed(i) = toc(tRun);
            tok(i)     = string(Si.seqToken);
            metrics{i} = Si.metrics;
        catch ME
            status(i)  = "fail";
            elapsed(i) = toc(tRun);
            tok(i)     = string(tki);
            metrics{i} = struct();
            fprintf(2,'[Batch] #%d %s FAILED: %s\n', i, char(tok(i)), ME.message);
        end
    end
end

Tsec = toc(tAll);
fprintf('[Batch] Done in %.1f s. Success: %d / %d\n', Tsec, sum(status=="ok"), nRuns);

%% --------------------- SUMMARY CSV ---------------------
% Build a flat table and save under <projectRoot>/Plots/summary_*.csv
when  = datetime('now');                               % timestamp (datetime) :contentReference[oaicite:5]{index=5}
summ  = table((1:nRuns).', tok, status, elapsed(:), ...
              'VariableNames', {'Run','SeqToken','Status','Elapsed_s'});

% include a few key metrics if present
addField = @(name,getter)  add_metric_col(summ, name, getter(metrics));
summ = addField('MAE3_off', @(M) getfield_safe(M,'MAE3_off'));
summ = addField('MAE3_on',  @(M) getfield_safe(M,'MAE3_on'));
summ = addField('POV_MAE_on', @(M) getfield_safe(M,'MAE_pov_on'));
summ = addField('Imp_Pos_%',  @(M) getfield_safe(M,'improvements','pos'));

sumName = fullfile(basePlotsDir, sprintf('summary_%s.csv', datestr(when,'yyyymmdd_HHMMSS')));
writetable(summ, sumName);                               % write CSV (writetable) :contentReference[oaicite:6]{index=6}
fprintf('[Batch] Wrote summary: %s\n', sumName);

%% --------------------- HELPERS ---------------------
function out = add_metric_col(tbl, name, vals)
    v = nan(height(tbl),1);
    for ii=1:numel(vals), if ~isempty(vals{ii}), v(ii) = vals{ii}; end, end
    tbl.(name) = v; out = tbl;
end

function vals = getfield_safe(M, field, subfield)
    n = numel(M); vals = cell(n,1);
    for ii=1:n
        if isempty(M{ii}), vals{ii} = []; continue; end
        if nargin==2
            if isfield(M{ii}, field), vals{ii} = M{ii}.(field); else, vals{ii} = []; end
        else
            if isfield(M{ii}, field) && isfield(M{ii}.(field), subfield)
                vals{ii} = M{ii}.(field).(subfield);
            else
                vals{ii} = [];
            end
        end
    end
end

function tok = seq_to_token(seq)
    n = size(seq,1); s = strings(n,1);
    for i=1:n
        ch = seq(i,1); k = seq(i,2);
        signChar = "p"; if k < 0, signChar = "n"; end
        s(i) = string(ch) + signChar + string(abs(k));
    end
    tok = "seq_" + strjoin(s,"_");
end

function out = tern(tf,a,b), if tf, out=a; else, out=b; end
end