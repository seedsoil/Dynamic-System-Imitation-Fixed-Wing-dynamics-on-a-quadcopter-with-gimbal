%% run_60s_seq_noThrust.m
% Example: 60 s sim, no thrust channel, max ±5 scaling

% Make sure run_fw_qc_compare.m is on your path
% (same folder, or use addpath)

% % --- Build sequence ---
% % Channels: 2=elev, 3=ail, 4=rud
% % Magnitudes: integers in [-5,5], avoiding 0
% rng default  % reproducibility
% numPulses = 12;
% channels  = randi([2 4], numPulses, 1);
% magnitudes = randi([-15 15], numPulses, 1);
% 
% % eliminate any zeros (replace with ±1)
% magnitudes(magnitudes==0) = sign(randn(sum(magnitudes==0),1));
% 
% seq = [channels magnitudes];

% --- Build sequence(s): long ones only ---
% Channels: 2=elev, 3=ail, 4=rud
seqList = {
    [4 +2; 4 +1; 2 +2; 4 -2; 3 +1; 2 +2; 2 +1; 3 +2; 4 +1; 4 -2; 2 +2; 4 +2]    % seq_4p2_4p1_2p2_4n2_3p1_2p2_2p1_3p2_4p1_4n2_2p2_4p2
    [4 +5; 4 +1; 2 +3; 4 -4; 3 -1; 2 +5; 2 +3; 3 +5; 4 +2; 4 -5; 2 +4; 4 +5]    % seq_4p5_4p1_2p3_4n4_3n1_2p5_2p3_3p5_4p2_4n5_2p4_4p5
    %[4 +15; 4 -5; 2 +12; 4 -8; 3 -10; 2 +9; 2 +7; 3 +12; 4 +5; 4 -11; 2 +12; 4 +9]
};

% --- Run all sequences ---
for i = 1:numel(seqList)
    seq = seqList{i};
    S = run_fw_qc_compare(seq, ...
        TsimTotal=60, ...
        TAnalyzeTo=60, ...
        AutoFitSchedule=true, ...
        SettleMargin=2, ...
        Verbose=true);

    % --- Display summary ---
    disp("Seq " + string(i) + " results written to: " + S.outDir);
    disp("Evaluation metrics:");
    disp(S.metrics);
end


% --- Run comparison ---
S = run_fw_qc_compare(seq, ...
    TsimTotal=60, ...
    TAnalyzeTo=60, ...
    AutoFitSchedule=true, ...
    SettleMargin=2, ...
    Verbose=true);

% --- Display summary ---
disp("Results written to: " + S.outDir);
disp("Evaluation metrics:");
disp(S.metrics);
