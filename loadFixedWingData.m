% Load the fixed-wing data from Excel (already in quaternion form)
filename = 'FW20SecondsRun.xlsx';

try
    % Read the data (adjust sheet if needed)
    data = readtable(filename, 'Sheet', 1);

    % Expected column names in the file:
    % time, outputStates1..13 (=[x y z u v w q0 q1 q2 q3 p q r]), TEAROut1..4 (= [T E A R])
    colNames = [{'time'}, ...
                arrayfun(@(k)sprintf('outputStates%d',k),1:13,'uni',false), ...
                arrayfun(@(k)sprintf('TEAROut%d',k),1:4,'uni',false)];
    n = min(numel(colNames), width(data));
    data.Properties.VariableNames(1:n) = colNames(1:n);

    % Optional sanity check: quaternion near unit norm (no normalization performed)
    S = data{:, 2:14}; % outputStates1..13 in assumed order
    qn = sqrt(S(:,7).^2 + S(:,8).^2 + S(:,9).^2 + S(:,10).^2); % q0..q3 are cols 7..10
    if any(abs(qn - 1) > 1e-3)
        warning('FW trace quaternions deviate from unit norm by >1e-3 at %d samples.', sum(abs(qn-1) > 1e-3));
    end

    % Build FW_States in quaternion layout expected by your sim:
    % [time, x y z, u v w, q0 q1 q2 q3, p q r]
    FW_States = [ data.time, S ];

    % Build FW_Controls with time as first column (TEAROut = [Thrust Elevator Aileron Rudder])
    FW_Controls = [ data.time, data.TEAROut1, data.TEAROut2, data.TEAROut3, data.TEAROut4 ];

    % Export to base workspace
    assignin('base','FW_States',   FW_States);
    assignin('base','FW_Controls', FW_Controls);

    fprintf('Imported FW_States (quat) [%d×14] and FW_Controls [%d×5] to base workspace.\n', ...
            size(FW_States,1), size(FW_Controls,1));

catch ME
    warning('Unable to load "%s". Details: %s', filename, ME.message);
end
