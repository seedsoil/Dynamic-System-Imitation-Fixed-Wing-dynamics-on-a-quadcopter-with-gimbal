% token_to_seq.m
% Convert a token string into seq vector and print it in one line with semicolons.

% === define your token here ===
tok = "seq_4p2_3p1_2p2_4n2_3p1_2p2_4p1_2p1_3p2_4n1_2p2_3n1_4p2_2p2_3p1_4p1_3n2_2p2_4n1_2p1_3p2_4p1_2p2_3p1";

% --- Parse token ---
parts = split(extractAfter(tok,"seq_"), "_");
n = numel(parts);
seq = zeros(n,2);

for i = 1:n
    p = char(parts(i));
    m = regexp(p, '^(\d+)([pnPN])(\d+)$', 'tokens', 'once');
    if isempty(m)
        error("Bad token part: %s", p);
    end
    ch  = str2double(m{1});
    sgn = lower(m{2});
    if sgn == 'p'
        val = +str2double(m{3});
    else
        val = -str2double(m{3});
    end
    seq(i,:) = [ch, val];
end

% --- Print nicely in one line ---
fprintf('seq = [%s];\n', sprintf('%d %d;', seq.'));
