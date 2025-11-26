function [positive_sc_color, negative_sc_color] = normalized_sc_corr_colors()
%NORMALIZED_SC_CORR_COLORS Provide the shared palette for normalized SC corr plots.
%   Returns two RGB triplets (0-1 range) tailored for publication figures.
%   The positive channel uses a calm blue (#1F77B4) and the negative channel
%   uses a balanced crimson (#D62728). Both colors offer strong contrast
%   while remaining colorblind-friendly.
% light purple and light red
positive_sc_color = [142   139        254]/256;
negative_sc_color = [232   132        130]/256;

% light orange and blue
% positive_sc_color = [250        168     53      ]/256;
% negative_sc_color = [104        204     217             ]/256;

% positive_sc_color = [73   108        136]/256;
% negative_sc_color = [254   178        180]/256;
end

