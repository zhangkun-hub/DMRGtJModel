classdef ModernPalettes
    % ModernPalettes provides curated color schemes for scientific visualization
    % Designed for high contrast, accessibility, and modern aesthetics.
    
    properties (Constant)
        % === Primary Themes ===
        
        % Theme 1: Amber & Cyan (Current Hopping)
        % Warm, energetic vs Cool, tech-like. Great for kinetic/dynamic terms.
        AmberCyan = struct('Positive', [255, 171, 0]/255, ... % #FFAB00
                           'Negative', [0, 188, 212]/255);    % #00BCD4
                           
        % Theme 2: Teal & Orange (Current SC)
        % Classic complementary pair. Natural, balanced. Good for correlations.
        TealOrange = struct('Positive', [0, 150, 136]/255, ... % #009688 (Teal)
                            'Negative', [255, 112, 67]/255);   % #FF7043 (Deep Orange)

        % Theme 3: Indigo & Pink (Vibrant / Modern Science)
        % Very distinct, high energy. Good for distinguishing subtle phases.
        IndigoPink = struct('Positive', [63, 81, 181]/255, ...  % #3F51B5 (Indigo)
                            'Negative', [233, 30, 99]/255);     % #E91E63 (Pink)

        % Theme 4: Emerald & Violet (Soft / Nature)
        % Calming, distinct. Good for density or background-heavy plots.
        EmeraldViolet = struct('Positive', [76, 175, 80]/255, ... % #4CAF50 (Emerald)
                               'Negative', [156, 39, 176]/255);   % #9C27B0 (Violet)
        
        % Theme 5: Blue & Gold (Classic Academic / Colorblind Safe)
        % High visibility, very safe for print.
        BlueGold = struct('Positive', [255, 193, 7]/255, ...    % #FFC107 (Amber/Gold)
                          'Negative', [33, 150, 243]/255);      % #2196F3 (Blue)

        % Theme 6: Deep Purple & Lime (High Contrast / Dark Mode friendly)
        PurpleLime = struct('Positive', [205, 220, 57]/255, ...  % #CDDC39 (Lime)
                            'Negative', [103, 58, 183]/255);     % #673AB7 (Deep Purple)
    end
    
    methods (Static)
        function colors = Get(themeName)
            % Returns the color struct for a given theme name
            % Usage: c = ModernPalettes.Get('TealOrange');
            try
                colors = ModernPalettes.(themeName);
            catch
                error('ModernPalettes:UnknownTheme', ...
                    'Theme "%s" not found. Options: AmberCyan, TealOrange, IndigoPink, EmeraldViolet, BlueGold, PurpleLime', themeName);
            end
        end
        
        function Demo()
            % Displays a figure showing all available palettes
            p = properties('ModernPalettes');
            figure('Name', 'Modern Color Palettes', 'Color', 'w', 'Position', [100, 100, 600, 400]);
            axis off; hold on;
            
            y = 1;
            h = 0.8;
            w = 2;
            
            for i = 1:length(p)
                theme = ModernPalettes.(p{i});
                
                % Name
                text(0, y, p{i}, 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
                
                % Positive
                rectangle('Position', [2.5, y-h/2, w, h], 'FaceColor', theme.Positive, 'EdgeColor', 'none');
                text(2.5+w/2, y, 'Pos', 'Color', 'w', 'HorizontalAlignment', 'center');
                
                % Negative
                rectangle('Position', [5, y-h/2, w, h], 'FaceColor', theme.Negative, 'EdgeColor', 'none');
                text(5+w/2, y, 'Neg', 'Color', 'w', 'HorizontalAlignment', 'center');
                
                y = y - 1.5;
            end
            xlim([-1, 8]);
            ylim([y, 2]);
            title('Modern Scientific Color Palettes');
        end
    end
end

