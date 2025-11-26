function cfg = GetCommonConfig()
    % GETCOMMONCONFIG Returns the centralized configuration for simulation parameters.
    % Edit this file to change parameters for ALL plotting scripts.
    
    % ====== System Parameters ======
    cfg.Ly = 6; 
    cfg.Lx = 16;
    cfg.t2 = -0.4;
    cfg.J = 0.4;
    
    % ====== Filling ======
    % Specify 'Hole' or 'Mu'
    cfg.FillingType = 'Hole'; 
    cfg.FillingValue = 6;
    % Example for Mu:
    % cfg.FillingType = 'Mu'; 
    % cfg.FillingValue = 2.5;
    
    % ====== Simulation ======
    cfg.Geometry = 'OBC';   % 'OBC' or 'Cylinder'
    cfg.refer_site = 27;    % Reference site index (0-based)
    cfg.pairing_type = 's'; % 's' for singlet, 't' for triplet
    
    % ====== Bond Dimensions ======
    % Provide a list. 
    % Scripts supporting multi-line plots will show all.
    % Scripts plotting spatial maps will loop and generate one figure per D.
    cfg.D_values = [20000]; 
    
end
