function filename_postfix = GetSimulationParamsPostfix(Ly, Lx, t2, J, FillingType, FillingValue, D, Geometry)
% GETSIMULATIONPARAMSPOSTFIX Generates the standard filename suffix for simulation data.
%
% Usage:
%   postfix = GetSimulationParamsPostfix(Ly, Lx, t2, J, FillingType, FillingValue, D, Geometry)
%
% Arguments:
%   Ly           - System height (integer)
%   Lx           - System width (integer)
%   t2           - Next-nearest neighbor hopping (double)
%   J            - Exchange interaction (double)
%   FillingType  - String: 'Hole' or 'Mu'
%   FillingValue - Value for the filling parameter (double/int)
%   D            - Bond dimension (integer)
%   Geometry     - String: 'OBC', 'Cylinder', etc.
%
% Returns:
%   filename_postfix - A string string like '6x16t2=0J=0.4Hole6D20000OBC.json'

    % Validate inputs lightly (Linus would say: don't crash if I pass a float to int, just work)
    
    % Construct the filling part string
    if strcmpi(FillingType, 'Hole')
        filling_part = ['Hole', num2str(FillingValue)];
    elseif strcmpi(FillingType, 'Mu')
        filling_part = ['Mu', num2str(FillingValue)];
    else
        error('GetSimulationParamsPostfix:InvalidFillingType', ...
              'FillingType must be "Hole" or "Mu". Got: %s', FillingType);
    end

    % Construct the full string
    % Pattern: {Ly}x{Lx}t2={t2}J={J}{FillingPart}D{D}{Geometry}.json
    
    filename_postfix = [ ...
        num2str(Ly), 'x', num2str(Lx), ...
        't2=', num2str(t2), ...
        'J=', num2str(J), ...
        filling_part, ...
        'D', num2str(D), ...
        Geometry, '.json'];

end

