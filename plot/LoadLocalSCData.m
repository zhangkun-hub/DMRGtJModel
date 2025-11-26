function [sites_set, distance, sc_correlation] = LoadLocalSCData(FileNamePrefix, FileNamePostFix, Ly, pairing_spin_sym)
% FileNamePreFix is somthing like path + "scsyy" or "scsyx"
% FileNamePostfix contain the ".json" postfix
% return:
% sites_set: 4 by N matrix. N is # of data
% distance: N column matrix.
% sc_correlation: N column matrix.
% pairing_spin_sym = 's' or 't', represent singlet or triplet pairing.
% 't' only give the 0-th component of triplet pairing correlation

A = jsondecode(fileread([FileNamePrefix,'a',FileNamePostFix]));
B = jsondecode(fileread([FileNamePrefix,'b',FileNamePostFix]));
C = jsondecode(fileread([FileNamePrefix,'c',FileNamePostFix]));
D = jsondecode(fileread([FileNamePrefix,'d',FileNamePostFix]));

sites_set=zeros(numel(A), 4);
distance=zeros(numel(A), 1);
for i=1:numel(A)
    for j = 1:4
        sites_set(i, j) = A{i}{1}(j);
    end
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/Ly;
end



sign_of_spin_pairing_symmetry = 0;
if pairing_spin_sym == 's'
    sign_of_spin_pairing_symmetry = -1;
    sc_correlation = zeros(numel(A), 1);
elseif pairing_spin_sym == 't'
    sign_of_spin_pairing_symmetry = 1;
    sc_correlation = zeros(numel(A), 3);
else
    error('pairing_spin_sym must be either s or t');
end


if numel(A{1}{2}) == 1 %real number DMRG code
    fprintf("Use real number DMRG data!\n");
    for i=1:numel(A)
        %tJ like singlet: A-B
        if( pairing_spin_sym == 's')
            sc_correlation(i) = A{i}{2} ...
                + sign_of_spin_pairing_symmetry *B{i}{2};
        else
            sc_correlation(i,1) = (A{i}{2} ...
                + sign_of_spin_pairing_symmetry *B{i}{2})/sqrt(2);
            sc_correlation(i,2) = C{i}{2} ;
            sc_correlation(i,3) = D{i}{2} ;
        end
        %tJ like triplet: A+B, C, D


        % sc_correlation(i) = D{i}{2};
        %for hubbard like singlet pairing: sc_correlation(i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
else%complex number DMRG
    fprintf("Use complex number DMRG data!\n");
    for i=1:numel(A)
        %tJ like
        if( pairing_spin_sym == 's')
            sc_correlation(i) = A{i}{2}(1) + (1i)* A{i}{2}(2) ...
                + sign_of_spin_pairing_symmetry * (B{i}{2}(1) + (1i)* B{i}{2}(2));
        else
            sc_correlation(i,1) = A{i}{2}(1) + (1i)* A{i}{2}(2) ...
                + sign_of_spin_pairing_symmetry * (B{i}{2}(1) + (1i)* B{i}{2}(2));
            sc_correlation(i,2) = C{i}{2}(1);
            sc_correlation(i,3) = D{i}{2}(1) ;
        end
    end
end

end