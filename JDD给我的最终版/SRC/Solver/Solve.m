%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To solve finite element static equilibrium equations        *
%*                                                                 *
%* - Call procedures:                                              *
%*     ./LDLTFactor.m            - LDLTFactor()                    *
%*     Solve.m                   - Stiff2Sparse()                  *
%*     ./ColSol.m                - ColSol()                        *  
%*     Solve.m                   - WriteDis()                      *
%*     SRC/Mechanics/GetStress.m - GetStress()                     *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function Solve()

global cdata;
global sdata;

NEQ = sdata.NEQ;
NLCASE = cdata.NLCASE;
MODEX = cdata.MODEX;
sdata.DIS = zeros(NEQ, NLCASE, 'double');
sdata.V = zeros(NEQ, NLCASE, 'double');
sdata.STRAIN = zeros(NEQ, NLCASE, 'double');
sdata.STRESS = zeros(NEQ, NLCASE, 'double');

% The pre-process of Solution
% MODEX = 1, œ° Ë«ÛΩ‚     
% MODEX = 2, nnewmark


for L = 1:NLCASE
    if (MODEX == 1) 
    
        sdata.DIS(:,L) = sdata.STIFF \ sdata.R(:,L);
    %   Print displacements
        WriteDis(L);
    
    %   Calculation of stresses
        GetStress(L);
        
        vtkwrite('PARA',0);
    elseif (MODEX == 2)
       
%         beta = Opt_beta();
   beta = 0;     
        nnewmark(L,beta);
        
    elseif (MODEX == 3)
        
%         beta = Opt_beta();
beta = 0;
        newmark(L,beta);
        
    end
    
end


end

% ----------------------- Functions -----------------------------------


% Print Displacements
function WriteDis(NUM)

% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
NUMNP = cdata.NUMNP;
DIS = sdata.DIS(:, NUM); ID = sdata.ID;

fprintf(IOUT, '\n\n LOAD CASE %3d', NUM);
fprintf(IOUT, ['\n\n D I S P L A C E M E N T S\n' ...
    '\n       NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n']);

D = zeros(3, 1, 'double');
for II = 1:NUMNP
    D(:) = 0;
    if (ID(1, II) ~= 0) D(1) = DIS(ID(1, II)); end
    if (ID(2, II) ~= 0) D(2) = DIS(ID(2, II)); end
    if (ID(3, II) ~= 0) D(3) = DIS(ID(3, II)); end
    
    fprintf(IOUT, ' %10d        %18.6e%18.6e%18.6e\n', II, D(1), D(2), D(3));
end

end