%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of truss                       *
%*                                                                 *
%* - Call procedures:                                              *
%*     TrussStiff.m - InitTruss()                                  *
%*     ./ReadTruss.m - ReadTruss()                                 *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function AxisymStiff()

% Init variables of the element
InitStrain();

% Read Material and Elements
AxisymReadStrain();

fprintf('Solution phase ...\n\n');

<<<<<<< HEAD
% calculate addresses of diagonal elements
Addres();
=======
>>>>>>> 增添了三种质量矩阵生成方式

% Data check Or Solve
global cdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
AxisymAssemble();




end

% ----------------------- Functions -----------------------------------

% Init parameters of truss element
function InitStrain()
global sdata;
sdata.NNODE = 8;  %���ｫ4�ĳ���8
sdata.NDOF = 2;

end



