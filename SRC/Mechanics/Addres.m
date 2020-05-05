%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate addresses of diagonal elements in banded       *
%*     matrix whose column heights are known                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Truss/TrussStiff.m                                        *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     LeiYang Zhao, Yan Liu, Computational Dynamics Group,        *
%*     School of Aerospace Engineering, Tsinghua University,       *
%*     2019.02.22                                                  *
%*                                                                 *
%* *****************************************************************

function Addres()

% Get global data
global sdata;
global cdata;

NEQ = sdata.NEQ; MHT = sdata.MHT;
sdata.MAXA = zeros(NEQ+1, 1, 'int64');  %对角线元素的地址
MAXA = sdata.MAXA;

MAXA(1) = 1;
MAXA(2) = 2;
MK = 0;

if (NEQ > 1)
    for I = 2:NEQ
        if (MHT(I) > MK) MK = MHT(I); end
        MAXA(I+1) = MAXA(I) + MHT(I) + 1; %通过列高来计算对角线元素的地址
    end
end

sdata.MK = MK + 1;  %最大半带宽
sdata.NWK = MAXA(NEQ+1) - MAXA(1);  %K中skyline一共多少个元素
sdata.MAXA = MAXA;

% Write total system data
MM = round(sdata.NWK / NEQ);
fprintf(cdata.IOUT, ['\n\n  TOTAL SYSTEM DATA\n\n' ...
    '     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = %10d\n' ...
    '     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = %10d\n' ...
    '     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = %10d\n' ...
    '     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = %10d\n'], ...
    NEQ, sdata.NWK, sdata.MK, MM);

end