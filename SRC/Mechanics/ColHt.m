%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate column heights                                 *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Truss/ReadTruss.m                                         *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     LeiYang Zhao, Yan Liu, Computational Dynamics Group,        *
%*     School of Aerospace Engineering, Tsinghua University,       *
%*     2019.02.22                                                  *
%*                                                                 *
%* *****************************************************************

function ColHt(LM) %每个单元的LM，这仅仅是一列

% Get global data
global sdata;
MHT = sdata.MHT;%MHT中一共有K阶元素
LS = min(LM(LM ~= 0)); %LM中最小的元素
ND = sdata.NDOF * sdata.NNODE; %单元自由度，LM的长度
for I = 1:ND
    II = LM(I);
    if (II ~= 0)
        ME = II - LS;
        if (ME > MHT(II)) MHT(II) = ME; end
    end
end

sdata.MHT = MHT;

end