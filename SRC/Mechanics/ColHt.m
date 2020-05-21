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

function ColHt(LM) %ÿ����Ԫ��LM���������һ��

% Get global data
global sdata;
MHT = sdata.MHT;%MHT��һ����K��Ԫ��
LS = min(LM(LM ~= 0)); %LM����С��Ԫ��
ND = sdata.NDOF * sdata.NNODE; %��Ԫ���ɶȣ�LM�ĳ���
for I = 1:ND
    II = LM(I);
    if (II ~= 0)
        ME = II - LS;
        if (ME > MHT(II)) MHT(II) = ME; end
    end
end

sdata.MHT = MHT;

end