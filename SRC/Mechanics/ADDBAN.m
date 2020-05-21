%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     To assemble element stiffness into global stiffness         *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Truss/ReadTruss.m - Assemble()                            *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     LeiYang Zhao, Yan Liu, Computational Dynamics Group,        *
%*     School of Aerospace Engineering, Tsinghua University,       *
%*     2019.02.22                                                  *
%*                                                                 *
%* *****************************************************************

<<<<<<< HEAD
function ADDBAN(S, LM)  %S是单元的刚度阵，LM是单元的总体刚度的坐标对应

% Get global data
global sdata;
MAXA = sdata.MAXA; 
STIFF = sdata.STIFF;
ND = sdata.NDOF * sdata.NNODE;
for J = 1:ND
    JJ = LM(J);
    if (JJ > 0)
        for I = 1:J
            II = LM(I);
            if (II > 0)
                if (JJ > II) KK = MAXA(JJ) + JJ - II;
                else KK = MAXA(II) + II - JJ; end
                STIFF(KK) = STIFF(KK) + S(I, J);  %将单元刚度阵集成到总体刚度阵，STIFF只有上三角
            end
        end
    end
end

sdata.STIFF = STIFF;
=======
function ADDBAN(S, LM,N)

% Get global data
global sdata;
STIFF = sdata.STIFF;NUME = sdata.NUME;
ND = sdata.NDOF * sdata.NNODE;
N1 =ND*ND*(N-1);
for J = 1:ND
    JJ = LM(J);
    for I = 1:ND
        II = LM(I);
        STIFF(N1+ND*(J-1)+I,1) =  II;
        STIFF(N1+ND*(J-1)+I,2) =  JJ;
        STIFF(N1+ND*(J-1)+I,3) =  S(I, J);
    end
end
if (N==NUME)
    STIFF(STIFF(:,1)==0,:)=[];
    STIFF(STIFF(:,2)==0,:)=[];
    STIFF = sparse(STIFF(:,1),STIFF(:,2),STIFF(:,3),sdata.NEQ,sdata.NEQ);
    sdata.STIFF = STIFF;
else
    sdata.STIFF = STIFF;
end
>>>>>>> 澧炴坊浜嗕笁绉嶈川閲忕煩闃电敓鎴愭柟寮�

end