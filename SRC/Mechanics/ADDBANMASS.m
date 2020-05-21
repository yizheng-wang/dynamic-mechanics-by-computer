%* *****************************************************************
%* - Function of STAPMAT in MASSness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     To assemble element MASSness into global MASSness         *
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
function ADDBANMASS(M, LM)  %M «µ•‘™µƒ∏’∂»’Û£¨LM «µ•‘™µƒ◊‹ÃÂ∏’∂»µƒ◊¯±Í∂‘”¶

% Get global data
global sdata;
MAXA = sdata.MAXA; 
MASS = sdata.MASS;
ND = sdata.NDOF * sdata.NNODE;
for J = 1:ND
    JJ = LM(J);
    if (JJ > 0)
        for I = 1:J
            II = LM(I);
            if (II > 0)
                if (JJ > II) KK = MAXA(JJ) + JJ - II;
                else KK = MAXA(II) + II - JJ; end
                MASS(KK) = MASS(KK) + M(I, J);  %Ω´µ•‘™∏’∂»’ÛºØ≥…µΩ◊‹ÃÂ∏’∂»’Û£¨MASS÷ª”–…œ»˝Ω«
            end
        end
    end
end

sdata.MASS = MASS;
=======
function ADDBANMASS(M, LM,N)  %M «µ•‘™µƒ∏’∂»’Û£¨LM «µ•‘™µƒ◊‹ÃÂ∏’∂»µƒ◊¯±Í∂‘”¶

% Get global data
global sdata;
MASS = sdata.MASS;NUME = sdata.NUME;
ND = sdata.NDOF * sdata.NNODE;
N1 = ND*ND*(N-1);
for J = 1:ND
    JJ = LM(J);
    for I = 1:ND
        II = LM(I);
        MASS(N1+ND*(I-1)+J,1) =  II;
        MASS(N1+ND*(I-1)+J,2) =  JJ;
        MASS(N1+ND*(I-1)+J,3) =  M(I, J);
    end
end
if (N==NUME)
    MASS(MASS(:,1)==0,:)=[];
    MASS(MASS(:,2)==0,:)=[];
    MASS = sparse(MASS(:,1),MASS(:,2),MASS(:,3),sdata.NEQ,sdata.NEQ);
    sdata.MASS = MASS;
else
    sdata.MASS = MASS;
end
>>>>>>> Â¢ûÊ∑ª‰∫Ü‰∏âÁßçË¥®ÈáèÁü©ÈòµÁîüÊàêÊñπÂºè

end