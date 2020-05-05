%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of truss                       *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadTruss.m - ReadMaterial()                                *
%*     ReadTruss.m - ReadElements()                                *
%*                                                                 *
%* - Called by :                                                   *
%*     ./StrainStiff.m                                             *
%*                                                                 *
%* *****************************************************************

function AxisymReadStrain()

% Read Material information
ReadMaterial()

% Read Element information
ReadElements()

% the second time stamp
global cdata;
cdata.TIM(2,:) = clock;

end

% ----------------------- Functions -----------------------------------
% Read Material information
function ReadMaterial()

global cdata;
global sdata;
% Get file pointers
IIN = cdata.IIN;
IOUT = cdata.IOUT;

if (cdata.NPAR(3) == 0) 
    cdata.NPAR(3) = 1; 
end
fprintf(IOUT, '\n M A T E R I A L   D E F I N I T I O N\n');
fprintf(IOUT, '\n NUMBER OF DIFFERENT SETS OF MATERIAL\n');
fprintf(IOUT, ' AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . = %10d\n', ...
    cdata.NPAR(3));
fprintf(IOUT, '  SET       YOUNG''S       POISSON''S         GAUSS NUMERICAL  \n');
fprintf(IOUT, ' NUMBER     MODULUS         RATIO       INTEGRATION ORDER \n');
fprintf(IOUT, '               E              PR                  NINT \n');


% Read material datas
sdata.NUME = cdata.NPAR(2);
sdata.NUMMAT = cdata.NPAR(3);
NUMMAT = cdata.NPAR(3);
sdata.E = zeros(NUMMAT, 1, 'double');
sdata.AREA = zeros(NUMMAT, 1, 'double');
for I = 1:cdata.NPAR(3)
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E(N) = tmp(2);
    sdata.PR(N) = tmp(3);
    sdata.NINT(N) = tmp(4);  %��ƽ�浥Ԫ�ĺ��ɾ�������������˹���ֵ㣬�����õ��Ǽ�������
 %�����fprintfҲ�ͽ������޸ģ������������ļ���Ԫ�أ���tmp��5��ɾ����
    fprintf(IOUT, '%5d    %12.5e        %2.1f       %14.6e \n', N, tmp(2), tmp(3), tmp(4));
end

end   %

% Read elements information
function ReadElements()

global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
IOUT = cdata.IOUT;

fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n'); %��ԭ���Ļ����������4���ڵ�
fprintf(IOUT, '\n      ELEMENT          NODE        NODE        NODE        NODE              NODE        NODE        NODE        NODE      MATERIAL\n');
fprintf(IOUT, '      NUMBER-N           I           J           K           L          A          B           C         D       SET NUMBER\n');

% Get Position data
NUME = cdata.NPAR(2);
sdata.XYZ = zeros(16, NUME, 'double');   %��8�ĳ���16
sdata.MATP = zeros(NUME, 1, 'int64');                 % the type of material
sdata.LM = zeros(16, NUME, 'int64');                 %��8�ĳ���16 % connectivity matrix
sdata.MHT = zeros(sdata.NEQ, 1, 'int64');
X = sdata.X; Y = sdata.Y; ID = sdata.ID;
XYZ = sdata.XYZ; MATP = sdata.MATP; LM = sdata.LM;


for N = 1:NUME
    tmp = str2num(fgetl(IIN));
    I = round(tmp(2));
    J = round(tmp(3));
    K = round(tmp(4));
    L = round(tmp(5)); 
    A = round(tmp(6));%������4���ڵ�,ABCD
    B = round(tmp(7));
    C = round(tmp(8));
    D = round(tmp(9));
    MTYPE = round(tmp(10)); %��6�ĳ���10
    
%   Save element information
    XYZ(1, N) = X(I);
    XYZ(2, N) = Y(I);
    XYZ(3, N) = X(J);
    XYZ(4, N) = Y(J);
    XYZ(5, N) = X(K);
    XYZ(6, N) = Y(K);
    XYZ(7, N) = X(L);
    XYZ(8, N) = Y(L);
    XYZ(9, N) = X(A);   %��ԭ����8��ӳ���16
    XYZ(10, N) = Y(A);
    XYZ(11, N) = X(B);
    XYZ(12, N) = Y(B);
    XYZ(13, N) = X(C);
    XYZ(14, N) = Y(C);
    XYZ(15, N) = X(D);
    XYZ(16, N) = Y(D);
    MATP(N) = MTYPE;
    %�����fprintfҲ�������޸�
    fprintf(IOUT, '%10d      %10d  %10d  %10d  %10d %10d  %10d  %10d  %10d       %5d\n', N, I, J, K, L, A, B, C, D  ,MTYPE);

%   Compute connectivity matrix
    LM(1, N) = ID(1, I);
    LM(3, N) = ID(1, J);
    LM(5, N) = ID(1, K);
    LM(7, N) = ID(1, L);
    LM(9, N) = ID(1, A);   %ID�ĵ�һ����X��������ɶȶ�Ӧ������նȾ����λ�ã�����4��
    LM(11, N) = ID(1, B);
    LM(13, N) = ID(1, C);
    LM(15, N) = ID(1, D);
    
    LM(2, N) = ID(2, I);
    LM(4, N) = ID(2, J);
    LM(6, N) = ID(2, K);
    LM(8, N) = ID(2, L);  %ID�ĵڶ�����Y��������ɶȶ�Ӧ������նȾ����λ�ã�����4��
    LM(10, N) = ID(2, A);
    LM(12, N) = ID(2, B);
    LM(14, N) = ID(2, C);
    LM(16, N) = ID(2, D); 
%   Updata column heights and bandwidth
    ColHt(LM(:, N))
end
sdata.XYZ = XYZ; sdata.MATP = MATP; sdata.LM = LM;

% Clear the memory of X, Y, Z
sdata.X = double(0);
sdata.Y = double(0);
sdata.Z = double(0);

end