%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Assemble structure stiffness matrix                         *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Mechanics/ADDBAN.m - ADDBAN()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/StrainStiff.m                                 *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************
function AxisymAssemble()
global sdata;
global cdata;

sdata.STIFF = zeros(sdata.NWK, 1, 'double');
sdata.MASS = zeros(sdata.NWK, 1, 'double'); %M����һ���Ĵ���
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
EN = sdata.E; AREA = sdata.AREA; LM = sdata.LM; NINT = sdata.NINT;
PRN = sdata.PR; IOUT = cdata.IOUT;density = sdata.density;
masschoose = sdata.masschoose;

for N = 1:NUME
    S = zeros(16, 16, 'double');  %S�ǸնȾ��󣬽�8*8�ĳ�16*16
    M = zeros(16, 16, 'double');  %M�ǸնȾ��󣬽�8*8�ĳ�16*16
    XX = zeros(2, 8, 'double');    %XX�����꣬��2*4�ĳ�2*8
    B = zeros(4, 16, 'double');   %����Ӧ����󣬴�3*8�ĳ�4*16
    DB = zeros(3, 'double');     %û�������������沢û���õ�
    XG = zeros(4, 4, 'double');  %��˹���ֵ��λ�ã�ÿһ�ж�Ӧ�Ų�ͬ�Ļ��ַ�����������ԳƵ�Ԫ��2*2�ļ������־�����
    WGT = zeros(4, 4, 'double');   %��˹���ֵ��Ȩ�أ����ǲ�ȡ�������ַ�����ÿһ�����ֵ��Ȩ����1
    MTYPE = MATP(N);
    NU = PRN(MTYPE);
    E = EN(MTYPE);
    for i=1:2
        for j=1:8   %��4�ĳ���8����Ϊ��8���ڵ�
            XX(i,j) = XYZ((j-1)*2+i,N);  %xx���������꣬��һ�������꣬�ڶ�����y����
        end
    end   %�������XXû����
    
    %���ɵ�Ԫ��ƽ���뾶
    sumrow = sum(XX,2); %�н飬Ϊ�˼���ƽ���뾶����һ��
    rave = sumrow(1,1)/8;  %��ԳƵ�Ԫ�Ƚ����⣬��Ҫ��Ԫ��ƽ���뾶,�������ǳ�8�ˣ����ǲ�Ӱ����
    
if N==4 %�ϵ�  %�ϵ㣿�򲻳�������ܣ����Լ������£��һ���Ϊ�������˳�ѭ���ģ��������������ԵĹ��ܰ�
end
    
    XG=[0.0  -0.5773502691896  -0.7745966692415  -0.8611363115941;...
        0.0   0.5773502691896   0.0              -0.3399810435849;...
        0.0   0.0               0.7745966692415   0.3399810435849;...
        0.0   0.0               0.0               0.8611363115941];
    WGT=[2.0, 1.0, 0.5555555555555, 0.3478548451375;...
        0.0, 1.0, 0.8888888888888, 0.6521451548625;...
        0.0, 0.0, 0.5555555555555, 0.6521451548625;...
        0.0, 0.0, 0.0,             0.3478548451375];
    
    D=E/(1+NU)/(1-2*NU)*[1-NU NU 0 NU;...
                         NU 1-NU 0 NU;...
                         0 0 (1-2*NU)/2 0;...
                         NU NU 0 1-NU];
    point_position = [-1 -1;...    %���ǽڵ�ĵȲ�Ԫ����
                       1 -1;...
                       1  1;...
                      -1  1;...
                       0 -1;...
                       1  0;...
                       0  1;...
                      -1  0;];                  
    weight_point = [-0.3333, -0.3333, -0.3333, -0.3333, 1.3333, 1.3333, 1.3333, 1.3333] ;                 
    
    for LX=1:NINT  
        RI = XG(LX,NINT); %����˹���ֵ�ȡ��������Ĺ���������ɨ��
     
        for LY=1:NINT
            SI=XG(LY,NINT);
            
            [B,DET,rphy]=generateBAD(RI,SI,XX);  %ͨ�������˹���ֵ��λ�ã��õ�ÿ����˹���ֵ��B������ʽ�������XX���ж��������Bû����
            
            WT = WGT(LX,NINT)*WGT(LY,NINT)*DET*2*pi*rphy;%�õ���Ȩ��
            S=S+B'*D*B*WT;    
        end   
    end
   %��������Ҫ���ǻ��ֵ��Ȩ��
   if masschoose == 3; %�����ֽڵ���ֵķ�ʽ������������
    for L = 1 : 8  
          [ Nmatrix, DET] = generateN( point_position(L, 1), point_position(L, 2), XX ) ;
          M  = M + Nmatrix'* Nmatrix * DET * weight_point(L);
    end
    M = 2 * pi * rave * density * M;
    
    
   elseif masschoose == 1; %��һ������ӵķ�ʽ������������
    for LX=1:NINT  
        RI = XG(LX,NINT); %����˹���ֵ�ȡ��������Ĺ���������ɨ��
     
        for LY=1:NINT
            SI=XG(LY,NINT);
            
            [ Nmatrix, DET] = generateN( RI, SI, XX ) ;  %ͨ�������˹���ֵ��λ�ã��õ�ÿ����˹���ֵ��B������ʽ�������XX���ж��������Bû����
            WT = WGT(LX,NINT)*WGT(LY,NINT);
            M  = M + Nmatrix'* Nmatrix * DET * WT;%�õ���Ȩ��
        end   
    end
     M = 2 * pi * rave * density * M;
     M = diag(sum(M, 2)); 
   end
    
    
    %%��������ԳƵ�Ԫ�����Ե�Ԫ��K������Ҫ����2*pi*r��ѭ���г˹��ˣ��Ժ���Գ���2pi�������ٳ����ѭ����
    ADDBANMASS(M, LM(:, N));
    ADDBAN(S, LM(:, N));
end    

% The third time stamp
cdata.TIM(3, :) = clock;

end




