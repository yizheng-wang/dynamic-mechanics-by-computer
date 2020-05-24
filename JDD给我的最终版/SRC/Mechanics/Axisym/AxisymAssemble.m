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


sdata.STIFF = zeros(sdata.NDOF^2 * sdata.NNODE^2 * sdata.NUME, 3);
sdata.MASS = zeros(sdata.NDOF^2 * sdata.NNODE^2 * sdata.NUME, 3);
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
EN = sdata.E; LM = sdata.LM; NINT = sdata.NINT;
PRN = sdata.PR; density = sdata.density;
masschoose = sdata.masschoose;
    XX = zeros(2, 8, 'double');    %XX�����꣬��2*4�ĳ�2*8
    XG=[0.0  -0.5773502691896  -0.7745966692415  -0.8611363115941;...
        0.0   0.5773502691896   0.0              -0.3399810435849;... %��˹���ֵ��λ�ã�ÿһ�ж�Ӧ�Ų�ͬ�Ļ��ַ�����������ԳƵ�Ԫ��2*2�ļ������־�����
        0.0   0.0               0.7745966692415   0.3399810435849;...
        0.0   0.0               0.0               0.8611363115941];
    WGT=[2.0, 1.0, 0.5555555555555, 0.3478548451375;...%��˹���ֵ��Ȩ��
        0.0, 1.0, 0.8888888888888, 0.6521451548625;...
        0.0, 0.0, 0.5555555555555, 0.6521451548625;...
        0.0, 0.0, 0.0,             0.3478548451375];

  intiN(NINT); %��ʼ��N,BN
for N = 1:NUME
    S = zeros(16, 16, 'double');  %S�ǸնȾ��󣬽�8*8�ĳ�16*16
    M = zeros(16, 16, 'double');  %M�ǸնȾ��󣬽�8*8�ĳ�16*16

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

    D=E/(1+NU)/(1-2*NU)*[1-NU NU 0 NU;...
                         NU 1-NU 0 NU;...
                         0 0 (1-2*NU)/2 0;...
                         NU NU 0 1-NU];                
    
    for LX=1:NINT  
     
        for LY=1:NINT
            
            [B,DET,rphy]=generateBAD(LX,LY,XX);  %ͨ�������˹���ֵ��λ�ã��õ�ÿ����˹���ֵ��B������ʽ�������XX���ж��������Bû����
            
            WT = WGT(LX,NINT)*WGT(LY,NINT)*DET*2*pi*rphy;%�õ���Ȩ��
            S=S+B'*D*B*WT;    
        end   
    end
    
   %��������Ҫ���ǻ��ֵ��Ȩ��
    for LX=1:3  
        RI = XG(LX,3); %����˹���ֵ�ȡ��������Ĺ���������ɨ��
     
        for LY=1:3
            SI=XG(LY,3);
            
            [ Nmatrix, DET] = generateN( RI, SI, XX ) ;  %ͨ�������˹���ֵ��λ�ã��õ�ÿ����˹���ֵ��B������ʽ�������XX���ж��������Bû����
            WT = WGT(LX,3)*WGT(LY,3);
            M  = M + Nmatrix'* Nmatrix * DET * WT;%�õ���Ȩ��
        end   
    end
    
   if masschoose == 1
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
        
        
    
    
    
   elseif masschoose == 2;
    M = 2 * pi * rave * density * M;
    midm1 = sum(M(:));
    M = diag(diag(M)); 
    M = M/(sum(M(:))/midm1);
    
    
   elseif masschoose == 3;
    for L = 1 : 8  
          [ Nmatrix, DET] = generateN( point_position(L, 1), point_position(L, 2), XX ) ;
          M  = M + Nmatrix'* Nmatrix * DET * weight_point(L);
    end
    M = 2 * pi * rave * density * M;
       
    end
%     
    %%��������ԳƵ�Ԫ�����Ե�Ԫ��K������Ҫ����2*pi*r��ѭ���г˹��ˣ��Ժ���Գ���2pi�������ٳ����ѭ����
    ADDBANMASS(M, LM(:, N),N);
    ADDBAN(S, LM(:, N),N);

end    

% spy(sdata.STIFF)
end


function intiN(NINT)
global sdata;
sdata.N = zeros(8,NINT,NINT);
sdata.BN = zeros(2,8,NINT,NINT);
XG=[0.0  -0.5773502691896  -0.7745966692415  -0.8611363115941;...
    0.0   0.5773502691896   0.0              -0.3399810435849;...
    0.0   0.0               0.7745966692415   0.3399810435849;...
    0.0   0.0               0.0               0.8611363115941];
for LX=1:NINT
    RI = XG(LX,NINT); %����˹���ֵ�ȡ��������Ĺ���������ɨ��
    for LY=1:NINT
        SI=XG(LY,NINT);
        
        N(1) = 0.25*(1 - RI)*(1 - SI);
        N(2) = 0.25*(1 + RI)*(1 - SI);
        N(3) = 0.25*(1 + RI)*(1 + SI);
        N(4) = 0.25*(1 - RI)*(1 + SI);
        N(5) = 0.5*(1 - RI*RI)*(1 - SI); %��CQ4�������CQ8
        N(6) = 0.5*(1 + RI)*(1 - SI*SI);
        N(7) = 0.5*(1 - RI*RI)*(1 + SI);
        N(8) = 0.5*(1 - RI)*(1 - SI*SI);
        N(1) = N(1) - 0.5*N(5) - 0.5*N(8);%��ԭ����4���κ�����������
        N(2) = N(2) - 0.5*N(6) - 0.5*N(5);
        N(3) = N(3) - 0.5*N(7) - 0.5*N(6);
        N(4) = N(4) - 0.5*N(8) - 0.5*N(7);  %�����Nû�����⣬һ���ֳ��
        
        N11 = 0.25*(1-SI)*(2*RI+SI);%1��4�Žڵ��ǣ����½ǿ�ʼ��ʱ��ѭ����5-8���Ǵӵױߵ��е㿪ʼ��ʱ��ѭ��
        N12 = 0.25*(1-SI)*(2*RI-SI);
        N13 = 0.25*(1+SI)*(2*RI+SI);
        N14 = 0.25*(1+SI)*(2*RI-SI);
        
        N15 = -RI*(1-SI);  %����������4���ڵ��psi��
        N16 = 0.5*(1 - SI*SI);
        N17 = -RI*(1+SI);
        N18 = -0.5*(1 - SI*SI);
        
        N21 = 0.25*(1-RI)*(RI+2*SI);
        N22 = 0.25*(1+RI)*(-RI+2*SI);
        N23 = 0.25*(1+RI)*(RI+2*SI);
        N24 = 0.25*(1-RI)*(-RI+2*SI);
        
        N25 = -0.5*(1-RI*RI);  %����������4���ڵ��eta��
        N26 = -SI*(1+RI);
        N27 = 0.5*(1-RI*RI);
        N28 = -SI*(1-RI);
        %���ɵȲ�Ԫ���κ����ĵ���
        Bpsi = [N11, N12, N13, N14, N15, N16, N17, N18;...
            N21, N22, N23, N24, N25, N26, N27, N28];
        sdata.N(:,LX,LY) = N;
        sdata.BN(:,:,LX,LY) = Bpsi;
    end
end
end
