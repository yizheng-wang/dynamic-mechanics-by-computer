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
sdata.MASS = zeros(sdata.NWK, 1, 'double'); %M矩阵一样的处理
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
EN = sdata.E; AREA = sdata.AREA; LM = sdata.LM; NINT = sdata.NINT;
PRN = sdata.PR; IOUT = cdata.IOUT;density = sdata.density;
masschoose = sdata.masschoose;

for N = 1:NUME
    S = zeros(16, 16, 'double');  %S是刚度矩阵，将8*8改成16*16
    M = zeros(16, 16, 'double');  %M是刚度矩阵，将8*8改成16*16
    XX = zeros(2, 8, 'double');    %XX是坐标，从2*4改成2*8
    B = zeros(4, 16, 'double');   %这是应变矩阵，从3*8改成4*16
    DB = zeros(3, 'double');     %没看懂这里，好像后面并没有用到
    XG = zeros(4, 4, 'double');  %高斯积分点的位置，每一列对应着不同的积分方案，我们轴对称单元用2*2的减缩积分就行了
    WGT = zeros(4, 4, 'double');   %高斯积分点的权重，我们采取减缩积分方案，每一个积分点的权重是1
    MTYPE = MATP(N);
    NU = PRN(MTYPE);
    E = EN(MTYPE);
    for i=1:2
        for j=1:8   %将4改成了8，因为是8个节点
            XX(i,j) = XYZ((j-1)*2+i,N);  %xx是物理坐标，第一行是坐标，第二行是y坐标
        end
    end   %检查了下XX没问题
    
    %生成单元的平均半径
    sumrow = sum(XX,2); %中介，为了计算平均半径利用一下
    rave = sumrow(1,1)/8;  %轴对称单元比较特殊，需要单元的平均半径,这里忘记除8了，但是不影响秩
    
if N==4 %断点  %断点？打不成这个功能，我自己试了下，我还以为是用来退出循环的，可能是用来调试的功能把
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
    point_position = [-1 -1;...    %这是节点的等参元坐标
                       1 -1;...
                       1  1;...
                      -1  1;...
                       0 -1;...
                       1  0;...
                       0  1;...
                      -1  0;];                  
    weight_point = [-0.3333, -0.3333, -0.3333, -0.3333, 1.3333, 1.3333, 1.3333, 1.3333] ;                 
    
    for LX=1:NINT  
        RI = XG(LX,NINT); %将高斯积分点取出，这里的规律是先列扫描
     
        for LY=1:NINT
            SI=XG(LY,NINT);
            
            [B,DET,rphy]=generateBAD(RI,SI,XX);  %通过输入高斯积分点的位置，得到每个高斯积分点的B和行列式，这里的XX是有东西额，但是B没东西
            
            WT = WGT(LX,NINT)*WGT(LY,NINT)*DET*2*pi*rphy;%得到总权重
            S=S+B'*D*B*WT;    
        end   
    end
   %质量矩阵，要考虑积分点的权重
   if masschoose == 3; %第三种节点积分的方式生成质量矩阵
    for L = 1 : 8  
          [ Nmatrix, DET] = generateN( point_position(L, 1), point_position(L, 2), XX ) ;
          M  = M + Nmatrix'* Nmatrix * DET * weight_point(L);
    end
    M = 2 * pi * rave * density * M;
    
    
   elseif masschoose == 1; %第一种行相加的方式生成质量矩阵
    for LX=1:NINT  
        RI = XG(LX,NINT); %将高斯积分点取出，这里的规律是先列扫描
     
        for LY=1:NINT
            SI=XG(LY,NINT);
            
            [ Nmatrix, DET] = generateN( RI, SI, XX ) ;  %通过输入高斯积分点的位置，得到每个高斯积分点的B和行列式，这里的XX是有东西额，但是B没东西
            WT = WGT(LX,NINT)*WGT(LY,NINT);
            M  = M + Nmatrix'* Nmatrix * DET * WT;%得到总权重
        end   
    end
     M = 2 * pi * rave * density * M;
     M = diag(sum(M, 2)); 
   end
    
    
    %%这里是轴对称单元，所以单元的K矩阵需要乘以2*pi*r在循环中乘过了，以后可以乘以2pi，来减少程序的循环量
    ADDBANMASS(M, LM(:, N));
    ADDBAN(S, LM(:, N));
end    

% The third time stamp
cdata.TIM(3, :) = clock;

end




