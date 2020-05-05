%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function AxisymStress(NUM,NG)

% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
EN = sdata.E; LM = sdata.LM; NNU = sdata.PR;
U = sdata.DIS(:, NUM);
XX = zeros(2, 8, 'double');   %��4�ĳ���8
UV = zeros(1, 16, 'double');  %��Ӧ������������λ�Ƶ�һ���н�������8�ĳ���16
B = zeros(4, 16, 'double');   %��8�ĳ���16��3�ĳ���4
maxx = 0.0;
maxxn = 0;  %��ԳƵ�Ԫ��������Գ��غ�
%�����������ԳƵ�Ӧ��
fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT             SoildSR            SoildSZ'...
    '            SoildSRZ            SoildStheta\n' ...
    '       NUMBER\n'], NG);


for N = 1:NUME
    MTYPE = MATP(N);
    NU = NNU(MTYPE);
    E = EN(MTYPE);
    D=E/(1+NU)/(1-2*NU)*[1-NU NU 0 NU;NU 1-NU 0 NU;0 0 (1-2*NU)/2 0;NU NU 0 1-NU]; %�ĳ�����ԳƵ�Ԫ��D
    UV (1,1:16) =0;%��8�ĳ���16
    for i = 1:2
        for j = 1:8 %��4�ĳ���8
            XX(i,j) = XYZ(2*j+i-2,N);
        end
    end
    
    for i = 1:16 %��8�ĳ���16
        if LM(i, N) > 0
        UV(1,i) = U(LM(i, N));
        end
    end
   % [B1,~]=generateBAD(-1, -1,XX,B);   %JDD��˳����Ҳ�һ�����Ҵ����½ǿ�ʼѭ����5-8�ӵױ߿�ʼѭ����������ʱ��ѭ��
   % [B2,~]=generateBAD(1, -1,XX,B);
    %[B3,~]=generateBAD(1, 1,XX,B);
    %[B4,~]=generateBAD(-1, -1,XX,B);
   % [B5,~]=generateBAD(0, -1,XX,B);   %����4��
   % [B6,~]=generateBAD(1, 0,XX,B);
   % [B7,~]=generateBAD(0, 1,XX,B);
    %[B8,~]=generateBAD(-1, 0,XX,B); 
    [B1, ~, ~]=generateBAD (-3^0.5/3, -3^0.5/3, XX, B)
    [B2, ~, ~]=generateBAD (3^0.5/3, -3^0.5/3, XX, B)
    [B3, ~, ~]=generateBAD (3^0.5/3, 3^0.5/3, XX, B)
    [B4, ~, ~]=generateBAD (-3^0.5/3, -3^0.5/3, XX, B)
    % [B1,~]=STDM(-3^0.5/3,-3^0.5/3,XX,B);  %�����������Ӧ���ع����ø�˹�㾫�ȸ���
    % [B2,~]=STDM(3^0.5/3,-3^0.5/3,XX,B);  %JDD��˳����Ҳ�һ�����Ҵ����½ǿ�ʼѭ��
    % [B3,~]=STDM(3^0.5/3,3^0.5/3,XX,B);
    % [B4,~]=STDM(-3^0.5/3,-3^0.5/3,XX,B);
    STR1 = D*(B1*UV'); %ÿ������4��Ӧ��
    STR2 = D*(B2*UV');
    STR3 = D*(B3*UV');
    STR4 = D*(B4*UV');
    STR = (STR1+STR2+STR3+STR4)/4;
    fprintf(IOUT, ' %10d           %13.6e    %13.6e    %13.6e    %13.6e\n'...
        , N, STR(1), STR(2), STR(3), STR(4));
   %�����˸�˹���Ӧ���������������½���ʱ��ѭ����ʼ
    fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
        , STR1(1), STR1(2), STR1(3), STR1(4));
    fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
        , STR2(1), STR2(2), STR2(3), STR2(4));   
    fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
        , STR3(1), STR3(2), STR3(3), STR3(4));    
    fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
        , STR4(1), STR4(2), STR4(3), STR4(4));    
end
    %POINTSX(1) = STR1(1);
   % POINTSX(2) = STR2(1);
    %POINTSX(3) = STR3(1);
    %POINTSX(4) = STR4(1);
    %for ii = 1:4
      %  if(POINTSX(ii)>maxx)%���x����Ӧ��������֤С��Ӧ������
       %     maxx = max(POINTSX);
       %     maxcenx=XYZ(2*ii-1,N);
       %     maxceny=XYZ(2*ii,N);
      %  end
   % end
%end

%fprintf(IOUT, ['\n\n  MAX SolidSX  \n\n' ... 
%'       SolidSX            �ڵ�����\n\n']);
%fprintf(IOUT, '    %13.6e      (%5.4f , %5.4f) \n',...
%     maxx, maxcenx, maxceny);
    