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
<<<<<<< HEAD
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
EN = sdata.E; LM = sdata.LM; NNU = sdata.PR;
U = sdata.DIS(:, NUM);
XX = zeros(2, 8, 'double');   %½«4¸Ä³ÉÁË8
UV = zeros(1, 16, 'double');  %ÕâÓ¦¸ÃÊÇÓÃÀ´´¢´æÎ»ÒÆµÄÒ»¸öÖĞ½éÁ¿£¬½«8¸Ä³ÉÁË16
B = zeros(4, 16, 'double');   %½«8¸Ä³ÉÁË16£¬3¸Ä³ÉÁË4
maxx = 0.0;
maxxn = 0;  %Öá¶Ô³Æµ¥Ôª²¢ÇÒÊÇÖá¶Ô³ÆÔØºÉ
%Êä³öÔöÌíÁËÖá¶Ô³ÆµÄÓ¦Á¦
fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT             SoildSR            SoildSZ'...
    '            SoildSRZ            SoildStheta\n' ...
    '       NUMBER\n'], NG);


for N = 1:NUME
    MTYPE = MATP(N);
    NU = NNU(MTYPE);
    E = EN(MTYPE);
    D=E/(1+NU)/(1-2*NU)*[1-NU NU 0 NU;NU 1-NU 0 NU;0 0 (1-2*NU)/2 0;NU NU 0 1-NU]; %¸Ä³ÉÁËÖá¶Ô³Æµ¥ÔªµÄD
=======
MODEX = cdata.MODEX;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
EN = sdata.E; LM = sdata.LM; NNU = sdata.PR;
U = sdata.DIS(:, NUM);
sdata.STRAIN = zeros(NUME, 4, 'double');
sdata.STRESS = zeros(NUME, 5, 'double');
STRAIN = sdata.STRAIN;
STRESS = sdata.STRESS;
XX = zeros(2, 8, 'double');   %½«4¸Ä³ÉÁË8
UV = zeros(1, 16, 'double');  %ÕâÓ¦¸ÃÊÇÓÃÀ´´¢´æÎ»ÒÆµÄÒ»¸öÖĞ½éÁ¿£¬½«8¸Ä³ÉÁË16
 %Öá¶Ô³Æµ¥Ôª²¢ÇÒÊÇÖá¶Ô³ÆÔØºÉ
%Êä³öÔöÌíÁËÖá¶Ô³ÆµÄÓ¦Á¦
if (MODEX == 1)
    fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
        'E L E M E N T  G R O U P %4d\n\n' ...
        '       ELEMENT             SoildSR            SoildSZ'...
        '            SoildSRZ            SoildStheta  \n' ...
        '       NUMBER\n'], NG);
end
MTYPE = MATP(1);
NU = NNU(MTYPE);
E = EN(MTYPE);
D=E/(1+NU)/(1-2*NU)*[1-NU NU 0 NU;NU 1-NU 0 NU;0 0 (1-2*NU)/2 0;NU NU 0 1-NU]; %¸Ä³ÉÁËÖá¶Ô³Æµ¥ÔªµÄD
for N = 1:NUME
>>>>>>> å¢æ·»äº†ä¸‰ç§è´¨é‡çŸ©é˜µç”Ÿæˆæ–¹å¼
    UV (1,1:16) =0;%½«8¸Ä³ÉÁË16
    for i = 1:2
        for j = 1:8 %½«4¸Ä³ÉÁË8
            XX(i,j) = XYZ(2*j+i-2,N);
        end
    end
<<<<<<< HEAD
    
    for i = 1:16 %½«8¸Ä³ÉÁË16
        if LM(i, N) > 0
        UV(1,i) = U(LM(i, N));
        end
    end
   % [B1,~]=generateBAD(-1, -1,XX,B);   %JDDµÄË³ĞòºÍÎÒ²»Ò»Ñù£¬ÎÒ´Ó×óÏÂ½Ç¿ªÊ¼Ñ­»·£¬5-8´Óµ×±ß¿ªÊ¼Ñ­»·£¬¶¼ÊÇÄæÊ±ÕëÑ­»·
   % [B2,~]=generateBAD(1, -1,XX,B);
    %[B3,~]=generateBAD(1, 1,XX,B);
    %[B4,~]=generateBAD(-1, -1,XX,B);
   % [B5,~]=generateBAD(0, -1,XX,B);   %ÔöÌí4¸ö
   % [B6,~]=generateBAD(1, 0,XX,B);
   % [B7,~]=generateBAD(0, 1,XX,B);
    %[B8,~]=generateBAD(-1, 0,XX,B); 
    [B1, ~, ~] = generateBAD (-3^0.5/3, -3^0.5/3, XX);
    [B2, ~, ~] = generateBAD (3^0.5/3, -3^0.5/3, XX);
    [B3, ~, ~] = generateBAD (3^0.5/3, 3^0.5/3, XX);
    [B4, ~, ~] = generateBAD (-3^0.5/3, -3^0.5/3, XX);
    % [B1,~]=STDM(-3^0.5/3,-3^0.5/3,XX,B);  %ÕâÀï¿ÉÒÔÔöÌíÓ¦Á¦ÖØ¹¹£¬ÓÃ¸ßË¹µã¾«¶È¸ü¸ß
    % [B2,~]=STDM(3^0.5/3,-3^0.5/3,XX,B);  %JDDµÄË³ĞòºÍÎÒ²»Ò»Ñù£¬ÎÒ´Ó×óÏÂ½Ç¿ªÊ¼Ñ­»·
    % [B3,~]=STDM(3^0.5/3,3^0.5/3,XX,B);
    % [B4,~]=STDM(-3^0.5/3,-3^0.5/3,XX,B);
    STR1 = D*(B1*UV'); %Ã¿¸öµãÓĞ4¸öÓ¦Á¦
    STR2 = D*(B2*UV');
    STR3 = D*(B3*UV');
    STR4 = D*(B4*UV');
    STR = (STR1+STR2+STR3+STR4)/4;
    fprintf(IOUT, ' %10d           %13.6e    %13.6e    %13.6e    %13.6e\n'...
        , N, STR(1), STR(2), STR(3), STR(4));
   %ÔöÌíÁË¸ßË¹µãµÄÓ¦Á¦Êä³ö£¬Êä³ö´Ó×óÏÂ½ÇÄæÊ±ÕëÑ­»·¿ªÊ¼
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
      %  if(POINTSX(ii)>maxx)%×î´óx·½ÏòÓ¦Á¦ÓÃÒÔÑéÖ¤Ğ¡¿×Ó¦Á¦¼¯ÖĞ
       %     maxx = max(POINTSX);
       %     maxcenx=XYZ(2*ii-1,N);
       %     maxceny=XYZ(2*ii,N);
      %  end
   % end
%end

%fprintf(IOUT, ['\n\n  MAX SolidSX  \n\n' ... 
%'       SolidSX            ½Úµã×ø±ê\n\n']);
%fprintf(IOUT, '    %13.6e      (%5.4f , %5.4f) \n',...
%     maxx, maxcenx, maxceny);
    
=======
    for i = 1:16 %½«8¸Ä³ÉÁË16
        if LM(i, N) > 0
            UV(1,i) = U(LM(i, N));
        end
    end
    [B1, ~, ~] = generateBAD (1, 1, XX);
    [B2, ~, ~] = generateBAD (2, 1, XX);
    [B3, ~, ~] = generateBAD (2, 2, XX);
    [B4, ~, ~] = generateBAD (1, 1, XX);
    STI1 = (B1*UV');
    STI2 = (B2*UV');
    STI3 = (B3*UV');
    STI4 = (B4*UV');
    
    STR1 = D*STI1; %Ã¿¸öµãÓĞ4¸öÓ¦Á¦
    STR2 = D*STI2;
    STR3 = D*STI3;
    STR4 = D*STI4;
    STI = (STI1+STI2+STI3+STI4)/4;
    STR = (STR1+STR2+STR3+STR4)/4;
    Mises = (((STR(1)-STR(2))^2+(STR(2)-STR(4))^2+(STR(1)-STR(4))^2+6*STR(3)^2)^0.5)/sqrt(2);
    
    if (MODEX == 1)
        fprintf(IOUT, ' %10d           %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , N, STR(1), STR(2), STR(3), STR(4));
        %ÔöÌíÁË¸ßË¹µãµÄÓ¦Á¦Êä³ö£¬Êä³ö´Ó×óÏÂ½ÇÄæÊ±ÕëÑ­»·¿ªÊ¼
        fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , STR1(1), STR1(2), STR1(3), STR1(4));
        fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , STR2(1), STR2(2), STR2(3), STR2(4));
        fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , STR3(1), STR3(2), STR3(3), STR3(4));
        fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , STR4(1), STR4(2), STR4(3), STR4(4));
    end
        STRAIN(N,:) = STI;
        STRESS(N,1:4) = STR;
        STRESS(N,5) = Mises;
end
[maxMises,index] = max(STRESS(:,5));
sdata.maxMises(1) = maxMises;
sdata.maxMises(2) = index;
sdata.STRAIN = STRAIN;
sdata.STRESS = STRESS;
end

>>>>>>> å¢æ·»äº†ä¸‰ç§è´¨é‡çŸ©é˜µç”Ÿæˆæ–¹å¼
