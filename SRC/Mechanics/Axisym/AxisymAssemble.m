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
<<<<<<< HEAD
global cdata;

sdata.STIFF = zeros(sdata.NWK, 1, 'double');
sdata.MASS = zeros(sdata.NWK, 1, 'double'); %M¾ØÕóÒ»ÑùµÄ´¦Àí
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
EN = sdata.E; AREA = sdata.AREA; LM = sdata.LM; NINT = sdata.NINT;
PRN = sdata.PR; IOUT = cdata.IOUT;density = sdata.density;
masschoose = sdata.masschoose;

for N = 1:NUME
    S = zeros(16, 16, 'double');  %SÊÇ¸Õ¶È¾ØÕó£¬½«8*8¸Ä³É16*16
    M = zeros(16, 16, 'double');  %MÊÇ¸Õ¶È¾ØÕó£¬½«8*8¸Ä³É16*16
    XX = zeros(2, 8, 'double');    %XXÊÇ×ø±ê£¬´Ó2*4¸Ä³É2*8
    B = zeros(4, 16, 'double');   %ÕâÊÇÓ¦±ä¾ØÕó£¬´Ó3*8¸Ä³É4*16
    DB = zeros(3, 'double');     %Ã»¿´¶®ÕâÀï£¬ºÃÏñºóÃæ²¢Ã»ÓĞÓÃµ½
    XG = zeros(4, 4, 'double');  %¸ßË¹»ı·ÖµãµÄÎ»ÖÃ£¬Ã¿Ò»ÁĞ¶ÔÓ¦×Å²»Í¬µÄ»ı·Ö·½°¸£¬ÎÒÃÇÖá¶Ô³Æµ¥ÔªÓÃ2*2µÄ¼õËõ»ı·Ö¾ÍĞĞÁË
    WGT = zeros(4, 4, 'double');   %¸ßË¹»ı·ÖµãµÄÈ¨ÖØ£¬ÎÒÃÇ²ÉÈ¡¼õËõ»ı·Ö·½°¸£¬Ã¿Ò»¸ö»ı·ÖµãµÄÈ¨ÖØÊÇ1
=======


sdata.STIFF = zeros(sdata.NDOF^2 * sdata.NNODE^2 * sdata.NUME, 3);
sdata.MASS = zeros(sdata.NDOF^2 * sdata.NNODE^2 * sdata.NUME, 3);
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
EN = sdata.E; LM = sdata.LM; NINT = sdata.NINT;
PRN = sdata.PR; density = sdata.density;
masschoose = sdata.masschoose;
    XX = zeros(2, 8, 'double');    %XXÊÇ×ø±ê£¬´Ó2*4¸Ä³É2*8
    XG=[0.0  -0.5773502691896  -0.7745966692415  -0.8611363115941;...
        0.0   0.5773502691896   0.0              -0.3399810435849;... %¸ßË¹»ı·ÖµãµÄÎ»ÖÃ£¬Ã¿Ò»ÁĞ¶ÔÓ¦×Å²»Í¬µÄ»ı·Ö·½°¸£¬ÎÒÃÇÖá¶Ô³Æµ¥ÔªÓÃ2*2µÄ¼õËõ»ı·Ö¾ÍĞĞÁË
        0.0   0.0               0.7745966692415   0.3399810435849;...
        0.0   0.0               0.0               0.8611363115941];
    WGT=[2.0, 1.0, 0.5555555555555, 0.3478548451375;...%¸ßË¹»ı·ÖµãµÄÈ¨ÖØ
        0.0, 1.0, 0.8888888888888, 0.6521451548625;...
        0.0, 0.0, 0.5555555555555, 0.6521451548625;...
        0.0, 0.0, 0.0,             0.3478548451375];

  intiN(NINT); %³õÊ¼»¯N,BN
for N = 1:NUME
    S = zeros(16, 16, 'double');  %SÊÇ¸Õ¶È¾ØÕó£¬½«8*8¸Ä³É16*16
    M = zeros(16, 16, 'double');  %MÊÇ¸Õ¶È¾ØÕó£¬½«8*8¸Ä³É16*16

>>>>>>> å¢æ·»äº†ä¸‰ç§è´¨é‡çŸ©é˜µç”Ÿæˆæ–¹å¼
    MTYPE = MATP(N);
    NU = PRN(MTYPE);
    E = EN(MTYPE);
    for i=1:2
        for j=1:8   %½«4¸Ä³ÉÁË8£¬ÒòÎªÊÇ8¸ö½Úµã
            XX(i,j) = XYZ((j-1)*2+i,N);  %xxÊÇÎïÀí×ø±ê£¬µÚÒ»ĞĞÊÇ×ø±ê£¬µÚ¶şĞĞÊÇy×ø±ê
        end
    end   %¼ì²éÁËÏÂXXÃ»ÎÊÌâ
    
    %Éú³Éµ¥ÔªµÄÆ½¾ù°ë¾¶
    sumrow = sum(XX,2); %ÖĞ½é£¬ÎªÁË¼ÆËãÆ½¾ù°ë¾¶ÀûÓÃÒ»ÏÂ
    rave = sumrow(1,1)/8;  %Öá¶Ô³Æµ¥Ôª±È½ÏÌØÊâ£¬ĞèÒªµ¥ÔªµÄÆ½¾ù°ë¾¶,ÕâÀïÍü¼Ç³ı8ÁË£¬µ«ÊÇ²»Ó°ÏìÖÈ
<<<<<<< HEAD
    
if N==4 %¶Ïµã  %¶Ïµã£¿´ò²»³ÉÕâ¸ö¹¦ÄÜ£¬ÎÒ×Ô¼ºÊÔÁËÏÂ£¬ÎÒ»¹ÒÔÎªÊÇÓÃÀ´ÍË³öÑ­»·µÄ£¬¿ÉÄÜÊÇÓÃÀ´µ÷ÊÔµÄ¹¦ÄÜ°Ñ
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
    point_position = [-1 -1;...    %ÕâÊÇ½ÚµãµÄµÈ²ÎÔª×ø±ê
                       1 -1;...
                       1  1;...
                      -1  1;...
                       0 -1;...
                       1  0;...
                       0  1;...
                      -1  0;];                  
    weight_point = [-0.3333, -0.3333, -0.3333, -0.3333, 1.3333, 1.3333, 1.3333, 1.3333] ;                 
    
    for LX=1:NINT  
        RI = XG(LX,NINT); %½«¸ßË¹»ı·ÖµãÈ¡³ö£¬ÕâÀïµÄ¹æÂÉÊÇÏÈÁĞÉ¨Ãè
     
        for LY=1:NINT
            SI=XG(LY,NINT);
            
            [B,DET,rphy]=generateBAD(RI,SI,XX);  %Í¨¹ıÊäÈë¸ßË¹»ı·ÖµãµÄÎ»ÖÃ£¬µÃµ½Ã¿¸ö¸ßË¹»ı·ÖµãµÄBºÍĞĞÁĞÊ½£¬ÕâÀïµÄXXÊÇÓĞ¶«Î÷¶î£¬µ«ÊÇBÃ»¶«Î÷
=======

    D=E/(1+NU)/(1-2*NU)*[1-NU NU 0 NU;...
                         NU 1-NU 0 NU;...
                         0 0 (1-2*NU)/2 0;...
                         NU NU 0 1-NU];                
    
    for LX=1:NINT  
     
        for LY=1:NINT
            
            [B,DET,rphy]=generateBAD(LX,LY,XX);  %Í¨¹ıÊäÈë¸ßË¹»ı·ÖµãµÄÎ»ÖÃ£¬µÃµ½Ã¿¸ö¸ßË¹»ı·ÖµãµÄBºÍĞĞÁĞÊ½£¬ÕâÀïµÄXXÊÇÓĞ¶«Î÷¶î£¬µ«ÊÇBÃ»¶«Î÷
>>>>>>> å¢æ·»äº†ä¸‰ç§è´¨é‡çŸ©é˜µç”Ÿæˆæ–¹å¼
            
            WT = WGT(LX,NINT)*WGT(LY,NINT)*DET*2*pi*rphy;%µÃµ½×ÜÈ¨ÖØ
            S=S+B'*D*B*WT;    
        end   
    end
<<<<<<< HEAD
   %ÖÊÁ¿¾ØÕó£¬Òª¿¼ÂÇ»ı·ÖµãµÄÈ¨ÖØ
   if masschoose == 3; %µÚÈıÖÖ½Úµã»ı·ÖµÄ·½Ê½Éú³ÉÖÊÁ¿¾ØÕó
    for L = 1 : 8  
          [ Nmatrix, DET] = generateN( point_position(L, 1), point_position(L, 2), XX ) ;
          M  = M + Nmatrix'* Nmatrix * DET * weight_point(L);
    end
    M = 2 * pi * rave * density * M;
    
    
   elseif masschoose == 1; %µÚÒ»ÖÖĞĞÏà¼ÓµÄ·½Ê½Éú³ÉÖÊÁ¿¾ØÕó
=======
    
   %ÖÊÁ¿¾ØÕó£¬Òª¿¼ÂÇ»ı·ÖµãµÄÈ¨ÖØ
    for LX=1:3  
        RI = XG(LX,3); %½«¸ßË¹»ı·ÖµãÈ¡³ö£¬ÕâÀïµÄ¹æÂÉÊÇÏÈÁĞÉ¨Ãè
     
        for LY=1:3
            SI=XG(LY,3);
            
            [ Nmatrix, DET] = generateN( RI, SI, XX ) ;  %Í¨¹ıÊäÈë¸ßË¹»ı·ÖµãµÄÎ»ÖÃ£¬µÃµ½Ã¿¸ö¸ßË¹»ı·ÖµãµÄBºÍĞĞÁĞÊ½£¬ÕâÀïµÄXXÊÇÓĞ¶«Î÷¶î£¬µ«ÊÇBÃ»¶«Î÷
            WT = WGT(LX,3)*WGT(LY,3);
            M  = M + Nmatrix'* Nmatrix * DET * WT;%µÃµ½×ÜÈ¨ÖØ
        end   
    end
    
   if masschoose == 1
>>>>>>> å¢æ·»äº†ä¸‰ç§è´¨é‡çŸ©é˜µç”Ÿæˆæ–¹å¼
    for LX=1:NINT  
        RI = XG(LX,NINT); %½«¸ßË¹»ı·ÖµãÈ¡³ö£¬ÕâÀïµÄ¹æÂÉÊÇÏÈÁĞÉ¨Ãè
     
        for LY=1:NINT
            SI=XG(LY,NINT);
            
            [ Nmatrix, DET] = generateN( RI, SI, XX ) ;  %Í¨¹ıÊäÈë¸ßË¹»ı·ÖµãµÄÎ»ÖÃ£¬µÃµ½Ã¿¸ö¸ßË¹»ı·ÖµãµÄBºÍĞĞÁĞÊ½£¬ÕâÀïµÄXXÊÇÓĞ¶«Î÷¶î£¬µ«ÊÇBÃ»¶«Î÷
            WT = WGT(LX,NINT)*WGT(LY,NINT);
            M  = M + Nmatrix'* Nmatrix * DET * WT;%µÃµ½×ÜÈ¨ÖØ
        end   
    end
     M = 2 * pi * rave * density * M;
     M = diag(sum(M, 2)); 
<<<<<<< HEAD
   end
    
    
    %%ÕâÀïÊÇÖá¶Ô³Æµ¥Ôª£¬ËùÒÔµ¥ÔªµÄK¾ØÕóĞèÒª³ËÒÔ2*pi*rÔÚÑ­»·ÖĞ³Ë¹ıÁË£¬ÒÔºó¿ÉÒÔ³ËÒÔ2pi£¬À´¼õÉÙ³ÌĞòµÄÑ­»·Á¿
    ADDBANMASS(M, LM(:, N));
    ADDBAN(S, LM(:, N));
end    

% The third time stamp
cdata.TIM(3, :) = clock;

end




=======
        
        
    
    
    
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
    %%ÕâÀïÊÇÖá¶Ô³Æµ¥Ôª£¬ËùÒÔµ¥ÔªµÄK¾ØÕóĞèÒª³ËÒÔ2*pi*rÔÚÑ­»·ÖĞ³Ë¹ıÁË£¬ÒÔºó¿ÉÒÔ³ËÒÔ2pi£¬À´¼õÉÙ³ÌĞòµÄÑ­»·Á¿
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
    RI = XG(LX,NINT); %½«¸ßË¹»ı·ÖµãÈ¡³ö£¬ÕâÀïµÄ¹æÂÉÊÇÏÈÁĞÉ¨Ãè
    for LY=1:NINT
        SI=XG(LY,NINT);
        
        N(1) = 0.25*(1 - RI)*(1 - SI);
        N(2) = 0.25*(1 + RI)*(1 - SI);
        N(3) = 0.25*(1 + RI)*(1 + SI);
        N(4) = 0.25*(1 - RI)*(1 + SI);
        N(5) = 0.5*(1 - RI*RI)*(1 - SI); %½«CQ4ÔöÌí³ÉÁËCQ8
        N(6) = 0.5*(1 + RI)*(1 - SI*SI);
        N(7) = 0.5*(1 - RI*RI)*(1 + SI);
        N(8) = 0.5*(1 - RI)*(1 - SI*SI);
        N(1) = N(1) - 0.5*N(5) - 0.5*N(8);%¶ÔÔ­À´µÄ4¸öĞÎº¯Êı½øĞĞĞŞÕı
        N(2) = N(2) - 0.5*N(6) - 0.5*N(5);
        N(3) = N(3) - 0.5*N(7) - 0.5*N(6);
        N(4) = N(4) - 0.5*N(8) - 0.5*N(7);  %¼ì²éÁËNÃ»ÓĞÎÊÌâ£¬Ò»²¿·Ö³é²é
        
        N11 = 0.25*(1-SI)*(2*RI+SI);%1µ½4ºÅ½ÚµãÊÇ£¬×óÏÂ½Ç¿ªÊ¼ÄæÊ±ÕëÑ­»·£¬5-8ºÅÊÇ´Óµ×±ßµÄÖĞµã¿ªÊ¼ÄæÊ±ÕëÑ­»·
        N12 = 0.25*(1-SI)*(2*RI-SI);
        N13 = 0.25*(1+SI)*(2*RI+SI);
        N14 = 0.25*(1+SI)*(2*RI-SI);
        
        N15 = -RI*(1-SI);  %ÔöÌíÁËÆäÓà4¸ö½Úµã¶ÔpsiÇóµ¼
        N16 = 0.5*(1 - SI*SI);
        N17 = -RI*(1+SI);
        N18 = -0.5*(1 - SI*SI);
        
        N21 = 0.25*(1-RI)*(RI+2*SI);
        N22 = 0.25*(1+RI)*(-RI+2*SI);
        N23 = 0.25*(1+RI)*(RI+2*SI);
        N24 = 0.25*(1-RI)*(-RI+2*SI);
        
        N25 = -0.5*(1-RI*RI);  %ÔöÌíÁËÆäÓà4¸ö½Úµã¶ÔetaÇóµ¼
        N26 = -SI*(1+RI);
        N27 = 0.5*(1-RI*RI);
        N28 = -SI*(1-RI);
        %Éú³ÉµÈ²ÎÔªµÄĞÎº¯ÊıµÄµ¼Êı
        Bpsi = [N11, N12, N13, N14, N15, N16, N17, N18;...
            N21, N22, N23, N24, N25, N26, N27, N28];
        sdata.N(:,LX,LY) = N;
        sdata.BN(:,:,LX,LY) = Bpsi;
    end
end
end
>>>>>>> å¢æ·»äº†ä¸‰ç§è´¨é‡çŸ©é˜µç”Ÿæˆæ–¹å¼
