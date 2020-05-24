function newmark(L,~)
beta = 1/12;
global sdata;

NN = sdata.NEQ;
% X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID; NUMNP = cdata.NUMNP;%输出使用的数据
K = sdata.STIFF; M = sdata.MASS; 
C = sparse(NN, NN);

Qq = sdata.R(:,L)*2*3.14*30; %固定载荷 也可改为在循环中随时间变化的载荷
Q = Qq;
%固定载荷 也可改为在循环中随时间变化的载荷

TIM1 = clock;%计时

% beta = 0.25;%newmark参数
gamma = 1/2;

dt = 0.00001;
n = 200;%总时间步数，也可改为总时间控制


%计算参数
c0 = 1/(beta*dt^2);
c1 = gamma/(beta*dt);
c2 = 1/(beta*dt);
c3 = 1/(2*beta) - 1;
c4 = gamma/beta - 1;
c5 = dt*(gamma/(2*beta) - 1);
c6 = dt*(1 - gamma);
c7 = gamma*dt;

%形成有效刚度阵
K_eff = K + c0*M +c1*C;
 K_eff = sparse(K_eff);
%三角分解
[L1,D1,P] = ldl(K_eff);

%初始化
xd = zeros(NN,1, 'double');%迭代过程中使用的位移，速度，加速度，存储矩阵
d1 = xd;
v1 = zeros(NN,1, 'double');
a1 = M\(Q-K*d1); %V0 = 0


%t = zeros(n,1);
%t(1) = 0;



%迭代求解

 for i=1:n
%   Q = -Qq*sin(i*dt);
    Q_eff = Q + M*(c0*d1 + c2*v1 + c3*a1) + C*(c1*d1 + c4*v1 + c5*a1);
    
    %d2 = K_eff\Q_eff;
    d2 = P*(L1'\(D1\(L1\(P'*Q_eff))));                       
    a2 = c0*(d2 - d1) - c2*v1 - c3*a1;
    v2 = v1 + c6*a1 + c7*a2;
    
    %xd(:,i+1) = d2;
    
    %t(i+1) = t(i) + dt;
    d1 = d2; v1 = v2; a1 = a2;
    
    ddx(i) = d1(1921);
    ddz(i) = d1(1922);
    ttt(i) = i*dt;
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % 输出为vtk格式
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     if (rem(i,1000)==0)%每100次迭代完成就输出一次，输出迭代时间域的末值
%         sdata.DIS(:,1) = d2;
%         sdata.V(:,1) = v2;
%         GetStress(L);
%         timestap = i;
%         vtkwrite('new',timestap);
%     end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % 输出为vtk格式
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
 end
figure;
plot(ttt,ddx,'r'),hold on;
plot(ttt,ddz,'k'),hold on;
title('位移时程响应曲线')
xlabel('t')
ylabel('x')   

TIM2 = clock;
time(1) = etime(TIM2, TIM1);

end