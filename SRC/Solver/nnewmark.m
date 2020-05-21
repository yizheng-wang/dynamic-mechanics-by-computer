function nnewmark(L,~)
beta = 5/9;
global sdata;
NN = sdata.NEQ;
K = sdata.STIFF; M = sdata.MASS;
CHNOD = sdata.CHNOD;%变化载荷用（包含载荷节点和Y方向坐标，见READFILE）
CHNOD = sortrows(CHNOD,2);
C = sparse(NN, NN);
misesnum = 0;%输出最大米塞斯应力的参数


move = 1;%激波到达单元数
v = 3000;%激波速度 3000mm/ms;
underpower = (CHNOD(move*2+1,2)-CHNOD(move*2-1,2))*100*2*pi*10/6;%均匀力时两边节点力;%上一单元贡献的力


% C = C+diags(diag(C)+10);
Qq = sdata.R(:,L); %变载荷时做为初始载荷
%固定载荷 也可改为在循环中随时间变化的载荷
Qq = 100*2*pi*10*Qq/4.5;%100Mpa的均匀内压
Q = Qq*0;%初始取为零

TIM1 = clock;%计时

%beta = 0.25;%newmark参数（相同参数下此方法和newmark方法计算结果完全相同）
gamma = 5/6;

nt = 1;% 每次迭代的时间步数
dt = 1e-5;
n = 1000;%总时间步数，也可改为总时间控制

q1 = (dt*dt*(1-2*beta)/2);%提前计算系数，减少计算量
q2 = gamma*dt;
q3 = beta*dt*dt;

xd = zeros(NN,nt, 'double');%迭代过程中使用的位移，速度，加速度，存储矩阵
xv =zeros(NN,nt, 'double');
xa = zeros(NN,nt, 'double');

[C1,K1,M1,C2,K2,M2] = split(C,K,M);% 矩阵分解为对角阵和剩余部分之和的形式，M阵为集中质量阵不需要分解
MM1 = diag(M1+q2*C1+q3*K1);
xa(:,nt) = M\(Q-K*xd(:,1));%计算初始加速度

for II = 1:n
    
    
    dd0 = xd+dt*xv+q1*xa;
    vv0 = xv+(dt-q2)*xa;
%     Q = Qq*sin(dt*II);
    error11 =10;
    while (1e2>error11)&&(error11>1e-8)
        
        ddn = xd;%记录下来用于判断收敛
        
        Q1 = Q+M2*xa+C2*xv+K2*xd-K1*dd0-C1*vv0;
        
        xa = Q1./MM1;
        xd = dd0+q3*xa;
        xv = vv0+q2*xa;
        
        
        error11 = abs(norm(ddn,1)-norm(xd,1));
        if (error11 ~= 0)
            error11 = error11/abs(norm(xd,1));
        end
        
    end
    
    if (error11>1)||isnan(error11)
        error(' *** ERROR *** boom');
    end
    
    
    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    %更新Q
    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    if (move<=60)
        dispower = v*dt*II-15;
        if (dispower>CHNOD(move*2+1,2))
            DD = CHNOD(move*2-1,1);
            
            if (DD > 0)
                Q(DD) = Qq(DD);%激波经过单元后更新为均匀应力
            end
            DD = CHNOD(move*2,1);
            if (DD > 0)
                Q(DD) = Qq(DD);
            end
            DD = CHNOD(move*2+1,1);
            if (DD > 0)
                Q(DD) = Qq(DD);
            end
            move = move+1;
            
        end
        if (move<=60)
            xxueq = dispower-CHNOD(move*2-1,2);
            xxeq = -1+xxueq/(CHNOD(move*2,2)-CHNOD(move*2-1,2));%激波波头在单元内的等参位置
            W1 = xxeq^3/6-xxeq^2/4+5/12;
            W2 = xxeq^3/6+xxeq^2/4-1/12;
            W3 = xxeq-xxeq^3/3+2/3;
            power = xxueq*100*2*pi*10;%本单元上合力 *2piR
            power1 = power*W1/(W1+W2+W3);%各节点力
            power2 = power*W3/(W1+W2+W3);
            power3 = power*W2/(W1+W2+W3);
            DD = CHNOD(move*2-1,1);
            if (DD > 0)
                Q(DD) = power1+Qq(DD)-underpower;%更新各个节点
            end
            DD = CHNOD(move*2,1);
            if (DD > 0)
                Q(DD) = power2;
            end
            DD = CHNOD(move*2+1,1);
            if (DD > 0)
                Q(DD) = power3;
            end
        end
    end
    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    %更新Q
    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    
    
    ddx(II) = xd(1921);
    ddz(II) = xd(1922);
    ttt(II) = II*dt;
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % 输出为vtk格式
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rem(II,100)==0)%每10次迭代计算一次应力
        sdata.DIS(:,1) = xd;
        sdata.V(:,1) = xv;
        GetStress(L);
        misesnum = misesnum+1;
        Mises(1,misesnum) = sdata.maxMises(1);
        Mises(2,misesnum) = sdata.maxMises(2);
        Mises(3,misesnum) = II*dt;
        if (rem(II,1000)==0)%每100次迭代完成就输出一次，输出迭代时间域的末值
            timestap = II
            vtkwrite('nnew',timestap);
        end
    end
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
figure;
plot(Mises(3,:),Mises(1,:),'r'),hold on;
title('最大米塞斯应力')
xlabel('t')
ylabel('x')
[maxMisesall,index] = max(Mises(1,:));
maxMisesall
mixtime = Mises(2,index)
TIM2 = clock;
time(1) = etime(TIM2, TIM1)

end



function [C1,K1,M1,C2,K2,M2] = split(C,K,M)
C1 = diag(diag(C));
K1 = diag(diag(K));
M1 = diag(diag(M));
C2 = C1-C;
K2 = K1-K;
M2 = M1-M;
end
