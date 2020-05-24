function [beta] = Opt_beta()


%无阻尼情况下优化
%在sdata中添加了所需优化频率对应的阶数N_w


% Get global data
global sdata;


K = sdata.STIFF; M = sdata.MASS; %C = sdata.C; 
N_w = 20; %所需优化频率对应的阶数N_w
%beta = sdata.beta;



[w_m,w] = SSI(q,K,M,N_w);%子空间迭代法求解系统最高阶频率w_m和N_w阶频率w

Tmin = 2*pi/w_m;

%优化beta,dt不宜太小，最好取Tmin/10
if(dt<Tmin/10)
    c = 0;%阻尼比为0
    wd = w*(1 - c^2)^0.5;
    ta = (tan(wd*dt))^2;
    B = 1 + 2*c/(w*dt);
    x0 = (2/(B^2))*(1 + B*ta - (1 + 2*B*ta - B*B*ta)^0.5)/(1 + ta);
    beta = 1/x0 - 1/((w*dt)^2) - c/(w*dt);
    %sdata.beta = beta;
else
    printf('步长偏大，应满足dt<=%0.4e\n',Tmin/10);
end
    

