function nnewmark(L,~)
beta = 5/9;
global sdata;
NN = sdata.NEQ;
K = sdata.STIFF; M = sdata.MASS;
CHNOD = sdata.CHNOD;%�仯�غ��ã������غɽڵ��Y�������꣬��READFILE��
CHNOD = sortrows(CHNOD,2);
C = sparse(NN, NN);
misesnum = 0;%����������˹Ӧ���Ĳ���


move = 1;%�������ﵥԪ��
v = 3000;%�����ٶ� 3000mm/ms;
underpower = (CHNOD(move*2+1,2)-CHNOD(move*2-1,2))*100*2*pi*10/6;%������ʱ���߽ڵ���;%��һ��Ԫ���׵���


% C = C+diags(diag(C)+10);
Qq = sdata.R(:,L); %���غ�ʱ��Ϊ��ʼ�غ�
%�̶��غ� Ҳ�ɸ�Ϊ��ѭ������ʱ��仯���غ�
Qq = 100*2*pi*10*Qq/4.5;%100Mpa�ľ�����ѹ
Q = Qq*0;%��ʼȡΪ��

TIM1 = clock;%��ʱ

%beta = 0.25;%newmark��������ͬ�����´˷�����newmark������������ȫ��ͬ��
gamma = 5/6;

nt = 1;% ÿ�ε�����ʱ�䲽��
dt = 1e-5;
n = 1000;%��ʱ�䲽����Ҳ�ɸ�Ϊ��ʱ�����

q1 = (dt*dt*(1-2*beta)/2);%��ǰ����ϵ�������ټ�����
q2 = gamma*dt;
q3 = beta*dt*dt;

xd = zeros(NN,nt, 'double');%����������ʹ�õ�λ�ƣ��ٶȣ����ٶȣ��洢����
xv =zeros(NN,nt, 'double');
xa = zeros(NN,nt, 'double');

[C1,K1,M1,C2,K2,M2] = split(C,K,M);% ����ֽ�Ϊ�Խ����ʣ�ಿ��֮�͵���ʽ��M��Ϊ������������Ҫ�ֽ�
MM1 = diag(M1+q2*C1+q3*K1);
xa(:,nt) = M\(Q-K*xd(:,1));%�����ʼ���ٶ�

for II = 1:n
    
    
    dd0 = xd+dt*xv+q1*xa;
    vv0 = xv+(dt-q2)*xa;
%     Q = Qq*sin(dt*II);
    error11 =10;
    while (1e2>error11)&&(error11>1e-8)
        
        ddn = xd;%��¼���������ж�����
        
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
    
    
    %��������������������������������������������������������������
    %����Q
    %��������������������������������������������������������������
    if (move<=60)
        dispower = v*dt*II-15;
        if (dispower>CHNOD(move*2+1,2))
            DD = CHNOD(move*2-1,1);
            
            if (DD > 0)
                Q(DD) = Qq(DD);%����������Ԫ�����Ϊ����Ӧ��
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
            xxeq = -1+xxueq/(CHNOD(move*2,2)-CHNOD(move*2-1,2));%������ͷ�ڵ�Ԫ�ڵĵȲ�λ��
            W1 = xxeq^3/6-xxeq^2/4+5/12;
            W2 = xxeq^3/6+xxeq^2/4-1/12;
            W3 = xxeq-xxeq^3/3+2/3;
            power = xxueq*100*2*pi*10;%����Ԫ�Ϻ��� *2piR
            power1 = power*W1/(W1+W2+W3);%���ڵ���
            power2 = power*W3/(W1+W2+W3);
            power3 = power*W2/(W1+W2+W3);
            DD = CHNOD(move*2-1,1);
            if (DD > 0)
                Q(DD) = power1+Qq(DD)-underpower;%���¸����ڵ�
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
    %��������������������������������������������������������������
    %����Q
    %��������������������������������������������������������������
    
    
    ddx(II) = xd(1921);
    ddz(II) = xd(1922);
    ttt(II) = II*dt;
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % ���Ϊvtk��ʽ
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rem(II,100)==0)%ÿ10�ε�������һ��Ӧ��
        sdata.DIS(:,1) = xd;
        sdata.V(:,1) = xv;
        GetStress(L);
        misesnum = misesnum+1;
        Mises(1,misesnum) = sdata.maxMises(1);
        Mises(2,misesnum) = sdata.maxMises(2);
        Mises(3,misesnum) = II*dt;
        if (rem(II,1000)==0)%ÿ100�ε�����ɾ����һ�Σ��������ʱ�����ĩֵ
            timestap = II
            vtkwrite('nnew',timestap);
        end
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % ���Ϊvtk��ʽ
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end
figure;
plot(ttt,ddx,'r'),hold on;
plot(ttt,ddz,'k'),hold on;
title('λ��ʱ����Ӧ����')
xlabel('t')
ylabel('x')
figure;
plot(Mises(3,:),Mises(1,:),'r'),hold on;
title('�������˹Ӧ��')
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
