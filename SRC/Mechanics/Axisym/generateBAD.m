<<<<<<< HEAD
function [B,DET,rphy]=generateBAD (psi,eta,XY)  %��֮����B������ʽ,��������ҽ�Ҫ��Χ���޸���
=======
function [B,DET,rphy]=generateBAD (LX,LY,XY)  %��֮����B������ʽ
>>>>>>> 增添了三种质量矩阵生成方式
%rphy�����뵱ǰ��psi��eta����뾶  
B = zeros(4,16);
global sdata;
nnode = sdata.NNODE; 
<<<<<<< HEAD
sumrow = sum(XY,2); %�н飬Ϊ�˼���ƽ���뾶����һ��,���淢�ֺ�abaqus������ܶ࣬�����Ҵ��㲻��ƽ���뾶������
rave = sumrow(1,1)/8;  %��ԳƵ�Ԫ�Ƚ����⣬��Ҫ��Ԫ��ƽ���뾶
	N(1) = 0.25*(1 - psi)*(1 - eta);
	N(2) = 0.25*(1 + psi)*(1 - eta);
	N(3) = 0.25*(1 + psi)*(1 + eta);
	N(4) = 0.25*(1 - psi)*(1 + eta);
	N(5) = 0.5*(1 - psi*psi)*(1 - eta); %��CQ4�������CQ8
	N(6) = 0.5*(1 + psi)*(1 - eta*eta);
	N(7) = 0.5*(1 - psi*psi)*(1 + eta);
	N(8) = 0.5*(1 - psi)*(1 - eta*eta);
	N(1) = N(1) - 0.5*N(5) - 0.5*N(8);%��ԭ����4���κ�����������
	N(2) = N(2) - 0.5*N(6) - 0.5*N(5);
	N(3) = N(3) - 0.5*N(7) - 0.5*N(6);
	N(4) = N(4) - 0.5*N(8) - 0.5*N(7);  %�����Nû�����⣬һ���ֳ��
    
    N11 = 0.25*(1-eta)*(2*psi+eta);%1��4�Žڵ��ǣ����½ǿ�ʼ��ʱ��ѭ����5-8���Ǵӵױߵ��е㿪ʼ��ʱ��ѭ��
	N12 = 0.25*(1-eta)*(2*psi-eta);
	N13 = 0.25*(1+eta)*(2*psi+eta);
	N14 = 0.25*(1+eta)*(2*psi-eta);

	N15 = -psi*(1-eta);  %����������4���ڵ��psi��
	N16 = 0.5*(1 - eta*eta);
	N17 = -psi*(1+eta);
	N18 = -0.5*(1 - eta*eta);

	N21 = 0.25*(1-psi)*(psi+2*eta);
	N22 = 0.25*(1+psi)*(-psi+2*eta);
	N23 = 0.25*(1+psi)*(psi+2*eta);
	N24 = 0.25*(1-psi)*(-psi+2*eta);

	N25 = -0.5*(1-psi*psi);  %����������4���ڵ��eta��
	N26 = -eta*(1+psi);
	N27 = 0.5*(1-psi*psi); 
	N28 = -eta*(1-psi);
%���ɵȲ�Ԫ���κ����ĵ���
Bpsi = [N11, N12, N13, N14, N15, N16, N17, N18;...
        N21, N22, N23, N24, N25, N26, N27, N28];
J = Bpsi*XY'; %��������ʽ
DET=det(J);
Bphy = inv(J)* Bpsi ;  %��������������κ����ĵ���  %�������û����
rphy = [N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)]*XY(1,:)';
=======
N = sdata.N(:,LX,LY);
Bpsi = sdata.BN(:,:,LX,LY);
J = Bpsi*XY'; %��������ʽ
DET=det(J);
Bphy = inv(J)* Bpsi ;  %��������������κ����ĵ���  %�������û����
rphy = N'*XY(1,:)';
>>>>>>> 增添了三种质量矩阵生成方式

for i=1:nnode %����B
B(1, i*2-1) = Bphy(1, i);
B(2, i*2) = Bphy(2, i);  %������
B(3, i*2-1:i*2) = [Bphy(2, i),Bphy(1, i)];
B(4, i*2-1:i*2) = [N(i)/rphy, 0];
end  
end