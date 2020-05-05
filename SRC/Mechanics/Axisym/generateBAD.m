function [B,DET,rphy]=generateBAD (psi,eta,XY)  %总之生成B和行列式,这个程序我将要大范围的修改了
%rphy是输入当前的psi和eta输出半径  
B = zeros(4,16);
global sdata;
nnode = sdata.NNODE; 
sumrow = sum(XY,2); %中介，为了计算平均半径利用一下,后面发现和abaqus结果相差很多，所以我打算不用平均半径计算了
rave = sumrow(1,1)/8;  %轴对称单元比较特殊，需要单元的平均半径
	N(1) = 0.25*(1 - psi)*(1 - eta);
	N(2) = 0.25*(1 + psi)*(1 - eta);
	N(3) = 0.25*(1 + psi)*(1 + eta);
	N(4) = 0.25*(1 - psi)*(1 + eta);
	N(5) = 0.5*(1 - psi*psi)*(1 - eta); %将CQ4增添成了CQ8
	N(6) = 0.5*(1 + psi)*(1 - eta*eta);
	N(7) = 0.5*(1 - psi*psi)*(1 + eta);
	N(8) = 0.5*(1 - psi)*(1 - eta*eta);
	N(1) = N(1) - 0.5*N(5) - 0.5*N(8);%对原来的4个形函数进行修正
	N(2) = N(2) - 0.5*N(6) - 0.5*N(5);
	N(3) = N(3) - 0.5*N(7) - 0.5*N(6);
	N(4) = N(4) - 0.5*N(8) - 0.5*N(7);  %检查了N没有问题，一部分抽查
    
    N11 = 0.25*(1-eta)*(2*psi+eta);%1到4号节点是，左下角开始逆时针循环，5-8号是从底边的中点开始逆时针循环
	N12 = 0.25*(1-eta)*(2*psi-eta);
	N13 = 0.25*(1+eta)*(2*psi+eta);
	N14 = 0.25*(1+eta)*(2*psi-eta);

	N15 = -psi*(1-eta);  %增添了其余4个节点对psi求导
	N16 = 0.5*(1 - eta*eta);
	N17 = -psi*(1+eta);
	N18 = -0.5*(1 - eta*eta);

	N21 = 0.25*(1-psi)*(psi+2*eta);
	N22 = 0.25*(1+psi)*(-psi+2*eta);
	N23 = 0.25*(1+psi)*(psi+2*eta);
	N24 = 0.25*(1-psi)*(-psi+2*eta);

	N25 = -0.5*(1-psi*psi);  %增添了其余4个节点对eta求导
	N26 = -eta*(1+psi);
	N27 = 0.5*(1-psi*psi); 
	N28 = -eta*(1-psi);
%生成等参元的形函数的导数
Bpsi = [N11, N12, N13, N14, N15, N16, N17, N18;...
        N21, N22, N23, N24, N25, N26, N27, N28];
J = Bpsi*XY'; %生成行列式
DET=det(J);
Bphy = inv(J)* Bpsi ;  %生成物理坐标的形函数的导数  %检查了下没问题
rphy = [N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)]*XY(1,:)';

for i=1:nnode %生成B
B(1, i*2-1) = Bphy(1, i);
B(2, i*2) = Bphy(2, i);  %改了下
B(3, i*2-1:i*2) = [Bphy(2, i),Bphy(1, i)];
B(4, i*2-1:i*2) = [N(i)/rphy, 0];
end  
end