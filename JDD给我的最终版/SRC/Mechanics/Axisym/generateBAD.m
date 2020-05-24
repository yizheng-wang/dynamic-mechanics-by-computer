function [B,DET,rphy]=generateBAD (LX,LY,XY)  %总之生成B和行列式
%rphy是输入当前的psi和eta输出半径  
B = zeros(4,16);
global sdata;
nnode = sdata.NNODE; 
N = sdata.N(:,LX,LY);
Bpsi = sdata.BN(:,:,LX,LY);
J = Bpsi*XY'; %生成行列式
DET=det(J);
Bphy = inv(J)* Bpsi ;  %生成物理坐标的形函数的导数  %检查了下没问题
rphy = N'*XY(1,:)';

for i=1:nnode %生成B
B(1, i*2-1) = Bphy(1, i);
B(2, i*2) = Bphy(2, i);  %改了下
B(3, i*2-1:i*2) = [Bphy(2, i),Bphy(1, i)];
B(4, i*2-1:i*2) = [N(i)/rphy, 0];
end  
end