function [B,DET,rphy]=generateBAD (LX,LY,XY)  %��֮����B������ʽ
%rphy�����뵱ǰ��psi��eta����뾶  
B = zeros(4,16);
global sdata;
nnode = sdata.NNODE; 
N = sdata.N(:,LX,LY);
Bpsi = sdata.BN(:,:,LX,LY);
J = Bpsi*XY'; %��������ʽ
DET=det(J);
Bphy = inv(J)* Bpsi ;  %��������������κ����ĵ���  %�������û����
rphy = N'*XY(1,:)';

for i=1:nnode %����B
B(1, i*2-1) = Bphy(1, i);
B(2, i*2) = Bphy(2, i);  %������
B(3, i*2-1:i*2) = [Bphy(2, i),Bphy(1, i)];
B(4, i*2-1:i*2) = [N(i)/rphy, 0];
end  
end