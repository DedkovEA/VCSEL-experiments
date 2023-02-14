close all;
image([10,100], [1,4], stab, "CDataMapping", "scaled");
colorbar;
set(gca,'YDir','normal', "XScale", "log")
figure;
image([10,100], [1,4], real(qp - qm), "CDataMapping", "scaled");
colorbar;
set(gca,'YDir','normal', "XScale", "log")
%figure;
%image(psi, "CDataMapping", "scaled");