function [XY,T] = MergeClosePoints(XYin,cutoff)
%merge points that are too close
%XYin: size Nx2
%cutoff: threshold, specified by distance
sz=size(XYin);
if sz(1)==2
    XYin = XYin';
end
%% group points by cluster
Z = linkage(XYin,'ward');
% dendrogram(Z,'ColorThreshold',cutoff)
% W = inconsistent(Z);
% cutoff = 10;
T = cluster(Z,'cutoff',cutoff,'Criterion','distance');
figure(10001000);
gscatter(XYin(:,1),XYin(:,2),T);
title('cluster of points, ouput of MergeClosePoints');
%% merge points within a cluster
XY = zeros(max(T),2);
for ci = 1:max(T)
    xy = XYin(T==ci,:);
    xym = mean(xy,1);
    XY(ci,:)=xym;
end