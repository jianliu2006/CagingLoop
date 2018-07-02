function [ dismap,D ] = DistanceMapByFastMarching( OUTPUTgrid, grid_on, index, startPointID )
%   Computation of distance map based on fast marching
%		Input:
%           startPointID: index of the source point
%   Output:
%           dismap: the distance map of grids on the surface
%           D: the distance map of all grids
W = double(OUTPUTgrid);
for i1 = 1:size(index.outer_index,1)
    W(index.outer_index(i1,1),index.outer_index(i1,2),index.outer_index(i1,3)) = 1;
end

for i3 = 1:size(index.inner_index,1)
    W(index.inner_index(i3,1),index.inner_index(i3,2),index.inner_index(i3,3)) =0;
end

for i2 = 1:size(index.on_index,1)
    W(index.on_index(i2,1),index.on_index(i2,2),index.on_index(i2,3)) = 1;
end

sourcePointID = index.on_index(startPointID,:);
tic
process = 'Time of processing fasting marching'
[D,Y]= msfm(W,sourcePointID');
toc
tree =nn_prepare(grid_on);
for i=1:size(index.on_index,1)
    dismap(i,1) = D(index.on_index(i,1),index.on_index(i,2),index.on_index(i,3));
    if dismap(i,1) >= 10000
        [local_index,dis] = nn_search(grid_on,tree,grid_on(i,:),9);
        tt = [index.on_index(local_index,1),index.on_index(local_index,2),index.on_index(local_index,3)];
        tem =[];
        for j=1:9
            tem = [tem;D(tt(j,1),tt(j,2),tt(j,3))];
        end
        dismap(i,1) = min(tem);
    end
end
%% Test: display the disntance field
ds = pointCloud(grid_on);
ds.A.alpha = dismap;
ds.plot('Color','A.alpha','MarkerSize', 13);
axis off;axis equal;movegui('northeast');view3d rot;
set(gcf,'color','white')
title('Distance field')
colormap jet
hold on
scatter3(grid_on(startPointID,1),grid_on(startPointID,2),grid_on(startPointID,3),183,'MarkerEdgeColor','r','MarkerFaceColor','r');
end

