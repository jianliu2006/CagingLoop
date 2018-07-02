function [ path ] = compute_shortestpath(D,OUTPUTgrid,gridCOx,gridCOy,gridCOz ,startPointId, sourcePointId,index)
%   Computation of shortest path from ending point to source point based on
%   fast marching
%   Input:
%           D: distance map of the whole voxels
%           startPointId: the index of ending point on the surface
%           sourcePointId: the index of source point onthe surface
%   Output:
%           path: shortest path from ending point to source point
%   Detailed explanation goes here
n=length(gridCOx);
path=[];
start_p = index.on_index(startPointId,:)';
source_p = index.on_index(sourcePointId,:)';
shortestLine = shortestpath(D,start_p,source_p);

x0 = gridCOx(1,1); y0 = gridCOy(1,1); z0 = gridCOz(1,1);
x_max = gridCOx(1,n); y_max = gridCOy(1,n); z_max = gridCOz(1,n);

r_x = (x_max-x0)/(n-1);r_y=(y_max-y0)/(n-1);r_z = (z_max-z0)/(n-1);
p_x = x0-r_x;p_y = y0-r_y;p_z = z0-r_z;

path_x = shortestLine(:,1)*r_x+repmat(p_x,[size(shortestLine,1),1]);
path_y = shortestLine(:,2)*r_y+repmat(p_y,[size(shortestLine,1),1]);
path_z = shortestLine(:,3)*r_z+repmat(p_z,[size(shortestLine,1),1]);

pathstart=[gridCOx(1,start_p(1,1)) gridCOy(1,start_p(2,1)) gridCOz(1,start_p(3,1))];
path = [path_x path_y path_z];
path = [pathstart;path];
end

