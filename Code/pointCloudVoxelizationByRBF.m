function [ OUTPUTgrid,gridCOx,gridCOy,gridCOz,index,grid_on,grids_inner,grids_outer,FV ] = pointCloudVoxelizationByRBF( ptCloud,ptnormals,voxel_xnum,voxel_ynum,voxel_znum )
%   Point cloud voxelization
%   Input:
%           ptCloud: original point cloud data (or vertices of triangle mesh) pointNum*3
%           ptnormals: normals of original point cloud (or vertex normal) pointNum*3
%   Output:
%           OUTPUTgrid: voxel_xnum*voxel_ynum*voxel_znum voxels
%                          1: grid on the surface
%                          0: grid inside the surface
%                          -1: grid outside the surface
x_max = max(ptCloud(:,1));x_min = min(ptCloud(:,1));
y_max = max(ptCloud(:,2));y_min = min(ptCloud(:,2));
z_max = max(ptCloud(:,3));z_min = min(ptCloud(:,3));

lx = (x_max-x_min)/(voxel_xnum-1);
ly = (y_max-y_min)/(voxel_ynum-1);
lz = (z_max-z_min)/(voxel_znum-1);

% offset = (lx+ly+lz)/3;%0.0005

%% Approximation based on RBF
num = size(ptCloud,1);
pos_constraint = ptCloud+ptnormals*0.0001;pos_v = ones(1,num);
neg_constraint = ptCloud+ptnormals*(-1.0*0.0001);neg_v = -1.0*ones(1,num);
normal_v = zeros(1,num);
v=[normal_v pos_v neg_v];
samples = [ptCloud' pos_constraint' neg_constraint'];
% save sampling points and their value to txt file format
tic
samplesAndValue=[samples;v];
fileID = fopen('samples.txt','w');
fprintf(fileID,'%f %f %f %f\n',samplesAndValue);
fclose(fileID);
% Create implict function based on fastRBF
rbfImportCommand = ['FastRBF.exe',' ','import',' ','-ascii',' ','-format=','"%x %y %z %v"',' ','samples.txt',' ','sample.aranz'];
system(rbfImportCommand);
process = 'the cost of creation of implitic function based on RBF'
tic
rbfFitCommand = ['FastRBF.exe',' ','fit',' ','-ascii',' ','-accuracy=1.0',' ','-rho=0.01',' ','sample.aranz',' ','samplerbf.rbf'];
system(rbfFitCommand);
toc
% Create implict function based on rbfinterp_v1.2
%tic
% op = rbfcreate(samples,v,'RBFFunction', 'linear', 'Stats', 'on');
% toc
% rbfcheck(op);
%% voxelization of point cloud
OUTPUTgrid = ones(voxel_xnum,voxel_ynum,voxel_znum);
OUTPUTgrid = -2.0*OUTPUTgrid;% default value -2.0
ADJgrid = {}; 

grid_on = [];
grids_outer = [];
grids_inner = [];

FV=[];
index = [];
on_index = [];
outer_index=[];
inner_index = [];

gridCOx = x_min:lx:x_max;
gridCOy = y_min:ly:y_max;
gridCOz = z_min:lz:z_max;

[x,y,z] = meshgrid(gridCOx,gridCOy,gridCOz);

x= reshape(x,[(voxel_xnum)*(voxel_xnum)*(voxel_xnum),1]);
y= reshape(y,[(voxel_ynum)*(voxel_ynum)*(voxel_ynum),1]);
z= reshape(z,[(voxel_znum)*(voxel_znum)*(voxel_znum),1]);
samplings=[x y z];

%% Evaluate rbf basded on FastRBF
tic
valids = samplings';
fileID = fopen('valids.txt','w');
fprintf(fileID,'%f %f %f\n',valids);
fclose(fileID);
rbfImportValidsCommand =['FastRBF.exe',' ','import',' ','-ascii',' ','-format=','"%x %y %z"',' ','valids.txt',' ','valids.aranz']; 
system(rbfImportValidsCommand)
%tic
process= 'the cost of evaluation of rbf value of points'
rbfPointEvalCommand = ['FastRBF.exe',' ','pointeval',' ','-ascii',' ','-accuracy=1.0',' ','valids.aranz',' ','samplerbf.rbf',' ','validsrbf.aranz'];
system(rbfPointEvalCommand)
%toc
rbfExportCommand = ['FastRBF.exe',' ','export',' ','-ascii',' ','-txt',' ','-format=','"%x %y %z %v"',' ','validsrbf.aranz',' ','validsResult.txt'];
system(rbfExportCommand)
fileID = fopen('validsResult.txt','r');
validsResult = fscanf(fileID,'%f %f %f %f',[4 inf]);
rbfv = validsResult(4,:)';
toc
% Delete generated file
delete valids.txt;delete valids.aranz;delete validsrbf.aranz;
delete validsResult.txt;

tree_boundary = nn_prepare(samplings);
[b_index,b_dis]=nn_search(samplings,tree_boundary,samplings,9);
rbf_adjs = rbfv(b_index(:,2:7));
[max_rbf,adjMax_index] = max(rbf_adjs,[],2);
[min_rbf,adjMin_index] = min(rbf_adjs,[],2);
ff = max_rbf.*min_rbf;
    grid_on = samplings(find(ff<=1.0),:);
% grid_on = samplings(find((rbfv*offset)<=(0.5*offset*0.25) & (rbfv*offset) >= (-0.5*offset*0.25)),:);
% grid_on = samplings(find((rbfv*offset)<=(0.8*offset) & (rbfv*offset) >= (-0.8*offset)),:);

x_max = max(grid_on(:,1));x_min = min(grid_on(:,1));
y_max = max(grid_on(:,2));y_min = min(grid_on(:,2));
z_max = max(grid_on(:,3));z_min = min(grid_on(:,3));

lx = (x_max-x_min)/(voxel_xnum-1);
ly = (y_max-y_min)/(voxel_ynum-1);
lz = (z_max-z_min)/(voxel_znum-1);

x_index = round((grid_on(:,1)-repmat(x_min,[size(grid_on,1),1]))/lx)+1;
y_index = round((grid_on(:,2)-repmat(y_min,[size(grid_on,1),1]))/ly)+1;
z_index = round((grid_on(:,3)-repmat(z_min,[size(grid_on,1),1]))/lz)+1;

on_index = [x_index y_index z_index];

for i=1:size(on_index,1)
    OUTPUTgrid(on_index(i,1),on_index(i,2),on_index(i,3))=1;
end

grid_on = [gridCOx(1,x_index(:,1))' gridCOy(1,y_index(:,1))' gridCOz(1,z_index(:,1))'];
%% compute inner and outer
index_inner_outer = [];
for j=1:voxel_znum
    [index_x,index_y]=find( OUTPUTgrid(:,:,j) == -2 );
    index_inner_outer = [index_inner_outer;index_x index_y repmat(j,[size(index_x,1),1])];
end
temp = [gridCOx(1,index_inner_outer(:,1))' gridCOy(1,index_inner_outer(:,2))' gridCOz(1,index_inner_outer(:,3))'];
% Evaluate rbf based on rbfinterp_v1.2
%rbfv = rbfinterp(temp',op);rbfv = rbfv';

% Evaluate rbf basded on FastRBF
tic
valids = temp';
fileID = fopen('valids.txt','w');
fprintf(fileID,'%f %f %f\n',valids);
fclose(fileID);
rbfImportValidsCommand =['FastRBF.exe',' ','import',' ','-ascii',' ','-format=','"%x %y %z"',' ','valids.txt',' ','valids.aranz']; 
system(rbfImportValidsCommand)
%tic
process= 'the cost of evaluation of rbf value of points'
rbfPointEvalCommand = ['FastRBF.exe',' ','pointeval',' ','-ascii',' ','-accuracy=1.0',' ','valids.aranz',' ','samplerbf.rbf',' ','validsrbf.aranz'];
system(rbfPointEvalCommand)
%toc
rbfExportCommand = ['FastRBF.exe',' ','export',' ','-ascii',' ','-txt',' ','-format=','"%x %y %z %v"',' ','validsrbf.aranz',' ','validsResult.txt'];
system(rbfExportCommand)
fileID = fopen('validsResult.txt','r');
validsResult = fscanf(fileID,'%f %f %f %f',[4 inf]);
rbfv = validsResult(4,:)';
toc
% Delete generate files
delete valids.txt;delete valids.aranz;delete validsrbf.aranz;
delete validsResult.txt;

grids_inner = temp(find(rbfv<0),:);inner_index = index_inner_outer(find(rbfv<0),:);
grids_outer = temp(find(rbfv>0),:);outer_index = index_inner_outer(find(rbfv>0),:);
%% Delete the boundary from grids_inner and add the boundary to grids_outer
% find the boundary index from inner_index
xmin_index = inner_index(find(inner_index(:,1)==1),:);
xmax_index = inner_index(find(inner_index(:,1)==voxel_xnum),:);

ymin_index = inner_index(find(inner_index(:,2)==1),:);
ymax_index = inner_index(find(inner_index(:,2)==voxel_ynum),:);

zmin_index = inner_index(find(inner_index(:,3)==1),:);
zmax_index = inner_index(find(inner_index(:,3)==voxel_znum),:);

% add boundary to grid_outer index&point
[tf_x1,loc_x1]=ismember(inner_index,xmax_index,'rows');
[tf_x2,loc_x2]=ismember(inner_index,xmin_index,'rows');

[tf_y1,loc_y1]=ismember(inner_index,ymax_index,'rows');
[tf_y2,loc_y2]=ismember(inner_index,ymin_index,'rows');

[tf_z1,loc_z1]=ismember(inner_index,zmax_index,'rows');
[tf_z2,loc_z2]=ismember(inner_index,zmin_index,'rows');

outer_index = [outer_index;inner_index(find(tf_x1==1),:);...
                inner_index(find(tf_y1==1),:);...
                inner_index(find(tf_z1==1),:);...
                inner_index(find(tf_x2==1),:);...
                inner_index(find(tf_y2==1),:)];
            %inner_index(find(tf_z2==1),:)

grids_outer = [grids_outer;grids_inner(find(tf_x1==1),:);...
                grids_inner(find(tf_y1==1),:);...
                grids_inner(find(tf_z1==1),:);...
                grids_inner(find(tf_x2==1),:);...
                grids_inner(find(tf_y2==1),:)];
             % grids_inner(find(tf_z2==1),:)
% delete boundary from grid_inner index&point
delete_index = [find(tf_x1==1);find(tf_y1==1);find(tf_z1==1);...
                find(tf_x2==1);find(tf_y2==1)];%;find(tf_z2==1)
inner_index(delete_index,:) = [];
grids_inner(delete_index,:) = [];

for i=1:size(outer_index,1)
    OUTPUTgrid(outer_index(i,1),outer_index(i,2),outer_index(i,3))=-1;
end

for i=1:size(inner_index,1)
    OUTPUTgrid(inner_index(i,1),inner_index(i,2),inner_index(i,3))=0;
end
index.on_index = on_index;
index.inner_index = inner_index;
index.outer_index =outer_index;

%% surface reconstruction basded FastRBF
[x,y,z] = meshgrid(gridCOx,gridCOy,gridCOz);
grids = [x(:)';y(:)';z(:)'];
x_min = min(grids(1,:));y_min = min(grids(2,:));z_min = min(grids(3,:));
x_max = max(grids(1,:));y_max = max(grids(2,:));z_max = max(grids(3,:));
process = 'the cost of interpolation grids'
tic
rbfGridEvalCommand = ['FastRBF.exe',' ','grideval',' ','-ascii',' ','-accuracy=1.0',' ','-size=',int2str(voxel_xnum),' '...
    '-min=',num2str(x_min),',',num2str(y_min),',',num2str(z_min),' '...
    '-max=',num2str(x_max),',',num2str(y_max),',',num2str(z_max),' '...
    ,'samplerbf.rbf',' ','gridfile.aranz'];
system(rbfGridEvalCommand);
toc
process = 'generating mesh and save mesh'
tic
rbfMcubesCommand = ['FastRBF.exe',' ','mcubes',' ','-ascii',' ','gridfile.aranz',' ','gridmesh.aranz'];
system(rbfMcubesCommand);
rbfExportMeshCommand = ['FastRBF.exe',' ','export',' ','-ascii',' ','gridmesh.aranz',' ','gridmesh.obj'];
system(rbfExportMeshCommand);
toc
delete gridfile.aranz;delete gridmesh.aranz
%delete samplerbf.rbf;
%% display the grid on and grid_inner and reconstructed surface
[verts,faces] = readOBJ('gridmesh.obj');
FV.faces = faces;
FV.verts =verts;

 p=pointCloud(grid_on);
pi=pointCloud(grids_inner);
% po = pointCloud(grids_outer);

pi.plot('Color','g','MarkerSize',8);hold on
p.plot('Color','r','MarkerSize',8);hold on
% po.plot('Color','b','MarkerSize',15)
  patch('vertices',verts,'faces',faces,'edgecolor','none',...
         'facecolor',[0 0 1],'FaceLighting', 'phong');hold off
  light
axis off;axis equal;movegui('northeast');view3d rot;
set(gcf,'color','white')
title('voxelization and reconstruction')
end

