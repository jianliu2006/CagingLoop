function [ cagingPath ] = generateCagingGrasp( D,OUTPUTgrid,gridCOx,gridCOy,gridCOz,sourcePointId,index,...
                                               grid_on, dismap, gb_saddlePointId)
%   Compute the caging loop
%   Input:
%						gb_saddlePointId: the index of the saddle point
%   Output: 
%						cagingPath: caging loop
tree = nn_prepare(grid_on);
k=9;
cagePointId = getCagePoint(dismap,grid_on,tree,gb_saddlePointId,k);
path1 = compute_shortestpath(D,OUTPUTgrid,gridCOx,gridCOy,gridCOz ,cagePointId(1,1), sourcePointId,index);
path2 = compute_shortestpath(D,OUTPUTgrid,gridCOx,gridCOy,gridCOz ,cagePointId(2,1), sourcePointId,index);

path1 = [grid_on(gb_saddlePointId,:);path1];
path2 = [grid_on(gb_saddlePointId,:);path2];

path2 = flipud(path2);path2 = path2(2:size(path2,1),:);
finalPath = [path1;path2];
%%Smoothing the finalPath based on Laplacian
frontPath = [0 0 0;finalPath(1:(size(finalPath,1)-1),:)];
backPath = [finalPath(2:size(finalPath,1),:);0 0 0];
finalPath1 = 0.2*finalPath+0.4*frontPath+0.4*backPath;
startP = finalPath(2,:)*0.4+finalPath(size(finalPath,1),:)*0.4+finalPath(1,:)*0.2;
finalPath = [startP;finalPath1(2:(size(finalPath1,1)-1),:);startP];
cagingPath.finalPath = finalPath;

cagingPath.path1 = finalPath(1:size(path1,1),:);
cagingPath.path2 = finalPath(size(path1,1):size(finalPath,1),:);
end

function [cagePointId] = getCagePoint(dismap,grid_on,tree,pointId,k)
    cagePointId = [];
    %% Compute local coordinate system to get a plane
    [index1,dis1] =  nn_search(grid_on,tree,grid_on(pointId,:),15);
    X=grid_on(index1(1,2:15),:)-repmat(grid_on(pointId,:),[14,1]);
    [feature_vec,X1,feature_val,Psi]=pca(X',3);
    n_direct = feature_vec(:,3);u_direct = feature_vec(:,1);v_direct=feature_vec(:,2);
    %% Transform the global coordinate to the local coordinate
    [index,dis] =  nn_search(grid_on,tree,grid_on(pointId,:),k);
    adj_point = [];
    for j=2:k
        B = grid_on(index(1,j),:)-grid_on(pointId,:);B=B';
        A = [u_direct v_direct n_direct];
        X_lin = linsolve(A,B);
        %orig_X_lin = A*X_lin+grid_on(pointId,:)';
        adj_point = [adj_point;X_lin(1:2,1)'];
    end
    %% Get the sequence adj point
    bou = boundary(adj_point,0);
    adj_dismap_seq = dismap(index(1,bou+1)');
    %% Find the direction with lowtest increase speed of distance map
    judge = adj_dismap_seq-repmat(dismap(pointId),[size(adj_dismap_seq,1),1]);
    bou_max_negative = find(judge == min(judge));
    % First extend contact point
    cagePointId = [cagePointId bou(bou_max_negative(1,1))+1];
    
    adj_max = adj_point(bou(bou_max_negative(1,1)),:);
    adj_max_normal = adj_max/norm(adj_max);
    
    bou_negative = find(judge<0);
    id_bou_negative = bou(bou_negative);
    
    adj_negative = adj_point(id_bou_negative,:);
    adj_negative_len = sqrt(sum(adj_negative.*adj_negative,2));
    adj_negative_len = repmat(adj_negative_len ,[1,size(adj_negative,2)]);
    adj_negative_normal = adj_negative./adj_negative_len;
    
    adj_max_normal = repmat(adj_max_normal,[size(adj_negative_normal,1),1]);
    prod = sum(adj_negative_normal.*adj_max_normal,2);
    prod = find(prod == min(prod));
    % Second extend contac point
    cagePointId = [cagePointId id_bou_negative(prod(1,1))+1];
    cagePointId = index(1,cagePointId)';
end
