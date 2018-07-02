function [ saddle ] = detectSaddlePoint( dismap, grid_on,sourcePointId, grid_on_normals )
%   Detecting saddle points
%   Input:
%           disamp: distance map of grid on the surface
%           grid_on: grids on the surface
%   Output：
%           the ids of saddle points
tree = nn_prepare(grid_on);
k=9;
saddle=[];
for i=1:size(grid_on,1)
    saddle_flag = calculateIterNum(dismap,grid_on,tree,i,k);
    if saddle_flag >= 4
        saddle=[saddle;i];
    end
end
%% Filtering illegal saddle point based on geometry property
score = [];
source_p = grid_on(sourcePointId,:);
source_normal = grid_on_normals(sourcePointId,:);
for i=1:size(saddle,1)
     p2 = grid_on(saddle(i,1),:);
     normal2 = grid_on_normals(saddle(i,1),:);
     score_temp = diversityEval(source_p,p2,source_normal,normal2);
     score = [score;score_temp];
 end
%% 取前10%个打分值高的的鞍点
[opt_Y,opt_I] = sort(score,'descend');
saddle = saddle(opt_I(1:round(size(score)/30)));
%% Test: display all the saddle points
ps = pointCloud(grid_on);
ps.plot('Color','r','MarkerSize', 8);
axis off;axis equal;movegui('northeast');view3d rot;
set(gcf,'color','white')
title('Detection of saddle points')
hold on
scatter3(grid_on(sourcePointId,1),grid_on(sourcePointId,2),grid_on(sourcePointId,3),83,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on
scatter3(grid_on(saddle,1),grid_on(saddle,2),grid_on(saddle,3),53,'MarkerEdgeColor','g','MarkerFaceColor','g');
end

function [ score ] = diversityEval(p1,p2,normal1,normal2)
    d = p2-p1;
    d=d/norm(d);
    angle1 = acos(dot(-1.0*normal1,d));
    angle2 = acos(dot(normal2,d));
    
    score = exp(-0.5*pow2(max(angle1,angle2),4));
end

function [saddle_flag] = calculateIterNum(dismap,grid_on,tree,pointId,k)
    saddle_flag = 0 ;
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
    %dlmwrite('adj_point.txt',adj_point);
    %% Get the sequence adj point
    bou = boundary(adj_point,0);
    [in on] = inpolygon(0,0,adj_point(bou,1),adj_point(bou,2));
    if in == 0 || on == 1
        saddle_flag = -1;
        return ;
    end
    bd2 = adj_point(bou,:);% 2d boundary point

    adj_dismap_seq = dismap(index(1,bou+1)');
    judge = adj_dismap_seq-repmat(dismap(pointId),[size(adj_dismap_seq,1),1]);
    
    % filter illegal boundary points;
    bou = filterBoundary(bd2,bou,judge);
  
    %% Judge saddle
    adj_dismap_seq = dismap(index(1,bou+1)');
    judge = adj_dismap_seq-repmat(dismap(pointId),[size(adj_dismap_seq,1),1]);
    judge(find(judge>0),1)=1;judge(find(judge<0))=0;
    t=1;
    while t<size(judge,1)
        if judge(t)~=judge(t+1)
            saddle_flag = saddle_flag+1;
        end
           t=t+1;
    end
end

function [newBou3] = filterBoundary(bd2,bou,judge)
    dis_toler = judge./sqrt(sum(bd2.*bd2,2));
    start_flag = find(abs(dis_toler) == max(abs(dis_toler)));
    
    bou = [bou(start_flag:size(bou,1),1);bou(2:start_flag,1)];
    bd2 = [bd2(start_flag:size(bd2,1),:);bd2(2:start_flag,:)];
    judge = [judge(start_flag:size(judge,1));judge(2:start_flag)];
    
    dis_toler = judge./sqrt(sum(bd2.*bd2,2));
    newBou1 = find(abs(dis_toler)<=50);% 距离比参数，比值小说明邻域点与源点距离差距小，这样的邻域点参考价值小
    
    a = bd2(1:(size(bd2)-1),:);
    b = bd2(2:size(bd2),:);
    p = sum(a.*b,2)./(sqrt(sum(a.*a,2)).*sqrt(sum(b.*b,2)));
    p = acos(p)/pi*180;
    newBou2 = find(p<=15)+1;% angle parameter
    newBou3 = bou;
    newBou3([newBou1;newBou2]) = []; 
end
