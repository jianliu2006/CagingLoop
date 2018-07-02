Code instruction for Caging Loops in Shape Embedding Space: Theory and Computation.
We only approvide the main interfaces of our approach.

===============================
        Main Interfaces
===============================
1. run pointCloudVoxelizationByRBF.m to discretize the grasping space based on RBF
2. run DistanceMapByFastMarching.m to compute the distance field rooted at source point
3. run detectSaddlePoint.m to detect the saddle points
4. run generateCagingGrasp.m to generate a caging loop connecting source point and saddle point

===============================
Note:
1. The source code is tested on Matlab R2015B
2. Before running the interfaces, you need to unzip some files and add Full Path of the tool folder in Matlab.
===============================
