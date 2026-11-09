clc;
cfg.inputFile = 'E:\02 论文撰写\01 城市遥感\01 原始数据\山西临汾-42.8w-城市建筑群.txt';
cfg.outputFolderName = 'UrbanDEM_Linfen_Paper_Reimplementation';
% 基础栅格与自适应垂向分层参数
cfg.gridResolution = 0.50;          % m；水平栅格单元边长
cfg.Tmin = 0.15;                    % m；最小垂向分层间隔
cfg.lambda = 1.50;                  % MAD自适应调节系数
cfg.madCorrection = 1.4826;        % 正态分布一致性尺度修正系数
% 候选高程层的空间覆盖、离散度与稳定性评分参数
cfg.M = 4;                         % 将每个栅格细分为4×4个正方形子区域
cfg.layerWeights = [0.35,0.35,0.30]; % 点数比例、水平覆盖率、高程稳定性三项权重
cfg.sigma0 = 0.21;                 % m；高程离散度归一化尺度
cfg.Nmin = 3;                      % 每个候选层的最少有效点数
cfg.Cmin = 0.13;                   % 最小水平覆盖比例
cfg.penalty = 0.30;                % 对点数或覆盖率不足的候选层施加的惩罚值
% 邻域高程连续性校正参数
cfg.continuityWeights = [0.60,0.40];
cfg.Tc = 0.40;                     % m；邻域高程连续性尺度
cfg.neighbourhoodSize = 3;         % 采用3×3邻域
% 异常低值点剔除参数
cfg.Tb = 0.30;                     % m；底部噪声剔除阈值
% 多尺度形态学开运算与地形自适应阈值参数
cfg.radii = [1,2,4,8,12,18,26,38,54,72,96]; % 尺度值以栅格数表示；
cfg.orientations = [0,45,90,135];  % 线形结构元素的四个方向（角度）
cfg.T0 = 0.12;                     % m；动态阈值的基础值
cfg.alpha = 0.10;                  % 每增加一个尺度级次的阈值增量，而非按半径递增
cfg.beta = 0.30;                   % 局部坡度影响系数
cfg.gamma = 0.30;                  % 局部地形起伏影响系数
cfg.reliefWindow = 5;              % 5×5地形起伏计算窗口
% 地形先验、地面点恢复及DEM平滑参数
cfg.D0 = 1.50;                     % m；地形先验可信度的距离控制尺度
cfg.Tf = 0.25;                     % m；高置信地面点高程差阈值
cfg.Tfc = 0.50;                    % m；允许补偿恢复的高程差上限
cfg.rho0 = 0.70;                   % 邻域内可靠地面点比例阈值
cfg.groundNeighbourhoodRadius = 1.0; % m；地面点邻域距离，后续转换为栅格窗口
cfg.refineIterations = 3;         % 地面点补偿迭代次数上限
cfg.smoothingWeight = 0.80;       % 最终DEM平滑混合权重
cfg.smoothingWindow = 3;          % 3×3局部均值平滑窗口
cfg.referenceLabelFile = '';
cfg.maxScatterPoints = 100000;   % 单幅点云散点图最多显示的点数
cfg.saveMat = true;
cfg.saveFigures = true;
cfg.makeHoleSummary = true;
% ======================= 开始执行 =======================
if exist(cfg.inputFile,'file') ~= 2
    error(['Input point cloud not found:\n%s\n' ...
'Please verify the Windows E: drive and TXT filename.'],cfg.inputFile);
end
if exist('strel','file') ~= 2 || exist('imopen','file') ~= 2
    error('The Image Processing Toolbox is required (strel/imopen missing).');
end
checkConfiguration(cfg);
[inputDir,~,~] = fileparts(cfg.inputFile);
outDir = fullfile(inputDir,cfg.outputFolderName);
if exist(outDir,'dir') ~= 7, mkdir(outDir); end
diary(fullfile(outDir,'run_log.txt'));
diaryGuard = onCleanup(@() diary('off'));
fprintf('\n=========================================================\n');
fprintf('Urban DEM - manuscript-based MATLAB implementation\n');
fprintf('Input: %s\nOutput: %s\n',cfg.inputFile,outDir);
fprintf('=========================================================\n');
totalClock = tic;
% 第一步：读取原始点云并建立固定分辨率的水平栅格。
ticRead = tic;
raw = readPointCloud(cfg.inputFile);
xyz = raw(:,1:3);
N = size(xyz,1);
[xCoordinates,yCoordinates,cellID,nRows,nCols,xGrid,yGrid] = ...
    makeGrid(xyz,cfg.gridResolution);
nCells = nRows*nCols;
occupiedCount = accumarray(cellID,1,[nCells,1],@sum,0);
occupiedGrid = reshape(occupiedCount>0,nRows,nCols);
fprintf('[Input] %d points; %d columns; grid = %d columns x %d rows; cell = %.2f m\n', ...
    N,size(raw,2),nCols,nRows,cfg.gridResolution);
fprintf('[Input] Occupied raster cells: %d; read + grid %.2f s\n', ...
    nnz(occupiedGrid),toc(ticRead));
% 第二步：逐栅格执行高程分层和底层选择，
% 再利用邻域连续性校正并检测异常低值点。
ticLayer = tic;
[noisePoint,initialBottom,correctedBottom,cellMAD,stage1] = ...
    stageBottomIdentification(xyz,cellID,nRows,nCols,xCoordinates,yCoordinates,cfg);
validPoint = ~noisePoint;
fprintf('[Stage 1] Rejected anomalously low points: %d (%.2f%%); %.2f s\n', ...
    nnz(noisePoint),100*mean(noisePoint),toc(ticLayer));
% 第三步：计算多方向、多尺度形态学结果，
% 引入地形阈值并逐级更新候选地面栅格。
ticPMF = tic;
[Z,validCellAfterNoise] = pointMinimumGrid(xyz(:,3),cellID,validPoint,nRows,nCols);
Zfilled = interpolateGridSurface(Z,validCellAfterNoise,xGrid,yGrid,'natural');
[candidateCell,pmfResidual,localSlope,localRelief,thresholdFinal,levelStat] = ...
    stageTerrainAwarePMF(Zfilled,validCellAfterNoise,cfg);
fprintf('[Stage 2] Candidate ground cells: %d; %.2f s\n', ...
    nnz(candidateCell),toc(ticPMF));
% 第四步：构建地形先验、计算可信度并精筛候选地面点。
ticRefine = tic;
[candidatePoint,groundPoint,terrainPrior,confidenceGrid,pointHeightResidual,rhoGrid] = ...
    stageTerrainPriorAndRefinement(xyz,cellID,validPoint, ...
                                  candidateCell,Zfilled,xGrid,yGrid,cfg);
fprintf('[Stage 3] Ground: %d (%.2f%%); non-ground: %d (%.2f%%); %.2f s\n', ...
    nnz(groundPoint),100*mean(groundPoint),nnz(~groundPoint), ...
    100*mean(~groundPoint),toc(ticRefine));
% 第五步：依据最终地面点构建DEM，执行空值补全与平滑。
ticDem = tic;
[DEM,groundObserved,DEMunsmoothed] = stageDEM( ...
    xyz(:,3),cellID,groundPoint,terrainPrior,xGrid,yGrid,cfg);
if cfg.makeHoleSummary
    holeSummary = summarizeHoles(groundObserved);
else
    holeSummary = struct('count',NaN,'small',NaN,'elongated',NaN,'large',NaN);
end
fprintf('[Stage 4] DEM: %d x %d; internal gaps = %g (small %g, elongated %g, large %g); %.2f s\n', ...
    nCols,nRows,holeSummary.count,holeSummary.small,holeSummary.elongated, ...
    holeSummary.large,toc(ticDem));
% 第六步：导出数值结果；存在独立参考标签时计算分类指标。
classLabel = uint8(groundPoint);
outputTable = array2table(raw);
colNames = cell(1,size(raw,2));
colNames(1:3) = {'X','Y','Z'};
for j = 4:size(raw,2), colNames{j} = sprintf('Attribute%d',j); end
outputTable.Properties.VariableNames = colNames;
outputTable.Class = classLabel;            % 1表示地面点，0表示非地面点
writetable(outputTable,fullfile(outDir,'classified_points.csv'));
writematrix(xyz(groundPoint,:),fullfile(outDir,'ground_points_XYZ.txt'),'Delimiter','tab');
writematrix(xyz(~groundPoint,:),fullfile(outDir,'non_ground_points_XYZ.txt'),'Delimiter','tab');
writematrix(xyz(noisePoint,:),fullfile(outDir,'bottom_noise_points_XYZ.txt'),'Delimiter','tab');
[XX,YY] = meshgrid(xGrid,yGrid);
writematrix([XX(:),YY(:),DEM(:)],fullfile(outDir,'DEM_XYZ.csv'));
writematrix(DEM,fullfile(outDir,'DEM_matrix.csv'));
writematrix(double(candidateCell),fullfile(outDir,'candidate_ground_cells.csv'));
writematrix(pmfResidual,fullfile(outDir,'PMF_residual.csv'));
writematrix(terrainPrior,fullfile(outDir,'terrain_prior.csv'));
writematrix(confidenceGrid,fullfile(outDir,'terrain_confidence.csv'));
writematrix(correctedBottom,fullfile(outDir,'corrected_bottom_elevation.csv'));
accuracy = struct('available',false);
if ~isempty(cfg.referenceLabelFile)
    refLabels = readReferenceLabels(cfg.referenceLabelFile,N);
    accuracy = evaluateClassification(groundPoint,refLabels);
    fprintf('[Reference] Type I %.4f%%; Type II %.4f%%; Total %.4f%%; Kappa %.4f%%\n', ...
        accuracy.typeI,accuracy.typeII,accuracy.total,100*accuracy.kappa);
    metricsTable = struct2table(accuracy);
    writetable(metricsTable,fullfile(outDir,'classification_metrics.csv'));
else
    fprintf('[Reference] No independent labels supplied: Type I/II, total error and Kappa NOT calculated.\n');
end
% 第七步：导出原始点云、残差、候选栅格与分类结果组成的诊断图。
if cfg.saveFigures
    try
        saveResultFigures(xyz,groundPoint,pmfResidual,candidateCell,DEM,xGrid,yGrid,cfg,outDir);
    catch figErr
        % ========== 修复warning报错：使用MException标准格式 ==========
        warning(figErr.identifier, 'Figures could not be saved: %s', figErr.message);
    end
end
elapsed = toc(totalClock);
summary = table(N,nCols,nRows,cfg.gridResolution,nnz(noisePoint), ...
    nnz(candidateCell),nnz(groundPoint),nnz(~groundPoint),elapsed, ...
'VariableNames',{'InputPoints','GridColumns','GridRows','GridResolution_m', ...
'BottomNoisePoints','CandidateGroundCells','GroundPoints','NonGroundPoints','ElapsedSeconds'});
writetable(summary,fullfile(outDir,'run_summary.csv'));
results = struct();
results.config = cfg;
results.outputDirectory = outDir;
results.inputPointCount = N;
results.gridSize = [nRows,nCols];
results.groundPointCount = nnz(groundPoint);
results.nonGroundPointCount = nnz(~groundPoint);
results.noisePointCount = nnz(noisePoint);
results.candidateCellCount = nnz(candidateCell);
results.elapsedSeconds = elapsed;
results.accuracy = accuracy;
results.holeSummary = holeSummary;
results.scaleHistory = levelStat;
results.groundMask = groundPoint;
results.bottomNoiseMask = noisePoint;
results.candidateCellMask = candidateCell;
results.initialBottom = initialBottom;
results.correctedBottom = correctedBottom;
results.cellMAD = cellMAD;
results.DEM = DEM;
results.DEMunsmoothed = DEMunsmoothed;
results.terrainPrior = terrainPrior;
results.terrainConfidence = confidenceGrid;
results.PMFresidual = pmfResidual;
results.localSlope = localSlope;
results.localRelief = localRelief;
results.lastThreshold = thresholdFinal;
results.pointHeightResidual = pointHeightResidual;
results.neighbourhoodGroundRatio = rhoGrid;
results.stage1 = stage1;
if cfg.saveMat
    save(fullfile(outDir,'UrbanDEM_Linfen_results.mat'),'results','-v7.3');
end
fprintf('\nDone. Total elapsed time %.2f seconds (%.2f minutes).\n',elapsed,elapsed/60);
fprintf('Products saved to: %s\n',outDir);
fprintf('NOTE: This implementation is based on the equations and parameters\n');
fprintf('of the manuscript; it is not the original V3.4.1 source code.\n');
 
% ======================== 辅助函数 ========================
function checkConfiguration(c)
assert(c.gridResolution>0 && c.M>=1 && c.M==round(c.M));
assert(c.neighbourhoodSize>=1 && mod(c.neighbourhoodSize,2)==1);
assert(c.reliefWindow>=1 && mod(c.reliefWindow,2)==1);
assert(c.smoothingWindow>=1 && mod(c.smoothingWindow,2)==1);
assert(numel(c.layerWeights)==3 && abs(sum(c.layerWeights)-1)<1e-10);
assert(numel(c.continuityWeights)==2 && abs(sum(c.continuityWeights)-1)<1e-10);
assert(c.Tc>0 && c.sigma0>0 && c.D0>0 && c.Tfc>c.Tf);
assert(all(c.radii>=1) && all(diff(c.radii)>0));
assert(c.refineIterations>=1);
end
function A = readPointCloud(fname)
% 保留所有数值列
try
    A = readmatrix(fname,'FileType','text');
catch
    A = readmatrix(fname);
end
if isempty(A) || size(A,2)<3
    error('TXT file must contain at least three numeric columns: X Y Z.');
end
if ~isnumeric(A), error('Input data could not be read as numbers.'); end
validRows = all(isfinite(A(:,1:3)),2);
if nnz(~validRows)>0
    warning('Discarding %d rows with invalid X/Y/Z (e.g. text headers).',nnz(~validRows));
    A = A(validRows,:);
end
if size(A,1)<3, error('Need at least three valid XYZ points.'); end
A = double(A);
end
function [gx,gy,ids,nRows,nCols,xg,yg] = makeGrid(xyz,res)
gx = min(xyz(:,1));
gy = min(xyz(:,2));
nCols = floor((max(xyz(:,1))-gx)/res)+1;
nRows = floor((max(xyz(:,2))-gy)/res)+1;
if nRows*nCols>3e7
    error('Grid exceeds 30 million cells. Check units / input coordinates.');
end
column = min(nCols,max(1,floor((xyz(:,1)-gx)/res)+1));
row = min(nRows,max(1,floor((xyz(:,2)-gy)/res)+1));
ids = sub2ind([nRows,nCols],row,column);
xg = gx + ((0:nCols-1)+0.5)*res;
yg = gy + ((0:nRows-1)+0.5)*res;
end
function [noise,initialBottom,finalBottom,madRaster,diag] = ...
    stageBottomIdentification(xyz,ids,nRows,nCols,x0,y0,c)
% 每个栅格单独处理高程分层，先按层内统计量筛选，
% 再通过邻域参考高程修正可靠底层。
N = size(xyz,1);
numCells = nRows*nCols;
[sortedCells,order] = sort(ids);
blockStarts = [1; find(diff(sortedCells)~=0)+1; N+1];
numGroups = numel(blockStarts)-1;
layerCache = cell(numCells,1);
selectedInitial = nan(numCells,1);
initialLower = nan(numCells,1);
allMad = nan(numCells,1);
validLayers = 0;
for g = 1:numGroups
    ind = order(blockStarts(g):blockStarts(g+1)-1);
    cellNumber = sortedCells(blockStarts(g));
    Z = xyz(ind,3);
    [zSorted,perm] = sort(Z,'ascend');
    indSorted = ind(perm);
    n = numel(indSorted);
    midZ = median(zSorted);
    madValue = c.madCorrection*median(abs(zSorted-midZ));
    allMad(cellNumber) = madValue;
    gap = max(c.Tmin,c.lambda*madValue);
    cuts = [0;find(diff(zSorted)>gap);n];
    [rr,cc] = ind2sub([nRows,nCols],cellNumber);
    cellStartX = x0+(cc-1)*c.gridResolution;
    cellStartY = y0+(rr-1)*c.gridResolution;
    K = numel(cuts)-1;
    support = zeros(K,1);
    coverage = zeros(K,1);
    layerMedian = zeros(K,1);
    layerMin = zeros(K,1);
    sigma = zeros(K,1);
    score = zeros(K,1);
    for k=1:K
        kk = cuts(k)+1:cuts(k+1);
        groupInd = indSorted(kk);
        zLocal = zSorted(kk);
        nk = numel(groupInd);
        layerMedian(k) = median(zLocal);
        layerMin(k) = zLocal(1);
        support(k) = nk/n;
        subX = 1+floor((xyz(groupInd,1)-cellStartX)/c.gridResolution*c.M);
        subY = 1+floor((xyz(groupInd,2)-cellStartY)/c.gridResolution*c.M);
        subX = max(1,min(c.M,subX));
        subY = max(1,min(c.M,subY));
        coverage(k) = numel(unique(subX+(subY-1)*c.M))/(c.M^2);
        sigma(k) = c.madCorrection*median(abs(zLocal-layerMedian(k)));
        q = exp(-sigma(k)/c.sigma0);
        penalty = c.penalty*(nk<c.Nmin || coverage(k)<c.Cmin);
        score(k) = c.layerWeights*[support(k);coverage(k);q]-penalty;
    end
    [~,chosen] = max(score); % 评分相同时优先选择排序靠前的低高程层
    selectedInitial(cellNumber) = layerMedian(chosen);
    initialLower(cellNumber) = layerMin(chosen);
    layerCache{cellNumber} = struct('median',layerMedian,'minimum',layerMin, ...
'score',score,'originalIndices',ind);
    validLayers = validLayers+K;
end
initialBottom = reshape(selectedInitial,[nRows,nCols]);
madRaster = reshape(allMad,[nRows,nCols]);
referenceHeight = nanNeighbourMedian(initialBottom,c.neighbourhoodSize);
missingRef = isnan(referenceHeight);
referenceHeight(missingRef) = initialBottom(missingRef);
correctedLower = nan(numCells,1);
correctedRepresentative = nan(numCells,1);
for g=1:numGroups
    cellNumber = sortedCells(blockStarts(g));
    info = layerCache{cellNumber};
    hRef = referenceHeight(cellNumber);
    if ~isfinite(hRef)
        hRef = selectedInitial(cellNumber);
    end
    d = abs(info.median-hRef);
    correctedScore = c.continuityWeights(1)*info.score + ...
        c.continuityWeights(2)*exp(-d/c.Tc);
    [~,best] = max(correctedScore);
    correctedLower(cellNumber) = info.minimum(best);
    correctedRepresentative(cellNumber) = info.median(best);
end
finalBottom = reshape(correctedRepresentative,[nRows,nCols]);
noise = xyz(:,3)<correctedLower(ids)-c.Tb;
noise = logical(noise);
diag = struct('nonemptyCells',numGroups,'totalCandidateLayers',validLayers, ...
'meanCandidateLayersPerOccupiedCell',validLayers/max(numGroups,1), ...
'initialLayerMinGrid',reshape(initialLower,[nRows,nCols]), ...
'finalLayerMinGrid',reshape(correctedLower,[nRows,nCols]));
end
function med = nanNeighbourMedian(A,window)
r = floor(window/2);
[nR,nC] = size(A);
stack = nan(nR,nC,window*window);
j=0;
for dr=-r:r
    for dc=-r:r
        j=j+1;
        rTarget = max(1,1-dr):min(nR,nR-dr);
        cTarget = max(1,1-dc):min(nC,nC-dc);
        rSource = rTarget+dr;
        cSource = cTarget+dc;
        sl = nan(nR,nC);
        sl(rTarget,cTarget)=A(rSource,cSource);
        stack(:,:,j)=sl;
    end
end
med = median(stack,3,'omitnan');
end
function [Z,mask] = pointMinimumGrid(z,ids,allowed,nR,nC)
minVals = accumarray(ids(allowed),z(allowed),[nR*nC,1],@min,NaN);
Z = reshape(minVals,nR,nC);
mask = isfinite(Z);
end

function filled = interpolateGridSurface(Z,mask,xg,yg,method)
if ~any(mask(:)), error('Cannot interpolate: no valid grid cells.'); end
if all(mask(:)), filled=Z; return; end
[XX,YY] = meshgrid(xg,yg);
if nnz(mask)<3 || nnz(mask & isfinite(Z))<3
    filled=Z;
    filled(~mask)=Z(find(mask,1,'first'));
    return;
end
try
    F = scatteredInterpolant(XX(mask),YY(mask),Z(mask),method,'nearest');
    filled = F(XX,YY);
catch
    F = scatteredInterpolant(XX(mask),YY(mask),Z(mask),'nearest','nearest');
    filled = F(XX,YY);
end
filled(mask)=Z(mask); % 实测栅格保留原始高程值
if any(~isfinite(filled(:)))
    filled(~isfinite(filled))=median(Z(mask));
end
end
function [mask,Dfinal,S,R,Tlast,history] = stageTerrainAwarePMF(Z,validGrid,c)
% 对四个方向分别进行灰度形态学开运算，取中位数融合。
% 每一级候选栅格与上一级结果取交集，实现渐进筛选。
[dZdX,dZdY] = gradient(Z,c.gridResolution,c.gridResolution);
S = hypot(dZdX,dZdY);         % 局部坡度梯度模长（高程变化量与水平距离之比）
seR = strel('square',c.reliefWindow);
R = imdilate(Z,seR)-imerode(Z,seR);
R=max(R,0);
mask = logical(validGrid);
K = numel(c.radii);
history = table('Size',[K 4], ...
'VariableTypes',{'double','double','double','double'}, ...
'VariableNames',{'ScaleLevel','RadiusCells','CandidateCells','ElapsedSeconds'});
Dfinal=zeros(size(Z));
Tlast = zeros(size(Z));
for k=1:K
    tick=tic;
    radius=c.radii(k);
    opening3d=zeros([size(Z),numel(c.orientations)]);
    for a=1:numel(c.orientations)
        % 将尺度值r解释为半径（栅格数），线形结构元素长度为2r+1。
        se = strel('line',2*radius+1,c.orientations(a));
        opening3d(:,:,a)=imopen(Z,se);
    end
    Zopen=median(opening3d,3);
    D=max(0,Z-Zopen);
    T=c.T0+c.alpha*k+c.beta*S+c.gamma*R;
    mask=mask & (D<=T);
    Dfinal=D;
    Tlast=T;
    history.ScaleLevel(k)=k;
    history.RadiusCells(k)=radius;
    history.CandidateCells(k)=nnz(mask);
    history.ElapsedSeconds(k)=toc(tick);
    fprintf('  PMF radius %3d cells (%5.1f m): %6d candidate cells, %7.2f s\n', ...
        radius,radius*c.gridResolution,nnz(mask),history.ElapsedSeconds(k));
end
if ~any(mask(:))
    warning(['No candidate cells survived all scales. Reverting to cells ' ...
'passing the final scale only. Verify PMF thresholds.']);
    mask=validGrid & (Dfinal<=Tlast);
end
if nnz(mask)<3
    warning('Fewer than 3 candidate cells: using occupied valid cells to construct a prior.');
    mask=validGrid;
end
end
function [candidatePoints,ground,prior,conf,delta,rhoGrid] = ...
    stageTerrainPriorAndRefinement(xyz,ids,validPoints,candidateCells,Z,xg,yg,c)
nR=numel(yg);nC=numel(xg);
% 候选栅格采用该格内有效点的最低高程作为地形参考。
% 栅格高程代表值还可以采用其他统计方式；这里保持既有实现不变。
% 修改此规则可能改变地面点分类和DEM构建结果。
candidatePoints=validPoints & candidateCells(ids);
priorSeed=nan(nR,nC);
mask=candidateCells & isfinite(Z);
priorSeed(mask)=Z(mask);
prior=interpolateGridSurface(priorSeed,mask,xg,yg,'natural');
% 以候选地面栅格的最近距离近似候选地面点的最近距离。
% 距离按栅格分辨率计算，不是精确的点级欧氏距离。
distance=bwdist(mask)*c.gridResolution;
conf=exp(-distance/c.D0);
xq=min(max(xyz(:,1),xg(1)),xg(end));
yq=min(max(xyz(:,2),yg(1)),yg(end));
zPrior=interp2(xg,yg,prior,xq,yq,'linear');
delta=abs(xyz(:,3)-zPrior);
% 按地形先验高程偏差识别高置信候选地面点。
baseGround=candidatePoints & isfinite(delta) & (delta<c.Tf);
ground=baseGround;
radiusCells=max(1,round(c.groundNeighbourhoodRadius/c.gridResolution));
window=ones(2*radiusCells+1);
validRasterCount=reshape(accumarray(ids(validPoints),1,[nR*nC,1],@sum,0),nR,nC);
validNeighbourCount=imfilter(double(validRasterCount),window,0,'same');
for k=1:c.refineIterations
    grdCount=reshape(accumarray(ids(ground),1,[nR*nC,1],@sum,0),nR,nC);
    nearbyGround=imfilter(double(grdCount),window,0,'same');
    rhoGrid=nearbyGround./max(validNeighbourCount,1);
    rhoAtPoints=rhoGrid(ids);
    % 对未通过初筛但满足高程差和邻域地面点比例要求的候选点进行恢复。
    recovered=candidatePoints & ~baseGround & (delta>=c.Tf) & ...
        (delta<c.Tfc) & (rhoAtPoints>=c.rho0);
    nextGround=baseGround | recovered;
    changed=nnz(xor(nextGround,ground));
    ground=nextGround;
    fprintf('  Ground refinement iteration %d: %d ground points (%d changed)\n', ...
        k,nnz(ground),changed);
    if changed==0,break;end
end
% 迭代结束后重新计算邻域地面点比例，以保持导出栅格与最终分类一致。
grdCount=reshape(accumarray(ids(ground),1,[nR*nC,1],@sum,0),nR,nC);
rhoGrid=imfilter(double(grdCount),window,0,'same')./max(validNeighbourCount,1);
end
function [DEM,observed,preSmoothed] = stageDEM(z,ids,ground,prior,xg,yg,c)
nR=numel(yg);nC=numel(xg);
medVals = nan(nR*nC,1);
if any(ground)
    used=ids(ground);
    medVals=accumarray(used,z(ground),[nR*nC,1],@median,NaN);
else
    warning('No final ground points. Final DEM will equal interpolated terrain prior.');
end
Zsample=reshape(medVals,nR,nC);
observed=isfinite(Zsample);
if nnz(observed)>=3
    Zdem=interpolateGridSurface(Zsample,observed,xg,yg,'natural');
else
    Zdem=prior;
end
% 对缺乏实测地面点的栅格使用地形先验高程补充，
% 不将补充的高程视为实际观测点。
preSmoothed=Zdem;
preSmoothed(~observed)=prior(~observed);
% 使用给定权重融合当前DEM和局部均值，实现平滑处理。
window=ones(c.smoothingWindow);
N=imfilter(double(isfinite(preSmoothed)),window,0,'same');
V=preSmoothed;
V(~isfinite(V))=0;
Zmean=imfilter(V,window,0,'same')./max(N,1);
DEM=c.smoothingWeight*preSmoothed+(1-c.smoothingWeight)*Zmean;
end
function stat=summarizeHoles(observed)
% 空洞分类仅用于附加统计，不影响地面点分类或DEM插值。
% 按面积和长宽比将内部空洞划分为三类。
% 仅识别被有效地面栅格包围的内部空洞。
% 与测区边界相连的缺失区域不计入内部空洞。
internal=imfill(observed,'holes') & ~observed;
conn=bwconncomp(internal,8);
stat=struct('count',conn.NumObjects,'small',0,'elongated',0,'large',0);
for i=1:conn.NumObjects
    [r,c]=ind2sub(size(internal),conn.PixelIdxList{i});
    h=max(r)-min(r)+1;w=max(c)-min(c)+1;
    area=numel(r);
    aspect=max(h,w)/max(min(h,w),1);
    if area<=12
        stat.small=stat.small+1;
    elseif aspect>=5
        stat.elongated=stat.elongated+1;
    else
        stat.large=stat.large+1;
    end
end
end
function ref=readReferenceLabels(path,N)
if exist(path,'file')~=2,error('Reference label file is missing: %s',path);end
v=readmatrix(path);
if isempty(v),error('Reference label file is empty.');end
v=v(:,end);
if numel(v)~=N
    error('Reference must have %d labels, corresponding to valid input rows; got %d.',N,numel(v));
end
if ~all(isfinite(v) & (v==0|v==1))
    error('Reference labels must be 0 (non-ground) or 1 (ground).');
end
ref=logical(v);
end
function a=evaluateClassification(pred,ref)
% 根据参考标签计算Ⅰ类误差、Ⅱ类误差、总误差和Kappa系数。
TP=nnz(pred & ref); TN=nnz(~pred & ~ref);
FP=nnz(pred & ~ref); FN=nnz(~pred & ref);
NG=TP+FN;NN=TN+FP; N=numel(pred);
assert(NG>0 && NN>0,'Both ground and non-ground reference labels required.');
a=struct('available',true,'TP',TP,'TN',TN,'FP',FP,'FN',FN);
a.typeI=100*FN/NG;
a.typeII=100*FP/NN;
a.total=100*(FN+FP)/N;
Po=(TP+TN)/N;
Pe=((TP+FP)*(TP+FN)+(FN+TN)*(FP+TN))/(N*N);
a.kappa=(Po-Pe)/(1-Pe);
end
function saveResultFigures(xyz,ground,D,cells,DEM,xg,yg,c,outDir)
% 图像用于过程检查与结果展示，不替代独立参考数据的定量评价。
N=size(xyz,1);
step=max(1,ceil(N/c.maxScatterPoints));
showIdx=1:step:N;
fig=figure('Visible','off','Color','w','Position',[100 100 1100 850]);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
nexttile;
scatter3(xyz(showIdx,1),xyz(showIdx,2),xyz(showIdx,3),1.5,xyz(showIdx,3),'filled');
axis equal;axis tight;view(35,32);grid on;
title('(a) Original point cloud');xlabel('X (m)');ylabel('Y (m)');zlabel('Z (m)');
nexttile;
imagesc(xg,yg,D);axis image;set(gca,'YDir','normal');
colorbar;title('(b) PMF residual (m)');xlabel('X (m)');ylabel('Y (m)');
nexttile;
imagesc(xg,yg,double(cells));axis image;set(gca,'YDir','normal');
colormap(gca,[0.23 0.20 0.58;0.97 0.84 0.27]);
cb=colorbar;cb.Ticks=[0,1];
title('(c) Candidate ground cells');xlabel('X (m)');ylabel('Y (m)');
nexttile;
hold on;
non=showIdx(~ground(showIdx));
gnd=showIdx(ground(showIdx));
scatter3(xyz(non,1),xyz(non,2),xyz(non,3),1.5,[.72 .72 .72],'filled');
scatter3(xyz(gnd,1),xyz(gnd,2),xyz(gnd,3),2.0,[.05 .55 .25],'filled');
hold off;axis equal;axis tight;view(35,32);grid on;
title('(d) Ground classification');xlabel('X (m)');ylabel('Y (m)');zlabel('Z (m)');
exportgraphics(fig,fullfile(outDir,'Figure14_Linfen_four_panels.png'),'Resolution',200);
close(fig);
fig=figure('Visible','off','Color','w','Position',[100 100 1000 700]);
imagesc(xg,yg,DEM);axis image;set(gca,'YDir','normal');colorbar;
title('Reconstructed urban DEM (m)');xlabel('X (m)');ylabel('Y (m)');
exportgraphics(fig,fullfile(outDir,'DEM_final.png'),'Resolution',200);
close(fig);
end
