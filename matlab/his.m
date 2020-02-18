function [matHIS,matHISPerAssay] = his(matData,intMin,intMax,intSteps,strTail,boolDraw)
% Help for HIS()
%
% calculates the Hierarchical Interaction Score between each row of
% matData.
%
% Usage:
%    
%    matHIS = his(matData,intMin,intMax,intSteps)
%
% Only matData is required input. The rest defaults to the following
% values:
%
% intMin   (default = 1.5) is the lower (inner) threshold. Note that
% thresholds are mirrored and that non-phenotypes are supposed to be
% between -intMin and intMin.
%
% intMax   (default = 5) is the upper (outer) threshold.
%
% intSteps (default = 200) is the total number of steps over which the
% calculation runs. As the thresholds are mirrored, 200 steps means 100
% steps from -intMax to -intMin, and 100 steps from intMin to intMax.
%
% Note that the his assumes your data to be centered around 0. This is
% important to get the direction of the thresholding right. The direction
% of the thresholds are defined by: 
% 
%           (x >= t, if x > 0) and (x <= t, if x <= 0).
%
% Optionall, strTail can be 'both' (default) or 'single', to make the steps
% strictly go from intMin to intMax in intSteps steps.
% 
% matHIS = his(matData,intMin,intMax,intSteps,strTail)
% 
%
% (c) 2012, Berend Snijder.
%           bsnijder@gmail.com

if nargin<1
    matData = randn(2000,20);
end

if nargin<2
    intMin=1.5;
end
if nargin<3
    intMax=5;
end
if nargin<4
    intSteps=200;
end
if nargin<5
    strTail = 'both';
end
if nargin==0
    boolDraw = false;
elseif nargin<6
    boolDraw = false;
end

% Limit the data to those rows and columns with at least one phenotype
% (abs(value)) above lower threshold.
matOkRowIX = any(abs(matData)>=intMin,2);
matOkColIX = any(abs(matData)>=intMin,1);
matData = matData(matOkRowIX,matOkColIX);
 
% If there is no data left, we are done here. exit
if isempty(matData)
    % return zeros
    matHIS = sparse(numel(matOkRowIX),numel(matOkRowIX));
    return
end
 
% Calculate thresholds, going from -intMax to -intMin in round(intSteps/2))
% steps, and from intMin to intMax in another round(intSteps/2)) steps.
if strcmpi(strTail,'single')
    matTs = linspace(intMin,intMax,round(intSteps));
else
    matTs = [linspace(-intMax,-intMin,round(intSteps/2)),linspace(intMin,intMax,round(intSteps/2))];
end
    

% This map contains the intermediate statistic for pairwise interactions for
% rows (after data discarding), of the lower amount of hits present in each
% nested hit relation between variables, summed over all thresholds. The
% final statistic is this value divided by (numel(matTs)/2).
boolDoPerEdge = nargout==2;
iMaxVal = numel(matTs);
[iNodes,iAssays] = size(matData);
if boolDoPerEdge
    if iMaxVal < 255%intmax('uint8')
        matSumMinHitCounts = zeros(iNodes,iNodes,iAssays,'uint8');
    elseif iMaxVal < 65535%intmax('uint16')
        matSumMinHitCounts = zeros(iNodes,iNodes,iAssays,'uint16');
    else
        matSumMinHitCounts = zeros(iNodes,iNodes,iAssays,'single');
    end    
else
    matSumMinHitCounts = zeros(iNodes,iNodes,'uint16');
end
% matSumMinHitCounts = sparse(sum(matOkRowIX),sum(matOkRowIX));
 
if boolDraw
    hF = figure();
    matTStats = NaN(2,numel(matTs));
    [matDataHist,matDataHistEdges] = hist(matData(:),30);
end

% loop over all (negative and positive) thresholds
for t = 1:numel(matTs);
 
    if boolDraw
        matPreviousResult = matSumMinHitCounts;
    end
    
    % find hits at current threshold
    if matTs(t)<=0
        matHits = matData<=matTs(t);
    else
        matHits = matData>=matTs(t);
    end
 
    % look for all present hit-patterns (note that this is dependent on the
    % legacy behaviour of unique).
    if ~verLessThan('matlab', '7.14')
        [matPatterns,~,matGenePatternMapping]=unique(matHits,'rows','legacy');
    else
        [matPatterns,~,matGenePatternMapping]=unique(matHits,'rows');
    end
 
    % keep track of number of hits in each pattern
    matPatternHitCount = sum(matPatterns,2);
 
    % find children & parent relationships between each unique hit-pattern
    boolChildren = false(size(matPatterns,1));
    boolParents = false(size(matPatterns,1));
    for i = 1:size(matPatterns,1)
        x2 = matPatterns(i,:);
        matChildCount = sum(matPatterns(:,x2),2);
        matNotChildCount = matPatternHitCount - matChildCount;
        
        % find children (i.e. subsets of current hit-pattern) (note that
        % this includes the pattern itself)
        matCandidateChildren = matChildCount>0 & matNotChildCount==0;
        % find parents
        matCandidateParents = matChildCount==sum(x2) & matNotChildCount>0;
 
        % store parent & child relations
        boolChildren(i,matCandidateChildren) = true;
        boolParents(i,matCandidateParents) = true;
    end
    
    % Store the minimum-hit-count per hit-pattern relation (after
    % transitive reduction) for all genes present with those hit-patterns.
    
    % (Is there a matrix operation that achieves the same?)
    for i = 1:size(matPatterns,1)
        % a link is drawn from current set to candidate child sets that do not
        % have parents that are also children of the current set
        c = find(boolChildren(i,:));
 
        % for each child-set, see if it has parents that are children of
        % current set
        for c2 = c
            % genes with the same hit pattern are both child and parent of
            % eachother, so make sure the same pattern class does not
            % exclude a linkage from being drawn.
            
            % Say "hello" to the slowest line in current code:
            p = find(boolParents(c2,:));
            p(p==i) = [];
            
            if ~any(ismembc(p,c))
                % intMinPatternHitCount = min(matPatternHitCount([i,c2]));
                xi = matGenePatternMapping==i;
                xj = matGenePatternMapping==c2;
                if boolDoPerEdge
                    matSumMinHitCounts(xi,xj,matPatterns(c2,:)) = matSumMinHitCounts(xi,xj,matPatterns(c2,:)) + 1;
                else
                    matSumMinHitCounts(xi,xj) = matSumMinHitCounts(xi,xj) + matPatternHitCount(c2);
                end
            end
        end
    end
    
    if boolDraw
        matSumMinHitCounts(logical(eye(size(matSumMinHitCounts,1)))) = 0;
        matTStats(1,t) = sum(matSumMinHitCounts(:)) - sum(matPreviousResult(:));
        matTStats(2,t) = size(matPatterns,1);
        
        figure(hF);
        clf
        subplot(2,3,1)
        imagesc(matPatterns)
        title('matPatterns')
        subplot(2,3,2)
        imagesc(matSumMinHitCounts)
        title('intermediate HIS')
        subplot(2,3,3)
        bar(matDataHistEdges,matDataHist,'hist')
        matYLim = get(gca,'YLim');
        line([matTs(t),matTs(t)],matYLim)
        title('matData')
        subplot(2,3,4)
        imagesc(boolParents)
        title('boolParents')
        subplot(2,3,5)
        imagesc(boolChildren)
        title('boolChildren') 
        subplot(2,3,6)
        plotyy(matTs,matTStats(1,:),matTs,matTStats(2,:));
        xlabel('threshold')
        drawnow
        
    end
end

if boolDoPerEdge
    matHISPerAssay = zeros(iNodes,iNodes);
    matHISPerAssay(matOkRowIX,matOkRowIX,matOkColIX) = matSumMinHitCounts;
    
    % and sum over all assays
    matSumMinHitCounts = sum(matSumMinHitCounts,3);
end


% remove HIS self links
matSumMinHitCounts(logical(eye(size(matSumMinHitCounts,1)))) = 0;

if ~boolDoPerEdge
    % check if everything is ok
    if any(matSumMinHitCounts(:)>=intmax(class(matSumMinHitCounts)))
        warning('bs:Bla','We should increase data-type of matSumMinHitCounts. Max value reached...')
    end
end 

% clear up some memory...
clear matData matPatterns matGenePatternMapping
 
% Get indices (a & b) and corresponding values (c) of summed min hit counts
[a,b,c] = find(matSumMinHitCounts);
clear matSumMinHitCounts
 
% Remap back to the full dataset, before discarding of rows.
d = find(matOkRowIX);
 
% We can return a sparse interaction matrix.
matHIS = sparse(d(a),d(b),double(c)/numel(matTs),numel(matOkRowIX),numel(matOkRowIX),numel(c));
 
% Or we return a list of interactions.
% matHIS = cat(2,d(a),d(b),double(c)/numel(matTs));
end
