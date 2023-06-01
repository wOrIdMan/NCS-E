function [recordedBestX,FeasibleObj1st,FeasibleFes1st] = NCS_E(MaxFes,evaluateFunc,funcNum,lb,ub,m,dim,input)
%% 
epsim = 1e-4;
FeasibleX1st = nan(1,dim);
FeasibleFes1st = nan;
FeasibleObj1st = nan;
populationX = rand(m,dim).*(ub-lb)+lb;
[populationObj, populationG, populationH] = evaluateFunc(populationX,funcNum);
populationConV = calVio(m,populationG,populationH,epsim);
Fes = 1;
sigmaLb = ones(1,dim ).*1e-10;
sigmaUb = (ub-lb)./5;
sigma = sigmaUb.*ones(m,dim);
% sigma = ones(m,dim);
r = input.r;
epoch = input.epoch;
generation = 1;
recordedBesty = inf(1,MaxFes);
recordedBestConV = inf(1,MaxFes);
bestX=nan;
recordedBestX = nan(10,dim);
[recordedBesty,recordedBestConV,bestX,Fes,FeasibleX1st,FeasibleObj1st,FeasibleFes1st] = record(populationX,populationObj,populationConV ,recordedBesty,recordedBestConV,bestX,Fes,FeasibleX1st,FeasibleObj1st,FeasibleFes1st);
successCount = zeros(m,1);
maxGen = floor(MaxFes/m);
recordedSigma = zeros(m,maxGen);
recordedPenalty = zeros(1,maxGen);
theta=0.9;
[~,sortedVioIndex] = sort(populationConV);
v0 = mean(populationConV(sortedVioIndex(1:floor(theta*m))));
T = maxGen*0.35;
% T = 1000;
gama = max(3, (-5-log(v0))/log(0.05));
archiveMaxLen = 10*m;
archiveLen = m;
archiveX(1:archiveLen,:) = populationX;
archiveObj(1:archiveLen,1) = populationObj;
archiveG(1:archiveLen,:) = populationG;
archiveH(1:archiveLen,:) = populationH;
vBestX = nan(1,dim);
vBestObj = inf;
vBestG = inf;
vBestH = inf;

TC=maxGen*0.5;
% TC = 1000;
n=ceil(0.9*m);   
sortedVio = sort(populationConV);
EPSILON = mean(sortedVio(1:n));
CP=max(3,(-5-log(EPSILON))/log(0.05));
voterNum = 3;
voteW = ones(m ,voterNum)/voterNum;
recordVoteW = nan(maxGen,voterNum,m );
recordVoteW(generation,:,:) = voteW';
recordSuccVote = zeros(m ,voterNum);
updateDelta = 0.01;
newVoteW = zeros(m ,voterNum);
updateT = ones(1,m);
lr = 0.01;
while(Fes<MaxFes)
    recordedBestX(ceil(Fes/MaxFes*10),:) = bestX;
    lambda = 1 + randn*max(0.1*(1-Fes/MaxFes),0);
    sigmaLb = ones(1,dim).*10^(-8*sqrt(Fes/MaxFes));
    temper = 10^(-2*(Fes/MaxFes));
    if generation<T
        v = max(0,v0*(1-generation/T).^gama);
    else
        v=0;
    end
    if generation>TC
        epsilon = 0;
    else
        epsilon = EPSILON*((1-generation/TC)^CP);
    end
    
    offspringX = NCSoffs(populationX,sigma,lb,ub,Fes,MaxFes);
    [offspringObj, offspringG, offspringH] = evaluateFunc(offspringX,funcNum);
    offspringConV = calVio(m, offspringG, offspringH, epsim);
    [recordedBesty,recordedBestConV,bestX,Fes,FeasibleX1st,FeasibleObj1st,FeasibleFes1st] = record(offspringX,offspringObj,offspringConV,recordedBesty,recordedBestConV,bestX,Fes,FeasibleX1st,FeasibleObj1st,FeasibleFes1st);
    
    
    
    
    [archiveX, archiveObj, archiveG, archiveH, archiveLen] = updateArchive(archiveX, archiveObj, archiveG, archiveH, offspringX, offspringObj, offspringG, offspringH, archiveLen, archiveMaxLen);
    archiveRelaxedConV = RelaxConV(archiveLen, archiveG, archiveH, epsim, v);
    vBestRelaxedConV = RelaxConV(1, vBestG, vBestH, epsim, v);
    for i = 1:archiveLen
        if vBestRelaxedConV > archiveRelaxedConV(i) || (vBestRelaxedConV == archiveRelaxedConV(i) && vBestObj>archiveObj(i))
            vBestX = archiveX(i,:);
            vBestObj = archiveObj(i);
            vBestG = archiveG(i,:);
            vBestH = archiveH(i,:);
            vBestRelaxedConV = archiveRelaxedConV(i);
        end
    end
    
    
    PenaltyF = max(0, max((vBestObj-archiveObj)./(archiveRelaxedConV-vBestRelaxedConV)));
    recordedPenalty(generation) = PenaltyF;
    populationRelaxedConV = RelaxConV(m, populationG, populationH, epsim, v);
    offspringRelaxedConV = RelaxConV(m, offspringG, offspringH, epsim, v);
    populationFit = populationObj + PenaltyF.*populationRelaxedConV; 
    offspringFit = offspringObj + PenaltyF.*offspringRelaxedConV;
    populationFitNorm = populationFit-min([populationFit;offspringFit]);
    offspringFitNorm = offspringFit-min([populationFit;offspringFit]);
    offspringFitNorm = (offspringFitNorm+1e-9) ./ (populationFitNorm+offspringFitNorm+2*1e-9);

    
    populationObjNorm = (populationObj-min([populationObj;offspringObj]));
    offspringObjNorm = (offspringObj-min([populationObj;offspringObj]));
    [populationCorr,offspringCorr] = Corr(populationX,offspringX,sigma);

    offspringObjNorm = (offspringObjNorm+1e-9) ./ (populationObjNorm+offspringObjNorm+2*1e-9);
    offspringConVNorm = (offspringConV+1e-9) ./ (populationConV+offspringConV+2*1e-9);
    offspringCorrNorm = (offspringCorr+1e-9)./(populationCorr+offspringCorr+2*1e-9);
    
    vote = zeros(m ,voterNum);
    for i = 1:m
        % Feasibility rule
        if populationConV(i) == offspringConV(i)
            vote(i,1) = (offspringObjNorm(i) < lambda * offspringCorrNorm(i))*voteW(i,1);
        else
            vote(i,1) = (offspringConVNorm(i) < lambda * offspringCorrNorm(i))*voteW(i,1);
        end
        % epsilon method
        if (populationConV(i)<=epsilon) && (offspringConV(i)<=epsilon)
            vote(i,2) = (offspringObjNorm(i) < lambda * offspringCorrNorm(i)) * voteW(i,2);
        else
            vote(i,2) = (offspringConVNorm(i) < lambda * offspringCorrNorm(i)) * voteW(i,2);
        end
        % v-level penalty
        vote(i,3) = (offspringFitNorm(i) < lambda * offspringCorrNorm(i)) * voteW(i,3);

        isreplace = sum(vote(i,:))>sum(voteW(i,:)-vote(i,:));

        if isreplace
            successCount(i) = successCount(i)+1;
            populationX(i,:) = offspringX(i,:);
            populationObj(i) = offspringObj(i);
            populationG(i,:) = offspringG(i,:);
            populationH(i,:) = offspringH(i,:);
            populationConV(i) = offspringConV(i);
        end
        voteCorr = zeros(1,voterNum);
        for j = 1:voterNum
            if vote(i,j)>0
                voteCorr(j) = offspringCorr(i);
            else
                voteCorr(j) = populationCorr(i);
            end
        end
        newVoteW(i,:) = newVoteW(i,:) + voteCorr/sum(voteCorr);
    end
    
    %% 1/5 rule
    if mod(generation,epoch)==0
        for i = 1:m
            if successCount(i)/epoch > 0.2
                sigma(i,:) = sigma(i,:)/r;
            elseif successCount(i)/epoch < 0.2
                sigma(i,:) = sigma(i,:)*r;
            end
            sigma(i,sigma(i,:)>sigmaUb) = sigmaUb(sigma(i,:)>sigmaUb);
            sigma(i,sigma(i,:)<sigmaLb) = sigmaLb(sigma(i,:)<sigmaLb);
        end
        successCount = zeros(m,1);
        recordedSigma(:,generation) = sigma(:,1);

        
        for i = 1:m
            voteW(i,:) = lr*newVoteW(i,:)/epoch + (1-lr)*voteW(i,:);
            recordSuccVote(i,:) = zeros(1,3);
            newVoteW(i,:) = zeros(1,3);
            
        end
    end

    recordVoteW(generation,:,:) = voteW';
    

    generation = generation+1;
end
%%
recordedBestX(end,:) = bestX;
recordedBestX = recordedBestX';
end


function [popCorr,offsCorr] = Corr(pop,offs,sigma)

[m,~] = size(pop);
popCorr = zeros(m,1);
offsCorr = zeros(m,1);
for i = 1:m
    popDist = zeros(1,m);
    offsDist = zeros(1,m);
    for j = [1:i-1,i+1:m]
        popDist(j) = BD(pop(i,:),sigma(i,:),pop(j,:),sigma(j,:)); % Bhattacharyya distance
        offsDist(j) = BD(offs(i,:),sigma(i,:),pop(j,:),sigma(j,:));
    end
    popDist(i) = [];
    offsDist(i) = [];
    popCorr(i) = min(popDist);
    offsCorr(i) = min(offsDist);
end
end
function dist = BD(x1,sigma1,x2,sigma2)
% dist = sum((x1-x2).^2./(sigma1.^2+sigma2.^2+1e-9))/4+log((prod((sigma1.^2+sigma2.^2)./2)+1e-100)/(sqrt(prod(sigma1.^2)*prod(sigma2.^2))+1e-100))/2;
dist = sum((x1-x2).^2./(sigma1.^2+sigma2.^2+1e-9))/4 + 0.5*sum(log((sigma1.^2+sigma2.^2)/2) - log(sigma1) - log(sigma2));
end





function [recordedBesty,recordedBestConV,bestx,Fes,FeasibleX1st,FeasibleObj1st,FeasibleFes1st] = record(offs,offsf,offsVio,recordedBesty,recordedBestConV,bestx,Fes,FeasibleX1st,FeasibleObj1st,FeasibleFes1st)
m = length(offsf);
for i = 1:m
    if isnan(FeasibleFes1st)&&offsVio(i)==0
        FeasibleFes1st = Fes;
        FeasibleObj1st = offsf(i);
        FeasibleX1st = offs(i,:);
    end
    if (offsf(i)<recordedBesty(Fes)&&offsVio(i)==recordedBestConV(Fes))||(offsVio(i)<recordedBestConV(Fes))
        bestx = offs(i,:);
        recordedBesty(Fes) = offsf(i);
        recordedBestConV(Fes) = offsVio(i);
    end
    Fes = Fes + 1;
    recordedBesty(Fes) = recordedBesty(Fes-1);
    recordedBestConV(Fes) = recordedBestConV(Fes-1);
end
end

function [archiveX, archiveObj, archiveG, archiveH, archiveLen] = updateArchive(archiveX, archiveObj, archiveG, archiveH, offspringX, offspringObj, offspringG, offspringH, archiveLen, archiveMaxLen)
m = size(offspringX,1);
for i = 1:m
    if archiveLen<archiveMaxLen
        archiveLen = archiveLen+1;
        replaceedI = archiveLen;
    else
        replaceedI = randi(archiveMaxLen);
    end
    archiveX(replaceedI,:) = offspringX(i,:);
    archiveObj(replaceedI) = offspringObj(i);
    archiveG(replaceedI,:) = offspringG(i,:);
    archiveH(replaceedI,:) = offspringH(i,:);
end

end
function RelaxedConV = RelaxConV(m,g,h,epsim,v)
RelaxedConV = zeros(m,1);
g(g<v) = 0;
if ~isempty(g)
    RelaxedConV = RelaxedConV + sum(g,2);
end
h = abs(h);
h(h<=epsim+v) = 0;
if ~isempty(h)
    RelaxedConV = RelaxedConV + sum(h,2);
end
end