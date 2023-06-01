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