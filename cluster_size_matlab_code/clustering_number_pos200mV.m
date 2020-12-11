clc;clear all;close all;

load('E:\LuoLAB\cx26\cAMP\27cAMP_strip\clustering_number\pos200mv_cAMP_COM.mat')

thr = 20 % dist= 10.55; 2nd well of cAMP atom N3, rdf distance

Tall = zeros(9919, 1);
csize = zeros(9919, 1);
count = zeros(9919, 1);
for i = 1 : 9919
    i
    space = zeros(27, 3);
    for j = 0 : 26
        if j < 10
            eval(['space(j+1, :) = pos200mv_cAMP_COM{1}.Vec_0000' num2str(j) '(i, :);']);
        else
            eval(['space(j+1, :) = pos200mv_cAMP_COM{1}.Vec_000' num2str(j) '(i, :);']);
        end
    end
    space(:, 3) = space(:, 3) - 74; % offset by protein COM z =74A
    
    index = find(space(:, 3) >= -30 & space(:, 3) < -20);
    %  for each z-axis 10A-width window, 
    %  change here [30,40] to the target z-axis 10A-width window
    if length(index) < 2
        continue,
    end
    subspace = space(index, :);
    
    
    Y = pdist(subspace);
    Ys = squareform(Y);
    Z = linkage(Y);
    T = cluster(Z,'Cutoff',thr,'Criterion','distance');
    uT = unique(T);
    thisd = 0;
    cnt = 0;
    for j = 1 : length(uT)
        thisc = find(T == uT(j));            
        if length(thisc) >= 2
            thisd = thisd + mean(mean(Ys(thisc, thisc)));
            cnt = cnt + 1;
        end
    end
    if cnt ~= 0
        csize(i) = thisd / cnt;
    end
            
    count(i)=cnt;        
    Tall(i) = length(unique(T));
end
cavg=mean(csize);
avg=mean(count);
['average cluster size in Angstrom =:' num2str(cavg) ]
['average number of clusters =:' num2str(avg) ]
