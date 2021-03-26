clc;clear all;close all;

datapath = 'E:\LuoLAB\cx26\milestoning_cAMP';
list = dir(fullfile(datapath, 'cav_wid*.csv'));

ref = readtable(fullfile(datapath, 'cav_nonV.csv'));

zmin = min(ref.Var1);
zmax = max(ref.Var1);

zmin = -40;
zmax = 40;

cohend = zeros(length(list), zmax-zmin+1);
for loop = 1 : length(list)
    loop
    dnow = readtable(fullfile(datapath, list(loop).name));
    adj = mean(dnow.Var1) - mean(ref.Var1);
    dnow.Var1 = dnow.Var1 - adj;
    a=size(dnow,1);
    writetable(dnow,[num2str(loop) 'window.csv']);
    flag = 0;
    for z = zmin : zmax-1
        flag = flag + 1;
        refz = find(ref.Var1>=z & ref.Var1<z+1);
        dz = find(dnow.Var1>=z & dnow.Var1<z+1);
        
        nref = length(refz);
        nd = length(dz);
        refv = ref.Var2(refz);
        dv = dnow.Var2(dz);
        
        psd = sqrt(((nref-1)*std(refv)^2 + (nd-1)*std(dv)^2) / (nref+nd-2));
        if isnan(psd)
            cohend(loop, flag) = 0;
        else
            cohend(loop, flag) = (mean(dv)-mean(refv)) / psd;
        end
    end
       
end

save('cohenD_mat_wid03_42.txt','cohend', '-ascii', '-double', '-tabs')
type 'cohenD_mat_wid03_42.txt'
