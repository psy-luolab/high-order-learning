% loading similarity analysis
% we fit and transform all data into one pca space
% then look into the loadings of different data

clear;clc;
subjects = {'sub-01';'sub-02';'sub-03'};
masks = {'？.nii'}; 
s1_path='H:\GJXX_2_reanalysis\C\First_level\S1\rsa';
s2_path='H:\GJXX_2_reanalysis\C\First_level\S2\rsa';

roi_path='H:\GJXX_2_reanalysis\mask';
n_subjects = numel(subjects);
n_masks = numel(masks);

alldata = [];
counterp = 0;

for EV = 70:90
    counter = 0;
for s = 1:length(subjects)
    
    sub = subjects{s};
    sub_s1_path=fullfile(s1_path,sub);
    sub_s2_path=fullfile(s2_path,sub);
    mask1_fn=fullfile(roi_path,masks{1});
    
    %%  
    S1data_HSC=fullfile(sub_s1_path,'glm_HSC_NA.nii');
    S1ds_HSC=cosmo_fmri_dataset(S1data_HSC,'mask',mask1_fn);
    S1data_LSC=fullfile(sub_s1_path,'glm_LSC_NA.nii');
    S1ds_LSC=cosmo_fmri_dataset(S1data_LSC,'mask',mask1_fn);
    S2data_HSC=fullfile(sub_s2_path,'glm_HSC_NA.nii');
    S2ds_HSC=cosmo_fmri_dataset(S2data_HSC,'mask',mask1_fn);
    S2data_LSC=fullfile(sub_s2_path,'glm_LSC_NA.nii');
    S2ds_LSC=cosmo_fmri_dataset(S2data_LSC,'mask',mask1_fn);
    %%
    H = size(S1ds_HSC.samples,1);
    L = size(S1ds_LSC.samples,1);
    %%
    allds = cosmo_stack({S1ds_HSC,S1ds_LSC,S2ds_HSC,S2ds_LSC});
    
    input = allds.samples';
    
    [coeff,score,~,~,explained,~]=pca(input);
     
    total_explained=0;
    
    for m = 1:length(explained)
        total_explained = total_explained + explained(m);
        if total_explained >= EV
            break
        end
    end
    
    loadings_hsc_s1 = coeff(1:H,:);
    loadings_lsc_s1 = coeff(H+1:H+L,:);
    loadings_hsc_s2 = coeff(H+L+1:H+L+H,:);
    loadings_lsc_s2 = coeff(H+L+H+1:H+L+H+L,:);
    
    %%
    r_hsc = 0;
    for x = 1:m 
        r = corr(loadings_hsc_s1(:,x),loadings_hsc_s2(:,x),'Type','Spearman');
        r_hsc = r_hsc + r;
       
    end
    
    %% 
    r_lsc = 0;
    for x = 1:m 
        r = corr(loadings_lsc_s1(:,x),loadings_lsc_s2(:,x),'Type','Spearman');
        r_lsc = r_lsc + r;
        
    end
    counter = counter + 1;
    HSC_data(counter,:) = r_hsc;
    LSC_data(counter,:) = r_lsc;
    
end

[h1,p1,ci1,stats1] = ttest(HSC_data,LSC_data);
counterp = counterp + 1;
allp(counterp,:) = p1;  
alldata = [alldata,HSC_data,LSC_data];

end

%save data.mat allp alldata
