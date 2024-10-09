% find and remove voxels with significant activation differences in HSC and LSC (p < .05, uncorrected) 
% then compare rd
% this method is based on Xue 2010 science paper

%%
clear;clc;
subjects = {'sub-02';'sub-03';'sub-04'};
masks = {'Hippocampus_L.nii'};
study_path='I:\FLXX1\Dimension\pattern'; 
roi_path='I:\FLXX1\mask';
n_subjects=numel(subjects);
msk = masks{1};

%% how many voxels?
for s = 1
    sub = subjects{s};
    sub_path=fullfile(study_path,sub);
    mask_fn=fullfile(roi_path,msk);
    %%
    data_HSC=fullfile(sub_path,'glm_T_stats_HSC.nii');
    ds_HSC=cosmo_fmri_dataset(data_HSC,'mask',mask_fn);
    d_HSC=ds_HSC.samples;

    nvoxels = size(d_HSC,2);
end

%%
mymask =[];
counter = 1;

%%
for n = 1: nvoxels
    allh=[];
    alll=[];
    count =1;
    for s = 1:n_subjects

        sub = subjects{s};
        sub_path=fullfile(study_path,sub);
        mask_fn=fullfile(roi_path,msk);
        %%
        data_HSC=fullfile(sub_path,'glm_T_stats_HSC.nii');
        ds_HSC=cosmo_fmri_dataset(data_HSC,'mask',mask_fn);
        d_HSC=ds_HSC.samples;
        hscdat = mean(d_HSC);
        h = hscdat(1,n);
        data_LSC=fullfile(sub_path,'glm_T_stats_LSC.nii');
        ds_LSC=cosmo_fmri_dataset(data_LSC,'mask',mask_fn);
        d_LSC=ds_LSC.samples;
        lscdat = mean(d_LSC);
        l = lscdat(1,n);

        allh(count,1)=h;
        alll(count,1)=l;
        count = count + 1;
    end
    [h,p,ci,stats]=ttest(allh,alll)

    if h==1
        mymask(1,counter)=0;
    else
        mymask(1,counter)=1;
    end
    counter =counter +1;

end

save voxelmask.mat mymask


%%



clear;clc;
subjects = {'sub-02';'sub-03';'sub-04';'sub-05';'sub-06';'sub-07';'sub-08';'sub-09';'sub-10';'sub-11';'sub-12';'sub-13';'sub-14';'sub-15';'sub-17';'sub-18';'sub-19';'sub-20';'sub-21';'sub-25';'sub-28';'sub-29';'sub-30'};
masks = {'Hippocampus_L.nii'};
study_path='I:\FLXX1\Dimension\pattern'; 
roi_path='I:\FLXX1\mask';
msk = masks{1};

alldata = [];
counterp = 0;
%%
load voxelmask.mat

%%
for ev = 70:90

counter=0;
counter1=0;

for s = 1:length(subjects)
    %%
    sub = subjects{s};
    sub_path=fullfile(study_path,sub);
    mask_fn=fullfile(roi_path,msk);
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data_HSC=fullfile(sub_path,'glm_T_stats_HSC.nii');
    ds_HSC=cosmo_fmri_dataset(data_HSC,'mask',mask_fn);
    
    d_HSC=ds_HSC.samples(:,mymask==1)';
    
    data_LSC=fullfile(sub_path,'glm_T_stats_LSC.nii');
    ds_LSC=cosmo_fmri_dataset(data_LSC,'mask',mask_fn);
    
    d_LSC=ds_LSC.samples(:,mymask==1)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    true_length = min(size(d_HSC,2),size(d_LSC,2));
    %% Or we can manually set a standard, which can sometimes increase the robustness of the results
    % true_length=25;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  true_length ~= size(d_HSC,2)
        
        allrd = 0;
        for x = 1:1000
            da_HSC = d_HSC(:,randperm(size(d_HSC,2),true_length)); %randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n.
            [coeff,score,latent,tsquared,explained,mu]=pca(da_HSC);
                total_explained_hsc=0;
                 RD_var_hsc=0;
            for n = 1:length(explained)
                total_explained_hsc=total_explained_hsc+explained(n);
                if total_explained_hsc > ev
                    break
                end
            end
            
            RD_var_hsc = n;
            allrd = allrd + RD_var_hsc;
            
        end
        
        counter=counter+1;
        RDhsc(counter,:)=allrd/1000;
        
    else
        %%%%%%
        %
        [coeff,score,latent,tsquared,explained,mu]=pca(d_HSC);
        total_explained_hsc=0;
         RD_var_hsc=0;
        for n = 1:length(explained)
            total_explained_hsc=total_explained_hsc+explained(n);
            if total_explained_hsc > ev
                break
            end
        end
        
        RD_var_hsc = n;
        counter=counter+1;
        RDhsc(counter,:)=RD_var_hsc;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% 
    if  true_length ~= size(d_LSC,2)
        
        allrd = 0;
        for x = 1:1000
            da_LSC = d_LSC(:,randperm(size(d_LSC,2),true_length)); %randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n.
            [coeff,score,latent,tsquared,explained,mu]=pca(da_LSC);
            total_explained_lsc=0;
             RD_var_lsc=0;
            for n = 1:length(explained)
                total_explained_lsc=total_explained_lsc+explained(n);
                if total_explained_lsc > ev
                    break
                end
            end
            
            RD_var_lsc = n;
            allrd = allrd + RD_var_lsc;
            
        end
        
        counter1=counter1+1;
        RDlsc(counter1,:)=allrd/1000;
        
    else
        %%%%%%
        [coeff,score,latent,tsquared,explained,mu]=pca(d_LSC);
         total_explained_lsc=0;
          RD_var_lsc=0;
        for n = 1:length(explained)
            total_explained_lsc=total_explained_lsc+explained(n);
            if total_explained_lsc > ev
                break
            end
        end
        
        RD_var_lsc = n;
        counter1=counter1+1;
        RDlsc(counter1,:)=RD_var_lsc;
    end
    
end

[h,p,ci,stats]=ttest(RDhsc,RDlsc);
counterp = counterp + 1;
allp(counterp,:) = p;
alldata = [alldata,RDhsc,RDlsc];

end
