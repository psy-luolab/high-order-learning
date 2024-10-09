% do pca dimension reduction based on sheng2022 science advances paper
%
%
% prepare your 4d.nii after single pattern estimation
% including activation patterns (trial*voxel) under HSC/LSC condition


% The toolbox should be installed:
%    1. cosmoMvpa;

%%
clear;clc;
subjects = {'sub-02';'sub-03';'sub-04';'sub-05';'sub-06'};
masks = {'?.nii'};
study_path='I:\FLXX1\Dimension\pattern';
roi_path='I:\FLXX1\mask';
msk = masks{1};
%%
alldata = [];
allp = [];
counterp = 0;
%%
for ev = 80

    counter_hsc=0;
    counter_lsc=0;

    for s = 1:length(subjects)
        total_explained_hsc=0;
        total_explained_lsc=0;
        %%
        sub = subjects{s};
        sub_path=fullfile(study_path,sub);
        mask_fn=fullfile(roi_path,msk);
        %%
        data_HSC=fullfile(sub_path,'glm_T_stats_HSC.nii');
        ds_HSC=cosmo_fmri_dataset(data_HSC,'mask',mask_fn);
        dsm_HSC=cosmo_pdist(ds_HSC.samples, 'correlation');
        RDM_HSC=cosmo_squareform(dsm_HSC);
        [~,~,explained]=pcacov(RDM_HSC);


        for n = 1:length(explained)
            total_explained_hsc=total_explained_hsc+explained(n);
            if total_explained_hsc > ev
                break
            end
        end

        RD_var_hsc = n;

        counter_hsc=counter_hsc+1;

        RDhsc(counter_hsc,:)=RD_var_hsc;

        %%
        data_LSC=fullfile(sub_path,'glm_T_stats_LSC.nii');
        ds_LSC=cosmo_fmri_dataset(data_LSC,'mask',mask_fn);
        dsm_LSC=cosmo_pdist(ds_LSC.samples, 'correlation');
        RDM_LSC=cosmo_squareform(dsm_LSC);

        [coeff,latent,explained]=pcacov(RDM_LSC);


        for j = 1:length(explained)
            total_explained_lsc=total_explained_lsc+explained(j);
            if total_explained_lsc > ev
                break
            end
        end

        RD_var_lsc = j;
        counter_lsc=counter_lsc+1;
        RDlsc(counter_lsc,:)=RD_var_lsc;

    end

    [h,p,ci,stats]=ttest(RDhsc,RDlsc);
    counterp = counterp + 1;
    allp(counterp,:) = p;
    alldata = [alldata,RDhsc,RDlsc];

end

%save sheng2022_ev80.mat allp alldata