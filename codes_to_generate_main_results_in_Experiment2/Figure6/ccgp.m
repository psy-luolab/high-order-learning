% ccgp(cross condition generalization performance) using svm
% need cosmoMvpa and svm


clear
clc
subject_ids={'sub-01';'sub-02';'sub-03';'sub-04';'sub-05';'sub-06';'sub-07';'sub-08';'sub-09';'sub-10';'sub-11';'sub-12';'sub-13';'sub-14';'sub-15';'sub-16';'sub-17';'sub-18';'sub-19';'sub-20';'sub-21';'sub-22';'sub-23';'sub-24';'sub-25';'sub-26';'sub-27';'sub-28';'sub-29';'sub-30';'sub-31';'sub-32';'sub-33';'sub-34';'sub-36';'sub-37';'sub-38';'sub-39';'sub-40'};

nsubjects=numel(subject_ids);

s2_path='H:\GJXX_2_reanalysis\N\First_level\S2\rsa';
s1_path='H:\GJXX_2_reanalysis\N\First_level\S1\rsa';

data_path='H:\GJXX_2_reanalysis\N\First_level\S1';

out_path='H:\GJXX_2_reanalysis\decoding\result_svm';

for i_subj=1:39
    subject_id=subject_ids{i_subj};
    sub_path=fullfile(data_path,subject_id);

    mask_fn=fullfile(sub_path,'mask.nii');
    sub_s2_path=fullfile(s2_path,subject_id);
    sub_s1_path=fullfile(s1_path,subject_id);

    %%
    output_path=fullfile(out_path,subject_id);

    if ~exist(output_path)
        mkdir(output_path);
    end
    %%

    data_HSC=fullfile(sub_s2_path,'glm_HSC_NA.nii');
    ds_hsc_s2=cosmo_fmri_dataset(data_HSC,'mask',mask_fn);
    %
    data_hsc_fn=fullfile(sub_s1_path,'glm_HSC_NA.nii');
    ds_hsc_s1=cosmo_fmri_dataset(data_hsc_fn,'mask',mask_fn);

    data_HSC=fullfile(sub_s2_path,'glm_LSC_NA.nii');
    ds_lsc_s2=cosmo_fmri_dataset(data_HSC,'mask',mask_fn);

    data_hsc_fn=fullfile(sub_s1_path,'glm_LSC_NA.nii');
    ds_lsc_s1=cosmo_fmri_dataset(data_hsc_fn,'mask',mask_fn);

    all_ds_train=cosmo_stack({ds_hsc_s1,ds_lsc_s1});
    all_ds_test=cosmo_stack({ds_hsc_s2,ds_lsc_s2});

    ds=cosmo_stack({all_ds_train,all_ds_test});

    %%
    measure_args = struct();
    measure = @searchlight_func;

    measure_args.s1_H=size(ds_hsc_s1.samples,1);

    measure_args.s1_L=size(ds_lsc_s1.samples,1);

    measure_args.s2_H=size(ds_hsc_s2.samples,1);

    measure_args.s2_L=size(ds_lsc_s2.samples,1);
    %%
    voxel_count=100;
    nbrhood=cosmo_spherical_neighborhood(ds,'count',voxel_count);
    p = cosmo_searchlight(ds,nbrhood,measure,measure_args);

    cosmo_map2fmri(p, ...
        fullfile(output_path,'svmmap.nii'));
end
