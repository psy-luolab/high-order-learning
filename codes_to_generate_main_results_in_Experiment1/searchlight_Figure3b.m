% do searchlight using sheng2022 method
% quick but not accurate

% The toolbox should be installed:
%    1. cosmoMvpa;

% Ensure that my_var_measure.m is located in the current folder

clear;clc;
subjects = {'sub-02';'sub-03'};
study_path='I:\FLXX1\Dimension\pattern';
voxel_count=100;
roi_path='L:\FLXX_1\roi';
masks = {'mask.nii'};
msk = masks{1};
%%
measure=@searchlight_Func;
%%
for s = 1:length(subjects)
    sub = subjects{s};
    
    sub_path=fullfile(study_path,sub);
    output_path=fullfile(study_path,sub);

    roi=fullfile(roi_path,sub);

    mask_fn=fullfile(roi,msk);

    ds_all=fullfile(sub_path,'glm_T_stats_HSC.nii');

    ds_all=cosmo_fmri_dataset(ds_all,'mask',mask_fn);

    nh=cosmo_spherical_neighborhood(ds_all,'count',voxel_count);

    my_map=cosmo_searchlight(ds_all,nh,measure);

    output_fn_all=fullfile(output_path,'rdvar_hsc_70.nii');

    cosmo_map2fmri(my_map,output_fn_all);
end

%%
for s = 1:length(subjects)
    sub = subjects{s};
    
    sub_path=fullfile(study_path,sub);
    output_path=fullfile(study_path,sub);

    roi=fullfile(roi_path,sub);

    mask_fn=fullfile(roi,msk);

    ds_all=fullfile(sub_path,'glm_T_stats_LSC.nii');

    ds_all=cosmo_fmri_dataset(ds_all,'mask',mask_fn);

    nh=cosmo_spherical_neighborhood(ds_all,'count',voxel_count);

    my_map=cosmo_searchlight(ds_all,nh,measure);

    output_fn_all=fullfile(output_path,'rdvar_lsc_70.nii');

    cosmo_map2fmri(my_map,output_fn_all);
end
