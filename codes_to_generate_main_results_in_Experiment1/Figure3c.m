% do pca dimension reduction using bootstrapping
%
%
% based on mack2020  NC paper
%
%
% prepare your 4d.nii after single pattern estimation
% including activation patterns (trial*voxel) under HSC/LSC condition


% The toolbox should be installed:
%    1. cosmoMvpa;

clear;clc;
subjects = {'sub-02';'sub-03';'sub-04'};
masks = {'.nii'};
study_path='I:\FLXX1\Dimension\pattern';
roi_path='I:\FLXX1\mask';
msk = masks{1};

alldata = [];
counterp = 0;

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
        d_HSC=ds_HSC.samples';

        data_LSC=fullfile(sub_path,'glm_T_stats_LSC.nii');
        ds_LSC=cosmo_fmri_dataset(data_LSC,'mask',mask_fn);
        d_LSC=ds_LSC.samples';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        true_length = min(size(d_HSC,2),size(d_LSC,2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  true_length ~= size(d_HSC,2)
            
            allrd = 0;
            for x = 1:1000
                da_HSC = d_HSC(:,randperm(size(d_HSC,2),true_length)); %randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n.
                [~,~,~,~,explained,~]=pca(da_HSC);
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
            [~,~,~,~,explained,~]=pca(d_HSC);
            total_explained_hsc=0;

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


        %%
        if  true_length ~= size(d_LSC,2)

            allrd = 0;
            for x = 1:1000
                da_LSC = d_LSC(:,randperm(size(d_LSC,2),true_length)); %randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n.
                [~,~,~,~,explained,~]=pca(da_LSC);
                total_explained_lsc=0;

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

            [~,~,~,~,explained,~]=pca(d_LSC);
            total_explained_lsc=0;

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

    [h,p,ci,stats]=ttest(RDhsc,RDlsc)
    counterp = counterp + 1;
    allp(counterp,:) = p;
    alldata = [alldata,RDhsc,RDlsc];
end

%save sametrialn_70_80_IFG.mat allp alldata
