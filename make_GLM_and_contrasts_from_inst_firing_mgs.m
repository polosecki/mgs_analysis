function [results]= make_GLM_and_contrasts_from_inst_firing_mgs(y,closest_surf_phi,tphi)

%% Define variables

angles=unique(tphi);
RF_angles_distance=abs(angles-closest_surf_phi);
RF_angles_distance(RF_angles_distance>180)=RF_angles_distance(RF_angles_distance>180)-180;
[closest_targ_phi,min_index]=min(RF_angles_distance); 


locations=circshift(unique(tphi),-min_index+1);%RF_centered; for computation purposes
non_shifted_locations=unique(tphi); %For nomenclature purposes

locations=[locations(1); locations(5)];
non_shifted_locations=[non_shifted_locations(1); non_shifted_locations(5)];

trials_with_good_stim_location=ismember(tphi,locations);
phys_targ=2*ismember(tphi,locations(1))-1;




%% RUN GLMs for effect sizes
%To obtain effect sizes am interested in:
%-Main effect of target postion

X1=double([ones(size(y,1),1) phys_targ]);



contrast.matrix={[0 1]/2; %main effect of target
                 };
             
contrast.name={'Target (main)';};
 crazy_mat1=inv(X1'*X1);
 results.GLM(1).ces_std=zeros(length(contrast.matrix),size(y,2));
 results.GLM(1).ces_covar=cell(length(contrast.matrix),size(y,2));
 results.GLM(1).F=zeros(length(contrast.matrix),size(y,2));
 results.GLM(1).Fsig=zeros(length(contrast.matrix),size(y,2));
 results.GLM(1).ces=zeros(length(contrast.matrix),size(y,2));
 results.GLM(1).ces_vector=cell(length(contrast.matrix),size(y,2));
 results.GLM(1).dof=zeros(length(contrast.matrix),size(y,2));
 for i=1:size(y,2)
     disp(['GLM(1) ' num2str(i)])
     temp_y=y(:,i);
     if rank(X1(trials_with_good_stim_location & ~isnan(temp_y),:))<rank(X1) || rank(X1(trials_with_good_stim_location & ~isnan(temp_y),:))==length(X1(trials_with_good_stim_location & ~isnan(temp_y),:)) 
         continue
     end
     [beta, rvar, ~, ~] = fast_glmfit(temp_y(trials_with_good_stim_location & ~isnan(temp_y)),X1(trials_with_good_stim_location & ~isnan(temp_y),:));

     for c=1:length(contrast.matrix)
         [F, Fsig, ces, edof] = fast_fratio(beta,X1(trials_with_good_stim_location & ~isnan(temp_y),:),rvar,contrast.matrix{c});
         ces_covar=rvar*contrast.matrix{c}*crazy_mat1*contrast.matrix{c}';
         ces_std=sqrt(mean(diag(ces_covar)))/sqrt(length(diag(ces_covar)));
         results.GLM(1).ces_std(c,i)=ces_std;
         results.GLM(1).ces_covar{c,i}=ces_covar;
         results.GLM(1).F(c,i)=F;
         results.GLM(1).Fsig(c,i)=Fsig;
         results.GLM(1).ces(c,i)=mean(ces);
         results.GLM(1).ces_vector{c,i}=ces;
         results.GLM(1).dof(c,i)=edof;
     end
     
     
     results.GLM(1).beta(:,i)=beta;
     results.GLM(1).rvar(1,i)=rvar;
     
 end
 results.GLM(1).X=X1;
 results.GLM(1).contrast.name=contrast.name;
 results.GLM(1).contrast.matrix=contrast.matrix;
 