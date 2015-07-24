function [results]= make_GLM_and_contrasts_from_inst_firing_mgs(y,closest_surf_phi,tphi,noise_type,mean_center,std_scale)

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

X1=double([ones(size(y,1),1) phys_targ/2]);

contrast.name={'Target (main)';};
% crazy_mat1=inv(X1'*X1);
 results.GLM(1).ces_std=zeros(size(X1,2),size(y,2));
 results.GLM(1).t=zeros(size(X1,2),size(y,2));
 results.GLM(1).Fsig=zeros(size(X1,2),size(y,2));
 results.GLM(1).ces=zeros(size(X1,2),size(y,2));
 results.GLM(1).dof=zeros(size(X1,2),size(y,2));
 rx=rank(X1);
 for i=1:size(y,2)
     disp(['GLM(1) ' num2str(i)])
     temp_y=y(:,i);
     tt=trials_with_good_stim_location & ~isnan(temp_y);
     if rank(X1(tt,:))<rx || rank(X1(tt,:))==length(X1(tt,:)) 
         continue
     end

         
        [b,dev,stats] = glmfit(X1(tt,:),temp_y(tt),noise_type,'link','identity','constant','off','estdisp','off');
        results.GLM(1).ces_std(:,i)=stats.se/std_scale;
        results.GLM(1).t(:,i)=stats.t;
        results.GLM(1).Fsig(:,i)=stats.p;
        results.GLM(1).ces(:,i)=stats.beta/std_scale;
        results.GLM(1).dof(:,i)=stats.dfe;

     
 end
 results.GLM(1).X=X1;
 results.GLM(1).contrast.name=contrast.name;
 