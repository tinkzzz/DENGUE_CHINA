load commutingdata.mat     %Commuting data
load population.mat        %city population

Nij=C(:,1);      %population movements on the first day

load dailyincidence.mat    %the number of peported cases in each city on each day
load parameter.mat         %identified model parameters
load mosmun.mat  %weighting of initial mosquito population in each city

m_rate=2/30;     %fixed mortality rate of mosquito populations
birth_rate= Birth_rate+m_rate;        %mosquito population birth rate
birth_rate=permute(birth_rate,[2,1]);  

%%%%%%%%%%%%%%%%%%%%Initial parameter range%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
num_ens=300;     %number of ensemble members
num_times=364;     %total length of data
T=num_times;     %number of days to keep track of reported cases
num_loc=size(part,1)-1;     %number of locations
num_mp=size(nl,1);     %number of metapopulation

%set OEV, observation error variance
OEV_case=zeros(size(obs_season));
for l=1:num_loc
    for t=1:num_times
        OEV_case(l,t)=max(1e-4,obs_season(l,t)^2/100);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize model
[Sm,Sh,Im,Ihr,Ihu,Seedc]=initialize(num_loc,nl,Nij,num_ens,obs_season,Cloc,mosnum);
obs_temp=zeros(num_loc,num_ens,num_times);%records of reported cases
%initialize parameters
[para,paramax,paramin,alphamap,betamap]=initializepara_eakf(obs_season,num_ens);

parastd=std(para,0,2);     %get ensemble spread of parameters
%%%%%%%%%%%%%%%%%%%%%
dailyIm_post_rec=zeros(num_loc,num_ens,T);     %posterior mosquito infection
dailyIhr_post_rec=zeros(num_loc,num_ens,T);     %posterior human infection
dailyIhu_post_rec=zeros(num_loc,num_ens,T);     %posterior human infection
dailyIm_prior_rec=zeros(num_loc,num_ens,T);     %prior mosquito infection
dailyIhr_prior_rec=zeros(num_loc,num_ens,T);     %prior human infection
dailyIhu_prior_rec=zeros(num_loc,num_ens,T);     %prior human infection

dailyIhr_post_rec1=zeros(num_loc,num_ens,T);     %posterior human infection
dailyIhu_post_rec1=zeros(num_loc,num_ens,T);     %posterior human infection


%initialize poseteriors (percentage of population)
Sm_post=zeros(num_loc,num_times,num_ens);
Sh_post=zeros(num_loc,num_times,num_ens);
Im_post=zeros(num_loc,num_times,num_ens);
Ihr_post=zeros(num_loc,num_times,num_ens);
Ihu_post=zeros(num_loc,num_times,num_ens);
%initialize cumulative reported and unreported infections
cumu_dailyIm_post=zeros(num_loc,num_ens);
cumu_dailyIhr_post=zeros(num_mp,num_ens);
cumu_dailyIhu_post=zeros(num_mp,num_ens);
%%%%%%%%%%%%%%%%%%%%%%%
lambda=1.2;     %inflation in EAKF to avoid divergence
num_para=size(para,1);     %number of parameters
para_post=zeros(num_para,num_ens,num_times);     %posterior parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% start EAKF of network model

for t=140:num_times     %start from May 20
    t
    tic
    %seeding
    if t<=size(Seedc,2)
        [Sm,Sh,Im,Ihr,Ihu]=seeding(Sm,Sh,Im,Ihr,Ihu,part,C(:,1),Seedc,t);
    end
    para=checkbound_para(para,paramax,paramin);
    %%%%%%%%%%%%%%%%%%
    %integrate forward one step
    dailyIm_prior=zeros(num_loc,num_ens);    
    dailyIhr_prior=zeros(num_mp,num_ens);
    dailyIhu_prior=zeros(num_mp,num_ens);

    for k=1:num_ens     %run for each ensemble member
        [Sm(:,k),Sh(:,k),Im(:,k),Ihr(:,k),Ihu(:,k),dailyIm_temp,dailyIhr_temp,dailyIhu_temp]=model_eakf(nl,part,C(:,1),Cave(:,1), ...
            Sm(:,k),Sh(:,k),Im(:,k),Ihr(:,k),Ihu(:,k),para(:,k),tao_season(:,t),betahtom_season(:,t),betamtoh_season(:,t),birth_rate(:,t),betamap,alphamap);
        dailyIhr_prior(:,k)=dailyIhr_temp;
        dailyIhu_prior(:,k)=dailyIhu_temp;
        dailyIm_prior(:,k)=dailyIm_temp;
    end
        dailyIhr_post=dailyIhr_prior;
        dailyIhu_post=dailyIhu_prior;
        dailyIm_post=dailyIm_prior;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for l=1:num_loc
        for k=1:num_ens
            for j=part(l):part(l+1)-1
                obs_temp(l,k,t)=obs_temp(l,k,t)+dailyIhr_prior(j,k);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    obs_ens=obs_temp(:,:,t);      %observation at t1, prior
    %loop through local observations
    for l=1:num_loc
        %%%%%%%%%%%%%%%%%%%case
        if (sum(obs_season(:,t))>0)
            %Get the variance of the ensemble
            obs_var = OEV_case(l,t);
            prior_var = var(obs_ens(l,:));
            post_var = prior_var*obs_var/(prior_var+obs_var);
            if prior_var==0       %if degenerate
                post_var=1e-3;
                prior_var=1e-3;
            end
            prior_mean = mean(obs_ens(l,:));
            post_mean = post_var*(prior_mean/prior_var + obs_season(l,t)/obs_var);
            %%%% Compute alpha and adjust distribution to conform to posterior moments
            alpha = (obs_var/(obs_var+prior_var)).^0.5;
            dy = post_mean + alpha*(obs_ens(l,:)-prior_mean)-obs_ens(l,:);
            %Loop over each state variable (connected to location l)
            %Im
            temp=Im(l,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            Im(l,:)=Im(l,:)+dx;
            %dailyIm_prior
            temp=dailyIm_post(l,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            dailyIm_post(l,:)=dailyIm_post(l,:)+dx;            
            %adjust related metapopulation
            neighbors=part(l):part(l+1)-1;    %metapopulation live in l
            for h=1:length(neighbors)
                j=neighbors(h);
                %Ihr
                temp=Ihr(j,:);
                A=cov(temp,obs_ens(l,:));
                rr=A(2,1)/prior_var;
                dx=rr*dy;
                Ihr(j,:)=Ihr(j,:)+dx;
                %Ihu
                temp=Ihu(j,:);
                A=cov(temp,obs_ens(l,:));
                rr=A(2,1)/prior_var;
                dx=rr*dy;
                Ihu(j,:)=Ihu(j,:)+dx;
                %dailyIhr
                temp=dailyIhr_post(j,:);
                A=cov(temp,obs_ens(l,:));
                rr=A(2,1)/prior_var;
                dx=rr*dy;
                dailyIhr_post(j,:)=round(max(dailyIhr_post(j,:)+dx,0));
                %dailyIhu
                temp=dailyIhu_post(j,:);
                A=cov(temp,obs_ens(l,:));
                rr=A(2,1)/prior_var;
                dx=rr*dy;
                dailyIhu_post(j,:)=round(max(dailyIhu_post(j,:)+dx,0));
            end
            %adjust beta
            temp=para(betamap(l),:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            para(betamap(l),:)=para(betamap(l),:)+dx;
            %inflation
            if t<=250
                if std(para(betamap(l),:))<parastd(betamap(l))
                    para(betamap(l),:)=mean(para(betamap(l),:),2)*ones(1,num_ens)+lambda*(para(betamap(l),:)-mean(para(betamap(l),:),2)*ones(1,num_ens));
                end
            end
        end
    end
    para=checkbound_para(para,paramax,paramin);

    %update posterior Ihr and Ihu
    cumu_dailyIhr_post=cumu_dailyIhr_post+dailyIhr_post;
    cumu_dailyIhu_post=cumu_dailyIhu_post+dailyIhu_post;
    cumu_dailyIm_post=cumu_dailyIm_post+dailyIm_post;
    

    %%%%%%%%%%%%%%%%update Sh Sm
    Sh=C(:,1)*ones(1,num_ens)-cumu_dailyIhr_post-cumu_dailyIhu_post;
    if sum(sum(dailyIm_post)) ~= 0
        Sm=Sm-dailyIm_post+dailyIm_prior;
    end

    %%%%%%%%%%%%%%%%
    [Sm,Sh,Im,Ihr,Ihu]=checkbound(Sm,Sh,Im,Ihr,Ihu);
    %%%%%%%%save statevariables
    for i=1:num_loc
        for j=1:num_ens
            Sh_post(i,t,j)=sum(Sh(part(i):part(i+1)-1,j));
            Ihr_post(i,t,j)=sum(Ihr(part(i):part(i+1)-1,j));
            Ihu_post(i,t,j)= sum(Ihu(part(i):part(i+1)-1,j));
            dailyIhr_post_rec(i,j,t)=sum(dailyIhr_post(part(i):part(i+1)-1,j));
            dailyIhu_post_rec(i,j,t)=sum(dailyIhu_post(part(i):part(i+1)-1,j));
            dailyIhr_prior_rec(i,j,t)=sum(dailyIhr_prior(part(i):part(i+1)-1,j));
            dailyIhu_prior_rec(i,j,t)=sum(dailyIhu_prior(part(i):part(i+1)-1,j));
            dailyIhr_post_rec1(i,j,t)=sum(dailyIhr_post((part(i)+1):part(i+1)-1,j));
            dailyIhu_post_rec1(i,j,t)=sum(dailyIhu_post((part(i)+1):part(i+1)-1,j));
        end
    end
    for i=1:num_loc
        for j=1:num_ens
            Sm_post(i,t,j)=Sm(i,j);
            Im_post(i,t,j)=Im(i,j);
            dailyIm_post_rec(i,j,t)=dailyIm_post(i,j);
            dailyIm_prior_rec(i,j,t)=dailyIm_prior(i,j);
        end
    end    
    para_post(:,:,t)=para;
    toc
end

save inference.mat para_post Sm_post Im_post Sh_post Ihr_post Ihu_post... 
 dailyIm_post_rec dailyIhr_post_rec dailyIhu_post_rec dailyIm_prior_rec...
 dailyIhr_prior_rec dailyIhu_prior_rec dailyIhr_post_rec1 dailyIhu_post_rec1;





