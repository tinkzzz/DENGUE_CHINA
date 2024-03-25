function [para,paramax,paramin,alphamap,betamap]=initializepara_eakf(obs_season,num_ens)
%initialize ensemble members
%para: Z,D,mu,theta,alpha1,alpha2,...,alpha337,beta1,...,beta337
Dlow=5;Dup=7;%infectious period
mulow=1;muup=1;%relative transmissibility
thetalow=0.25;thetaup=0.45;%movement factor
alphalow=0.25;alphaup=0.35;%reporting rate
betalow=0.3;betaup=0.5;%Relative transmission rate 

num_loc=size(obs_season,1);
U=0.25;
death=2/30;
%define alpha
alphamap=5+(1:num_loc)';
%define beta
betamap=5+num_loc+(1:num_loc)';

load weight.mat     %initial parameter weights
%D,theta,alpha1,alpha2,...,alpha337,beta1,...,beta337
paramin=[U;death;Dlow;mulow;thetalow;ones(num_loc,1)*alphalow;ones(num_loc,1)*betalow]';
paramax=[U;death;Dup;muup;thetaup;ones(num_loc,1)*alphaup;ones(num_loc,1)*betaup]';

parausemin=[U;death;Dlow;mulow;thetalow;ones(num_loc,1)*alphalow;ones(num_loc,1)*betalow]';
parausemax=[U;death;Dup;muup;thetaup;ones(num_loc,1)*alphaup;ones(num_loc,1)*betalow+0.4*(betaup-betalow)]';
para=lhsu(parausemin,parausemax,num_ens);

para=para';
paramin=paramin';
paramax=paramax';

