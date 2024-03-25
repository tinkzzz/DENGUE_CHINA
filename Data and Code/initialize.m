function [Sm,Sh,Im,Ihr,Ihu,Seedc]=initialize(num_loc,nl,Nij,num_ens,obs_season,Cloc,mosnum)
%compute initial conditions for all subpopulations
%
num_mp=size(nl,1);
Sm=round((mosnum.*Cloc)*ones(1,num_ens));
Sh=Nij.*ones(1,num_ens);
Im=zeros(num_loc,num_ens);
Ihr=zeros(num_mp,num_ens);
Ihu=zeros(num_mp,num_ens);

num_times=size(obs_season,2);
Seedc=zeros(num_loc,num_times);

reported=sum(obs_season,2);     %total incidence
reported(:,2)=(1:size(reported,1));
reported=sortrows(reported,-1);

seedlocs=reported(reported(:,1)>0,2);     %cases >0

for l=1:size(seedlocs)
    seedloc=seedlocs(l);
    temp=obs_season(seedloc,:); 
    %find peaks
    below=temp<=1;
    startpoints=[];
    cnt=0;
    last=1;
    for i=1:length(below)
        if below(i)==1
            cnt=cnt+1;
        else
            if cnt>=7
                startpoints=[startpoints,last]; 
                last=i;
            end
            cnt=0;
        end
    end
    for i=1:length(startpoints)
        index=find(temp>=5);
        index=index(index>=startpoints(i));
        if ~isempty(index)
            T0=index(1);%first reporting with at least 5 cases
            c=sum(temp(T0:min(T0+4,num_times)));
            Seedc(seedloc,max(1,T0-7))=c;
        end
    end
end
