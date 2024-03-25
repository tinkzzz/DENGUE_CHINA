function [Sm,Sh,Im,Ihr,Ihu]=seeding(Sm,Sh,Im,Ihr,Ihu,part,C,Seedc,t)
%set initial infection for all subpopulations
num_loc=size(Seedc,1);
num_ens=size(Sh,2);
for l=1:num_loc
    seedloc=l;
    seedIm=round(rand(1,num_ens)*12*Seedc(l,t));
    if t>=8
        seedIhr=round(rand(1,num_ens)*Seedc(l,t));
        seedIhu=round(rand(1,num_ens)*Seedc(l,t));
    end
    if Seedc(l,t)>0
        pop=sum(C(part(seedloc):part(seedloc+1)-1));
        for i=part(seedloc):part(seedloc+1)-1
            Ihr(i,:)=round(seedIhr*C(i)/pop);
            Ihu(i,:)=round(seedIhu*C(i)/pop);
            Sh(i,:)=Sh(i,:)-Ihr(i,:)-Ihu(i,:);
        end
        Im(l,:)=round(seedIm);
        Sm(l,:)=Sm(l,:)-Im(l,:);
    end
end