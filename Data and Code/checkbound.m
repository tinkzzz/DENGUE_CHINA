function [Sm,Sh,Im,Ihr,Ihu]=checkbound(Sm,Sh,Im,Ihr,Ihu)
%check variable bound
for k=1:size(Sm,2)
    Sm(Sm(:,k)<0,k)=0; Im(Im(:,k)<0,k)=0; 
    Sh(Sh(:,k)<0,k)=0; Ihr(Ihr(:,k)<0,k)=0; Ihu(Ihu(:,k)<0,k)=0;
end
