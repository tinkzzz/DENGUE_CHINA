function para = checkbound_para(para,paramax,paramin)
%check parameter bound
for i=1:size(para,1)
    para(i,para(i,:)<paramin(i))=paramin(i)*(1+0.2*rand(sum(para(i,:)<paramin(i)),1));
    para(i,para(i,:)>paramax(i))=paramax(i)*(1-0.2*rand(sum(para(i,:)>paramax(i)),1));
end

