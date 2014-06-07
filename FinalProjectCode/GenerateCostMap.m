function [ g,I] = GenerateCostMap( I,k,c,p1,p2,flag)
%This function generates the cost map to be used in the functional
%I - image
%k - size of dilation structure element
%p1,2 - percentile of threshold
%c - negative cost

global dispflag;

%I=I.^2;
I=abs(I);
T=prctile(I(1:end),p1);
I(I<T)=0;
I=(I-mean(mean(I)))./(std(std(I))+eps);


if dispflag==1
    figure;
    imagesc(I); title('Image before dilation');
end
se = strel('disk',k);
if flag==1
    Id = imdilate(I,se);
elseif flag==2
    Id = imclose(I,se);
end
if dispflag==1
    figure;
    imagesc(Id); title('Image after dilation');
end



g=(Id./max(max(Id)));
%g=Id/prctile(Id(1:end),90);


T=prctile(g(1:end),p2);
g(g<=T)=c;

if dispflag==1    
    figure;
    imagesc(g); title('g after dilation');
end



end

