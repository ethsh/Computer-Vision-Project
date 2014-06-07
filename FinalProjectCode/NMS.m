function [Inms] =NMS( I,n,T )
%n - size of window for mlocal maximum
%T - threshold for local maximum

 Inms=zeros(size(I)); %and this will hold the final result of the edges after the filtering
 Width=size(I,2);
 Height=size(I,1);
 
 
 R_left=I;
 ind=find(R_left>=T);
 not_zero=R_left(ind);
 R_left_not_zero=zeros(length(not_zero),3);
 R_left_not_zero(:,1)=not_zero;
 [y,x]=ind2sub(size(R_left),ind);
 R_left_not_zero(:,2:3)=[y,x];
 R_left_not_zero=sortrows(R_left_not_zero,-1);

 index=1;
for i=1:size(R_left_not_zero,1)
     row= R_left_not_zero( i,2);
     col= R_left_not_zero( i,3);
     if R_left(row,col)~=0
            Inms(row,col)=index;
            index=index+1;
            R_left(max((row-n),1):min((row+n),Height),max(1,col-n):min(Width,col+n))=0;
     end

end
