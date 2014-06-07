function [ F,gradF1,gradF2,gradF3] = CostFunc(I,Ix,Iy,m,n,R,x,y )
%Calculates the cost function and the gradient used in the FindBall script at the point (m,n,R);
%I=flipud(I);
 [theta,r] = cart2pol(x-m,y-n);
 Ipolar=[reshape(theta,numel(I),1),reshape(r,numel(I),1),reshape(I,numel(I),1),reshape(Ix,numel(I),1),reshape(Iy,numel(I),1)];
 
 %**************************************
 %{
 Ipolar=sortrows(Ipolar,1);
 phi0=Ipolar(1,1);
 ShiftSize=find(Ipolar(:,1)<=(phi0+phi),1,'last')-1;
 Ipolar(:,1)=circshift(Ipolar(:,1),-ShiftSize);
 
 
 [xrot,yrot]=pol2cart(Ipolar(:,1),Ipolar(:,2));
 xrot=xrot+m; yrot=yrot+n;
 idxrot=sub2ind(size(I),min(size(I,1),max(1,round(yrot))),min(size(I,2),max(1,round(xrot))));
 Itest=zeros(size(I));
 
 Itest(idxrot)=Ipolar(:,3);
% figure(1);imagesc(Itest+I);drawnow;
 %}
 %**************************************
 
 
 
 Ipolar=sortrows(Ipolar,2); % An array : [ theta , r , I(theta,r)]
%{
 figure;
 subplot(1,3,1);imagesc(Ipolar(:,1));colorbar;
 subplot(1,3,2);imagesc(Ipolar(:,2));colorbar;
 subplot(1,3,3);imagesc(Ipolar(:,3));colorbar;
 %}

 
 
 deltaR=0.1;
 Ridx_first=find(Ipolar(:,2)<=(R-deltaR),1,'last');
 Ridx_last=find(Ipolar(:,2)<=(R+deltaR),1,'last');
 F=-sum(Ipolar(1:Ridx_last,3).*Ipolar(1:Ridx_last,2));
%-sum(Ipolar(1:Ridx_last,6));
 
%Compute gradF

            %projection of gradI on the direction of rotation
            %{
            Im=[Ipolar(1:Ridx_last,4) ,Ipolar(1:Ridx_last,5)]*[cos(phi);-sin(phi)]; %dI/dm
            In=[Ipolar(1:Ridx_last,4) ,Ipolar(1:Ridx_last,5)]*[sin(phi);cos(phi)]; %dI/dn
            %}
            
            gradF1=sum(Ipolar(1:Ridx_last,2).*Ipolar(1:Ridx_last,4));
            gradF2=sum(Ipolar(1:Ridx_last,2).*Ipolar(1:Ridx_last,5));
            gradF3=-R.*sum(Ipolar(Ridx_first:Ridx_last,3));
            %{
            gradF1=sum(Ipolar(1:Ridx_last,2).*Im);
            gradF2=sum(Ipolar(1:Ridx_last,2).*In);
            gradF3=-R.*sum(Ipolar(Ridx_first:Ridx_last,3));
            %}


end

