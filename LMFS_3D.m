function [E,Ex,Ey,Ez,ix]=LMFS_3D(f,ns,N,dim)
[IX,BB]=knnsearch(f,f,'k',ns);
B=BB';
ix=IX';

B_m=max(B,[],1);
rt=10;
QR=(rt*B_m).';
QQ1=generateB(0.5*ns,1,'sphere').';
for i=1:N
   QQ(:,:,i)=QQ1.*QR(i,1)+ones(0.5*ns,1)*f(i,:);
   RA1(:,:,1:dim+1,i)=myDistanceMatrix(f(ix(:,i),:),QQ(:,:,i),dim);
end
RA(:,:,1:N)=RA1(:,:,dim+1,:);
RAx(:,:,1:N)=RA1(:,:,1,:);
RAy(:,:,1:N)=RA1(:,:,2,:);
RAz(:,:,1:N)=RA1(:,:,3,:);

DM1(:,:,1:N)=1./(4.*pi.*RA);
DM1x(:,:,1:N)=-RAx./(4.*pi.*RA.^3);
DM1y(:,:,1:N)=-RAy./(4.*pi.*RA.^3);
DM1z(:,:,1:N)=-RAz./(4.*pi.*RA.^3);

%% Direct way
for i=1:N
    paa=pinv(DM1(:,:,i),10^-8);
    E(:,:,i)=DM1(1,:,i)*paa;
    Ex(:,:,i)=DM1x(1,:,i)*paa;
    Ey(:,:,i)=DM1y(1,:,i)*paa;
    Ez(:,:,i)=DM1z(1,:,i)*paa;
end

%% MLS way
% ww=1-6*(B./B_m(1,:)).^2+8*(B./B_m(1,:)).^3-3*(B./B_m(1,:)).^4;
% w=kron(ones(0.5*ns,1),ww);
% w=reshape(w,ns,0.5*ns,N);
% ZZ=DM1.*w.^2;
% Y1=permute(ZZ,[2,1,3]);Z1=Y1;
% clear ZZ
% 
% for i=1:N
%     E(:,:,i)=DM1(1,:,i)*pinv(Z1(:,:,i)*DM1(:,:,i),10^-8)*Y1(:,:,i);
%     Ex(:,:,i)=DM1x(1,:,i)*pinv(Z1(:,:,i)*DM1(:,:,i),10^-8)*Y1(:,:,i);
%     Ey(:,:,i)=DM1y(1,:,i)*pinv(Z1(:,:,i)*DM1(:,:,i),10^-8)*Y1(:,:,i);
%     Ez(:,:,i)=DM1z(1,:,i)*pinv(Z1(:,:,i)*DM1(:,:,i),10^-8)*Y1(:,:,i);
% end

end