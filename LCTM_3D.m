function [E,ix]=LCTM_3D(f,ns,N)
[IX,BB]=knnsearch(f,f,'k',ns);
B=BB';
ix=IX';

for i=1:N
    Cord(i,1)=sum(f(ix(1:ns,i),1))/ns;
    Cord(i,2)=sum(f(ix(1:ns,i),2))/ns;
    Cord(i,3)=sum(f(ix(1:ns,i),3))/ns;
    Rx=f(ix(1:ns,i),1)-Cord(i,1);
    Ry=f(ix(1:ns,i),2)-Cord(i,2);
    Rz=f(ix(1:ns,i),3)-Cord(i,3);
    [wtheta(:,:,i),wrho(:,:,i),wz(:,:,i)] = cart2pol(Rx,Ry,Rz);
end

m1=1;
m2=2;
DM1(1:ns,1,1:N)=1;
DM1(1:ns,2,1:N)=wz;
for ii=1:m1
    DM1(1:ns,2+ii,1:N)=cosh(ii.*wz).*besselj(0,ii.*wrho);
    DM1(1:ns,2+m1+ii,1:N)=sinh(ii.*wz).*besselj(0,ii.*wrho);
    DM1(1:ns,2+2*m1+ii,1:N)=cos(ii.*wz).*besseli(0,ii.*wrho);
    DM1(1:ns,2+3*m1+ii,1:N)=sin(ii.*wz).*besseli(0,ii.*wrho);
end

for kk=1:m1
    for jj=1:m2
        DM1(1:ns,2+4*m1+(kk-1)*m2+jj,1:N)=cos(jj.*wtheta).*cosh(kk.*wz).*besselj(jj,kk.*wrho);
        DM1(1:ns,2+4*m1+m1*m2+(kk-1)*m2+jj,1:N)=sin(jj.*wtheta).*sinh(kk.*wz).*besselj(jj,kk.*wrho);
        DM1(1:ns,2+4*m1+2*m1*m2+(kk-1)*m2+jj,1:N)=cos(jj.*wtheta).*sinh(kk.*wz).*besselj(jj,kk.*wrho);
        DM1(1:ns,2+4*m1+3*m1*m2+(kk-1)*m2+jj,1:N)=sin(jj.*wtheta).*cosh(kk.*wz).*besselj(jj,kk.*wrho);
        DM1(1:ns,2+4*m1+4*m1*m2+(kk-1)*m2+jj,1:N)=cos(jj.*wtheta).*cos(kk.*wz).*besseli(jj,kk.*wrho);
        DM1(1:ns,2+4*m1+5*m1*m2+(kk-1)*m2+jj,1:N)=sin(jj.*wtheta).*sin(kk.*wz).*besseli(jj,kk.*wrho);
        DM1(1:ns,2+4*m1+6*m1*m2+(kk-1)*m2+jj,1:N)=cos(jj.*wtheta).*sin(kk.*wz).*besseli(jj,kk.*wrho);
        DM1(1:ns,2+4*m1+7*m1*m2+(kk-1)*m2+jj,1:N)=sin(jj.*wtheta).*cos(kk.*wz).*besseli(jj,kk.*wrho);  
    end
end

for ll=1:m2
    DM1(1:ns,2+4*m1+8*m1*m2+ll,1:N)=wrho.^ll.*cos(ll.*wtheta);
    DM1(1:ns,2+4*m1+8*m1*m2+m2+ll,1:N)=wrho.^ll.*sin(ll.*wtheta);
    DM1(1:ns,2+4*m1+8*m1*m2+2*m2+ll,1:N)=wz.*wrho.^ll.*cos(ll.*wtheta);
    DM1(1:ns,2+4*m1+8*m1*m2+3*m2+ll,1:N)=wz.*wrho.^ll.*sin(ll.*wtheta); 
end

%% Direct way
for i=1:N
    E(:,:,i)=DM1(1,:,i)*pinv(DM1(:,:,i),10^-10);
end
clear DM1

%% MLS way
% ww=1-6*(B./B(ns,:)).^2+8*(B./B(ns,:)).^3-3*(B./B(ns,:)).^4;
% w=kron(ones(size(DM1,2),1),ww);
% w=reshape(w,ns,size(DM1,2),N);
% ZZ=DM1.*w.^2;
% Y1=permute(ZZ,[2,1,3]);Z1=Y1;%Y1(:,1,:)=-sum(Y1,2);
% % Z1=[sum(ZZ1*DM1)];
% clear ZZ
% 
% for i=1:N
%     E(:,:,i)=DM1(1,:,i)*pinv(Z1(:,:,i)*DM1(:,:,i),10^-10)*Y1(:,:,i);
% end
% clear Z1 Y1 DM1
end