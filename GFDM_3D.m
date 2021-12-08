function [Ed,ix]=GFDM_3D(f,ns,N)
[IX,BB]=knnsearch(f,f,'k',ns);
B=BB';
ix=IX';

for i=1:N
    h1(:,i)=f(ix(:,i),1)-f(i,1);
    h2(:,i)=f(ix(:,i),2)-f(i,2);
    h3(:,i)=f(ix(:,i),3)-f(i,3);
end

w=ones(ns,length(f))-6*(B./B(ns,:)).^2+8*(B./B(ns,:)).^3-3*(B./B(ns,:)).^4;
h1=reshape(h1,ns,1,length(f));h2=reshape(h2,ns,1,length(f));h3=reshape(h3,ns,1,length(f));w=reshape(w,ns,1,length(f));w(1,:,:)=0;
ZZ=[h1,h2,h3,h1.^2/2,h2.^2/2,h3.^2/2,h1.*h2,h1.*h3,h2.*h3].*w.^2;
Y=permute(ZZ,[2,1,3]);Y(:,1,:)=-sum(Y,2);
Z=[sum(ZZ.*h1);sum(ZZ.*h2);sum(ZZ.*h3);sum(ZZ.*h1.^2/2);sum(ZZ.*h2.^2/2);sum(ZZ.*h3.^2/2);sum(ZZ.*h1.*h2);sum(ZZ.*h1.*h3);sum(ZZ.*h2.*h3);];
clear ZZ

for i=1:N
    E(:,:,i)=pinv(Z(:,:,i))*Y(:,:,i);
end
    Ed=E(4,:,:)+E(5,:,:)+E(6,:,:);
end