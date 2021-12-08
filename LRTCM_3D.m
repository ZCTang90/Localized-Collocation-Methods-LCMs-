function [E,Ex,Ey,Ez,ix]=LRTCM_3D(f,ns,N,dim,c)
[IX,BB]=knnsearch(f,f,'k',ns);
B=BB';
ix=IX';

for i=1:N
    DM1(:,:,1:dim+1,i)=myDistanceMatrix(f(ix(:,i),:), f(ix(:,i),:),dim);
end
DM1u(:,:,1:N)=exp(-c.*(DM1(:,:,1,:).^2-DM1(:,:,2,:).^2)).*cos(2.*c.*DM1(:,:,1,:).*DM1(:,:,2,:))+...
    exp(-c.*(DM1(:,:,2,:).^2-DM1(:,:,3,:).^2)).*cos(2.*c.*DM1(:,:,2,:).*DM1(:,:,3,:))+...
    exp(-c.*(DM1(:,:,3,:).^2-DM1(:,:,1,:).^2)).*cos(2.*c.*DM1(:,:,3,:).*DM1(:,:,1,:));
DM1x(:,:,1:N)=2.*c.*DM1(:,:,1,:).*cos(2.*c.*DM1(:,:,1,:).*DM1(:,:,3,:)).*exp(c.*(DM1(:,:,1,:).^2 - DM1(:,:,3,:).^2)) - 2.*c.*DM1(:,:,1,:).*cos(2.*c.*DM1(:,:,1,:).*DM1(:,:,2,:)).*exp(-c.*(DM1(:,:,1,:).^2 - DM1(:,:,2,:).^2)) - 2.*c.*DM1(:,:,2,:).*sin(2.*c.*DM1(:,:,1,:).*DM1(:,:,2,:)).*exp(-c.*(DM1(:,:,1,:).^2 - DM1(:,:,2,:).^2)) - 2.*c.*DM1(:,:,3,:).*sin(2.*c.*DM1(:,:,1,:).*DM1(:,:,3,:)).*exp(c.*(DM1(:,:,1,:).^2 - DM1(:,:,3,:).^2));
DM1y(:,:,1:N)=2.*c.*DM1(:,:,2,:).*cos(2.*c.*DM1(:,:,1,:).*DM1(:,:,2,:)).*exp(-c.*(DM1(:,:,1,:).^2 - DM1(:,:,2,:).^2)) - 2.*c.*DM1(:,:,2,:).*cos(2.*c.*DM1(:,:,2,:).*DM1(:,:,3,:)).*exp(-c.*(DM1(:,:,2,:).^2 - DM1(:,:,3,:).^2)) - 2.*c.*DM1(:,:,1,:).*sin(2.*c.*DM1(:,:,1,:).*DM1(:,:,2,:)).*exp(-c.*(DM1(:,:,1,:).^2 - DM1(:,:,2,:).^2)) - 2.*c.*DM1(:,:,3,:).*sin(2.*c.*DM1(:,:,2,:).*DM1(:,:,3,:)).*exp(-c.*(DM1(:,:,2,:).^2 - DM1(:,:,3,:).^2));
DM1z(:,:,1:N)=2.*c.*DM1(:,:,3,:).*cos(2.*c.*DM1(:,:,2,:).*DM1(:,:,3,:)).*exp(-c.*(DM1(:,:,2,:).^2 - DM1(:,:,3,:).^2)) - 2.*c.*DM1(:,:,3,:).*cos(2.*c.*DM1(:,:,1,:).*DM1(:,:,3,:)).*exp(c.*(DM1(:,:,1,:).^2 - DM1(:,:,3,:).^2)) - 2.*c.*DM1(:,:,1,:).*sin(2.*c.*DM1(:,:,1,:).*DM1(:,:,3,:)).*exp(c.*(DM1(:,:,1,:).^2 - DM1(:,:,3,:).^2)) - 2.*c.*DM1(:,:,2,:).*sin(2.*c.*DM1(:,:,2,:).*DM1(:,:,3,:)).*exp(-c.*(DM1(:,:,2,:).^2 - DM1(:,:,3,:).^2));

for i=1:ns
DM1x(i,i,:)=0;
DM1y(i,i,:)=0;
DM1z(i,i,:)=0;
end

%% Direct way
for i=1:N
    E(:,:,i)=DM1u(1,:,i)*pinv(DM1u(:,:,i),10^-10);
    Ex(:,:,i)=DM1x(1,:,i)*pinv(DM1u(:,:,i),10^-10);
    Ey(:,:,i)=DM1y(1,:,i)*pinv(DM1u(:,:,i),10^-10);
    Ez(:,:,i)=DM1z(1,:,i)*pinv(DM1u(:,:,i),10^-10);
end
clear Z1 Y1 DM1 DM1x DM1y DM1z

%% MLS way
% ww=1-6*(B./B(ns,:)).^2+8*(B./B(ns,:)).^3-3*(B./B(ns,:)).^4;
% w=kron(ones(ns,1),ww);
% w=reshape(w,ns,ns,N);
% ZZ=DM1u.*w.^2;
% Y1=permute(ZZ,[2,1,3]);Z1=Y1;
% clear ZZ
% 
% for i=1:N
%     E(:,:,i)=DM1u(1,:,i)*pinv(Z1(:,:,i)*DM1u(:,:,i),10^-10)*Y1(:,:,i);
%     Ex(:,:,i)=DM1x(1,:,i)*pinv(Z1(:,:,i)*DM1u(:,:,i),10^-10)*Y1(:,:,i);
%     Ey(:,:,i)=DM1y(1,:,i)*pinv(Z1(:,:,i)*DM1u(:,:,i),10^-10)*Y1(:,:,i);
%     Ez(:,:,i)=DM1z(1,:,i)*pinv(Z1(:,:,i)*DM1u(:,:,i),10^-10)*Y1(:,:,i);
% end
% clear Z1 Y1 DM1
end