function [Ed,ix]=LRBFCM_3D(f,ns,N,c)
[Ix,~]=knnsearch(f,f,'k',ns);
ix=Ix';

for i=1:N
    DM(:,:,1:4,i)=myDistanceMatrix(f(ix(:,i),:),f(ix(:,i),:),3);
    phi(:,:,i)=sqrt(DM(:,:,4,i).^2+c^2);
    phix(:,:,i)=DM(:,:,1,i)./sqrt(DM(:,:,4,i).^2+c^2);
    phiy(:,:,i)=DM(:,:,2,i)./sqrt(DM(:,:,4,i).^2+c^2);
    phiz(:,:,i)=DM(:,:,3,i)./sqrt(DM(:,:,4,i).^2+c^2);
    phixx(:,:,i)=(DM(:,:,2,i).^2+DM(:,:,3,i).^2+c^2)./sqrt(DM(:,:,4,i).^2+c^2).^3;
    phiyy(:,:,i)=(DM(:,:,1,i).^2+DM(:,:,3,i).^2+c^2)./sqrt(DM(:,:,4,i).^2+c^2).^3;
    phizz(:,:,i)=(DM(:,:,1,i).^2+DM(:,:,2,i).^2+c^2)./sqrt(DM(:,:,4,i).^2+c^2).^3;
end
    
for i=1:N
    E(1,:,i)=phi(:,:,i)\phix(:,1,i);
    E(2,:,i)=phi(:,:,i)\phiy(:,1,i);
    E(3,:,i)=phi(:,:,i)\phiz(:,1,i);
    E(4,:,i)=phi(:,:,i)\phixx(:,1,i);
    E(5,:,i)=phi(:,:,i)\phiyy(:,1,i);
    E(6,:,i)=phi(:,:,i)\phizz(:,1,i);
end
    Ed=E(4,:,:)+E(5,:,:)+E(6,:,:);
end