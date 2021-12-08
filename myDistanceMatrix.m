%% Distance Matrix
function DM = myDistanceMatrix(dsites, ctrs,dim)
[M,s] = size(dsites); N = length(ctrs); 
DM = zeros(M,N,dim+1);
for i=1:s
    [dr,cc]=ndgrid(dsites(:,i),ctrs(:,i));
    DM(:,:,i) = (dr-cc);
    DM(:,:,dim+1) = DM(:,:,dim+1) +(dr-cc).^2;
end
DM(:,:,dim+1)=sqrt(DM(:,:,dim+1));
end
