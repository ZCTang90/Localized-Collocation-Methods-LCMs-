%---------------------------------------------------------------------
% %产生boundary points
% 2009. 4. 25  Li Ming
% domian ('square','round','flower','sphere','doublesphere')
% nI: the No. of points
%-----------------------------------------------------------------------
function [bdpt,normalvec]=generateB(nB,part,domain)
switch lower(domain)
    case {'square'}     
        x=linspace(0,1,nB);
        y=linspace(0,1,nB);
        switch part
            case 1/4
               Bx=[x];
               By=[zeros(1,length(x))];
               normalx=zeros(1,length(x));
               normaly=-1*ones(1,length(y));
            case 1/2
                y=linspace(0,0.5,nB/2);
                Bx=[zeros(1,length(y)),x,10*ones(1,length(y))];
                By=[y,zeros(1,length(x)),y];
                normalx=[-1*ones(1,length(y)),zeros(1,length(x)),ones(1,length(y))];
                normaly=[zeros(1,length(y)),-1*ones(1,length(x)),zeros(1,length(y))];
            case 1.5
                Bx=[x,10*ones(1,length(y)-1)];
                By=[zeros(1,length(x)),y(2:nB)];
                normalx=[zeros(1,length(x)),10*ones(1,length(y)-1)];
                normaly=[ones(1,length(x)),zeros(1,length(y)-1)];
            case 3/4
                Bx=[zeros(1,length(y)-1),x,10*ones(1,length(y)-1)];
                By=[y(2:nB),zeros(1,length(x)),y(2:nB)];
                normalx=[-1*ones(1,length(y)-1),zeros(1,length(x)),ones(1,length(y)-1)];
                normaly=[zeros(1,length(y)-1),ones(1,length(x)),zeros(1,length(y)-1)];
            case 1
                Bx=[zeros(1,length(y)-1),x,ones(1,length(y)-1),x(2:nB-1)];
                By=[y(2:nB),zeros(1,length(x)),y(2:nB),ones(1,length(x)-2)];
                normalx=[-1*ones(1,length(y)-1),zeros(1,length(x)),ones(1,length(y)-1),zeros(1,nB-2)];
                normaly=[zeros(1,length(y)-1),-1*ones(1,length(x)),zeros(1,length(y)-1),ones(1,nB-2)];
            otherwise
               disp('There are 1/4,1/2,3/4,1 to be chosen in square domain')
        end
        bdpt=[Bx;By];
        normalvec=[normalx;normaly];
    case {'round'}    
        theta=linspace(1/(nB-1),2*pi,nB);
        Bx=cos(theta);
        By=sin(theta);
        bdpt=[Bx(1:nB*part);By(1:nB*part)];
        normalvec=[Bx(1:nB*part);By(1:nB*part)];
    %case {'peanut'}
    case {'flower'}
        theta=linspace(1/(nB-1),part*2*pi,nB*part);
        rho=(cos(3*theta)+sqrt(2-(sin(3*theta)).^2)).^(1/3);
        Bx=rho.*cos(theta);
        By=rho.*sin(theta);
        bdpt=[Bx;By];
        normalvec=[Bx;By]; %这个要修正
    %case {'star'}
    %case {'amoeba'}    
    case {'sphere'}
        dlong = pi*(3-sqrt(5));
        long = 0;
        dz = sqrt(1)*2.0/nB;
        z = sqrt(1) - dz/2;
        for k = 0 : part*nB-1
            r = sqrt(1-z.^2);
            bdpt(2,k+1)=cos(long)*r;
            bdpt(3,k+1)=sin(long)*r;
            bdpt(1,k+1)=z;
            z = z - dz;
            long = long + dlong;
        end
        normalvec=bdpt;
    case {'doublesphere'}
        dlong = pi*(3-sqrt(5));
        long = 0;
        dz = sqrt(1)*1.75*2.0/nB;
        z = 1.75*sqrt(1) - dz/2;
        if part>0.5
            for k = 0 : 0.5*nB-1
                r = sqrt(1-(z-0.75).^2);
                bdpt(2,k+1)=cos(long)*r;
                bdpt(3,k+1)=sin(long)*r;
                bdpt(1,k+1)=z;
                normalvec(2,k+1)=bdpt(2,k+1);%liwen
                normalvec(3,k+1)= bdpt(3,k+1);%liwen
                normalvec(1,k+1)=bdpt(1,k+1)-0.75;%liwen
                z = z - dz;
                long = long + dlong;
            end
            for k = 0.5*nB : part*nB-1
                r = sqrt(1-(z+0.75).^2);
                bdpt(2,k+1)=cos(long)*r;
                bdpt(3,k+1)=sin(long)*r;
                bdpt(1,k+1)=z;
                normalvec(2,k+1)=bdpt(2,k+1);%liwen
                normalvec(3,k+1)= bdpt(3,k+1);%liwen
                normalvec(1,k+1)=bdpt(1,k+1)+0.75;%liwen
                z = z - dz;
                long = long + dlong;
            end
           
        else
             for k = 0 : part*nB-1
                r = sqrt(1-(z-0.75).^2);
                bdpt(2,k+1)=cos(long)*r;
                bdpt(3,k+1)=sin(long)*r;
                bdpt(1,k+1)=z;
                normalvec(2,k+1)=bdpt(2,k+1);%liwen
                normalvec(3,k+1)= bdpt(3,k+1);%liwen
                normalvec(1,k+1)=bdpt(1,k+1)-0.75;%liwen
                z = z - dz;
                long = long + dlong;
             end
        end
        bdpt=bdpt';%liwen
        normalvec=normalvec';%liwen
        for k=1:length(normalvec) %liwen
        normalvec(k,:)= normalvec(k,:)./norm(normalvec(k,:)); %liwen
        end %liwen
    otherwise
        disp('The domain model must be one of square, round, flower,sphere,doublesphere')
end
