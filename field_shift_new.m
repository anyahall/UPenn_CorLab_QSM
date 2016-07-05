


R = 0.015; % m
N = 128;
[x,y]=meshgrid(linspace(-10*R,10*R,N));
B0 = 3; % [Tesla]
%chi_b=0.273*4*pi; %deox blood
chi_w = -8e-6; % susceptibility of water
r = sqrt(x.^2+y.^2);
r = r(:);
dx = 10*R/N;
dy = 10*R/N;

xv=x(:);
yv=y(:);

%create susc. distribution
chi_dist=[];
for i=1:length(r)
    if abs(r(i))<=R
        chi_dist(i)=chi_w;
    else
        chi_dist(i)=0;
    end
end

%calculate field shift -  long integral form of convolution
B = zeros(N*N,1);
for k=1:numel(xv)
    dB=zeros(N*N,1);
    if(mod(k,100)==0)
        disp(['Percent Done = ' num2str(k/N^2*100)]);
    end
    for j=1:numel(xv)
        if(k~=j)
            %mag_rrp = sqrt(x(j)^2+y(j)^2);
            a = [xv(j)-xv(k) yv(j)-yv(k)];  
            mag_rrp = sqrt(a(1)^2+a(2)^2);
            b = [0 1];  %B0 vector
            theta=findAngle(a,b);
            dB(j)=chi_dist(j)*(3*(cos(theta)^2)-1)/mag_rrp^3;
            %Bd(j)=Bd(j)*dx*dy;
        end
    end
    B(k)=B0*dx*dy*sum(dB);
end


%% analytical solution (Haacke textbook)

%dB=[]; dF=[]; yes=[]; no=[];
for j=1:numel(r)
    a=[0 1];
    b=[x(j) y(j)];
    th=findAngle(a,b);
   [Ba(j),Fa(j)]=sim_B(R,r(j),th);
end

%%

%B=B*N^2*dx*dy;
B=reshape(B,[N N]);

%%
gamma=42.57e6;
F1 = gamma*B;
F2 = gamma*Ba;
F2 = reshape(F2,[N N]);
imshow([F1 F2 F2-F1],[]); colormap jet

% %figure
% imshow(dF,[]); colormap jet; colorbar
% title('B-field shift [Hz]')