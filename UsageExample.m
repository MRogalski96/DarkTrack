%% Info
% This code demonstrates the exemplary usage of "MSHoloSim" engine for
% simulating a set of lensless digital in-line holographic microscopy 
% (DIHM) holograms and the usage of "DarkTrack" algorithm for
% reconstructing the objects locations from the set of DIHM holograms.
% 
% List of contents
% 1. Simulating hologram with MSHoloSim and reconstructing with AS method
% 2. Simulating a set of holograms by MSHoloSim and using DarkTrack to 
%    recover objects 4D positions
% 3. Using DarkTrack to recover simulated data (spiral microbead) object 4D 
%    positions
% 4. Using DarkTrack to recover real data objects (human sperm) 4D
%    positions
%
% Exemplary datasets, that are used in 2., 3. and 4., may be downloaded at:
% https://drive.google.com/drive/folders/1UNtZ3IeEX5ms_Vx85b0S685D-uipLe_l?usp=sharing
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland
%
% Last modified: 09.06.2022
% 
% See the https://github.com/MRogalski96/DarkTrack for more info
% 
% Cite as:
% [1] Mikołaj Rogalski, Jose Angel Picazo-Bueno, Julianna Winnik, Piotr 
% Zdańkowski, Vicente Micó, Maciej Trusiak. "DarkTrack: a path across the 
% dark-field for holographic 4D particle tracking under Gabor regime." 
% 2021. Submitted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Simulating hologram with MSHoloSim and reconstructing with AS method
clear
close all
clc

% System parameters (see MSHoloSim.m for description)
opts.n = [1.002, 1];
opts.A = 50; % (%)
opts.mag = 13;
opts.pixSize = 5.5; % (um)
opts.dist = 0.63; % (mm)
opts.lambda = 0.66; %(um)
opts.SNR = 20;
opts.imSize = 500; % (pix)
opts.b_sig = 0.4;
opts.ZF = 10; % in case of out of memory error increase the ZF parameter

% Simulating X,Y,Z positions of five microbeads (one in image center and
% four in corners) at differend distances from detector (Z+opts.dist*1000)
dx = opts.pixSize/opts.mag;
X = [100; 100; 250; 400; 400]*dx; % (um)
Y = [100; 400; 250; 100; 400]*dx; % (um)
Z = [50; -50; 0; -50; 50]; % (um)

% Microbeads diameters
D = [10; 8; 6; 12; 7];

% Simulate hologram
[holo,u0] = MSHoloSim(X,Y,Z,D,opts);

% Displaying simulation results
S = opts.imSize;
xx = (1:S)*dx; yy = (1:S)*dx;
figure; subplot(1,3,1); imagesc(xx,yy,holo); axis image;
title('Simulated hologram'); xlabel \mum; ylabel \mum
subplot(1,3,2); imagesc(xx,yy,abs(u0(S+1:2*S,S+1:2*S))); axis image;
title('Optical field at hologram plane - amplitude'); xlabel \mum;
ylabel \mum
subplot(1,3,3); imagesc(xx,yy,angle(u0(S+1:2*S,S+1:2*S))); axis image;
title('Optical field at hologram plane - phase'); xlabel \mum; ylabel \mum
colormap gray;

% Reconstruction with AS method
out1 = propagate_AS_2d(holo,-opts.dist*1000-Z(1),opts.n(2),opts.lambda,dx);
out1_wf = propagate_AS_2d(u0,-opts.dist*1000-Z(1),opts.n(2),opts.lambda,dx);
out2 = propagate_AS_2d(holo,-opts.dist*1000-Z(2),opts.n(2),opts.lambda,dx);
out2_wf = propagate_AS_2d(u0,-opts.dist*1000-Z(2),opts.n(2),opts.lambda,dx);
out3 = propagate_AS_2d(holo,-opts.dist*1000-Z(3),opts.n(2),opts.lambda,dx);
out3_wf = propagate_AS_2d(u0,-opts.dist*1000-Z(3),opts.n(2),opts.lambda,dx);

% Displaying reconstruction results
figure;
subplot(2,3,1); imagesc(xx,yy,abs(out1)); axis image;
title(['Amplitude at Z = ',num2str(Z(1)),' \mum']); xlabel \mum; ylabel \mum
subplot(2,3,2); imagesc(xx,yy,abs(out3)); axis image;
title(['Amplitude at Z = ',num2str(Z(3)),' \mum']); xlabel \mum; ylabel \mum
subplot(2,3,3); imagesc(xx,yy,abs(out2)); axis image;
title(['Amplitude at Z = ',num2str(Z(2)),' \mum']); xlabel \mum; ylabel \mum
subplot(2,3,4); imagesc(xx,yy,abs(out1_wf(S+1:2*S,S+1:2*S))); axis image;
title(['Amplitude from optical field at Z = ',num2str(Z(1)),' \mum'])
xlabel \mum; ylabel \mum
subplot(2,3,5); imagesc(xx,yy,abs(out3_wf(S+1:2*S,S+1:2*S))); axis image;
title(['Amplitude from optical field at Z = ',num2str(Z(3)),' \mum'])
xlabel \mum; ylabel \mum
subplot(2,3,6); imagesc(xx,yy,abs(out2_wf(S+1:2*S,S+1:2*S))); axis image;
title(['Amplitude from optical field at Z = ',num2str(Z(2)),' \mum'])
xlabel \mum; ylabel \mum; colormap gray
%
figure;
subplot(2,3,1); imagesc(xx,yy,angle(out1)); axis image;
title(['Phase at Z = ',num2str(Z(1)),' \mum']); xlabel \mum; ylabel \mum
subplot(2,3,2); imagesc(xx,yy,angle(out3)); axis image;
title(['Phase at Z = ',num2str(Z(3)),' \mum']); xlabel \mum; ylabel \mum
subplot(2,3,3); imagesc(xx,yy,angle(out2)); axis image;
title(['Phase at Z = ',num2str(Z(2)),' \mum']); xlabel \mum; ylabel \mum
subplot(2,3,4); imagesc(xx,yy,angle(out1_wf(S+1:2*S,S+1:2*S))); axis image;
title(['Phase from optical field at Z = ',num2str(Z(1)),' \mum'])
xlabel \mum; ylabel \mum
subplot(2,3,5); imagesc(xx,yy,angle(out3_wf(S+1:2*S,S+1:2*S))); axis image;
title(['Phase from optical field at Z = ',num2str(Z(3)),' \mum'])
xlabel \mum; ylabel \mum
subplot(2,3,6); imagesc(xx,yy,angle(out2_wf(S+1:2*S,S+1:2*S))); axis image;
title(['Phase from optical field at Z = ',num2str(Z(2)),' \mum'])
xlabel \mum; ylabel \mum; colormap gray

%% 2. Simulating a set of holograms by MSHoloSim and using DarkTrack to 
% recover objects 4D positions
clear
close all
clc

% Load exemplary X,Y,Z,D,opts parameters consisting 10 microbeads in 50
% frames
load('./Data/ExemplaryXYZDparameters_forSimulations.mat')

% System parameters
opts.n = [1.002, 1];
opts.A = 50; % (%)
opts.mag = 13;
opts.pixSize = 5.5; % (um)
opts.dist = 0.63; % (mm)
opts.lambda = 0.66; %(um)
opts.SNR = 20;
opts.imSize = 500; % (pix)
opts.b_sig = 0.4;
opts.ZF = 10; % in case of out of memory error increase the ZF parameter

% Simulating holograms (it may last dozen minutes)
[holo,~] = MSHoloSim(X,Y,Z,D,opts);

% Show generated holograms
dx = opts.pixSize/opts.mag;
xx = (1:size(holo,2))*dx; yy = (1:size(holo,1))*dx;
figure;
NoH = size(holo,3);
for tt = 1:NoH
    imagesc(xx,yy,holo(:,:,tt));
    title(['Simulated hologram ', num2str(tt),'/',num2str(NoH)])
    xlabel \mum; ylabel \mum; colormap gray; axis image;
    pause(0.1)
end

% DarkTrack algorithm
opts.propRange = [min(Z(:))-30, max(Z(:))+30]; % (um)
opts.propStep = 1; % (um)
adv.showTmpRes = 2; % additionally show the EDOF results during reconstruction

% DarkTrack algorithm
[outX,outY,outZ,EDOF,CR] = DarkTrack(holo,opts,adv);
outZ = outZ + 2; % DarkTrack is finding a focus position slightly above the
% microbead center. Slightly shift outZ positions to compensate this for
% better comparison with ground truth

% Displaying results - EDOF vs CR
miEDOF = min(EDOF(:)); maEDOF = max(EDOF(:));
miCR = min(CR(:)); maCR = max(CR(:));
figure('units','normalized','outerposition',[0 0 1 1]);
for tt = 1:NoH
    subplot(1,2,1); imagesc(xx,yy,EDOF(:,:,tt),[miEDOF,maEDOF]); axis image; xlabel \mum;
    ylabel \mum; colormap gray
    title(['EDOF reconstruction ', num2str(tt),'/',num2str(NoH)])
    set(gca,'fontsize',20)
    subplot(1,2,2); imagesc(xx,yy,CR(:,:,tt),[miCR,maCR]); axis image; xlabel \mum;
    ylabel \mum; colormap gray
    set(gca,'fontsize',20)
    title(['CR reconstruction at the middle of propRange ', num2str(tt),'/',num2str(NoH)])
    pause(0.1)
   
end
% Displaying results - Found X,Y,Z,T positions vs ground truth positions
figure; plot3(X(:),Y(:),Z(:),'.k'); hold on;
C = [X(1,5),Y(1,5),Z(1,5)] ;   % center of circle
R = D(1)/2;    % Radius of circle
teta=0:0.01:2*pi;
x=C(1)+R*cos(teta);
z=C(3)+R*sin(teta); y2=C(2)+R*sin(teta);
y=C(2)+zeros(size(x)); z2=C(3)+zeros(size(x));
plot3([x,x],[y,y2],[z,z2],'-r','Linewidth',2)
plot3(outX',outY',outZ','.-'); grid on; hold off
axis equal
xlim([dx,dx*size(holo,2)]); ylim([dx,dx*size(holo,1)]);
zlim([opts.propRange(1),opts.propRange(2)]); xlabel('x [\mum]')
ylabel('y [\mum]'); zlabel('z [\mum]');
title('Microbeads paths of movement')
nms{1} = 'ground truth positions';
nms{2} = 'microbead size';
for ss = 3:12; nms{ss} = ['microbead ',num2str(ss-2)]; end
legend(nms);
set(gca, 'YDir','reverse')

%% 3. Using DarkTrack to recover simulated data (spiral microbead) object 4D positions
clear
clc
close all

% Load holograms and opts
load('./Data/Data from article/SimulatedData/Data_Video1.mat'); clear u0

% DarkTrack
adv.showTmpRes = 2; % additionally show the EDOF results during reconstruction

% DarkTrack algorithm
[outX,outY,outZ,EDOF,CR] = DarkTrack(holo,opts,adv);
outZ = outZ + 1; % DarkTrack is finding a focus position slightly above the
% microbead center. Slightly shift outZ positions to compensate this for
% better comparison with ground truth

% Displaying results - EDOF vs CR
NoH = size(holo,3);
dx = opts.pixSize/opts.mag;
xx = (1:size(holo,2))*dx; yy = (1:size(holo,1))*dx;
miEDOF = min(EDOF(:)); maEDOF = max(EDOF(:));
miCR = min(CR(:)); maCR = max(CR(:));
figure('units','normalized','outerposition',[0 0 1 1]);
for tt = 1:NoH
    subplot(1,2,1); imagesc(xx,yy,EDOF(:,:,tt),[miEDOF,maEDOF]); axis image; xlabel \mum;
    ylabel \mum; colormap gray
    title(['EDOF reconstruction ', num2str(tt),'/',num2str(NoH)])
    subplot(1,2,2); imagesc(xx,yy,CR(:,:,tt),[miCR,maCR]); axis image; xlabel \mum;
    ylabel \mum; colormap gray
    title(['CR reconstruction at middle of propRange ', num2str(tt),'/',num2str(NoH)])
    pause(0.1)
end

% Displaying results - Found X,Y,Z,T positions vs ground truth positions
figure; plot3(X(:),Y(:),Z(:),'.k'); hold on;
C = [X(1,5),Y(1,5),Z(1,5)] ;   % center of circle
R = D(1)/2;    % Radius of circle
teta=0:0.01:2*pi;
x=C(1)+R*cos(teta);
z=C(3)+R*sin(teta); y2=C(2)+R*sin(teta);
y=C(2)+zeros(size(x)); z2=C(3)+zeros(size(x));
plot3([x,x],[y,y2],[z,z2],'-r','Linewidth',2)
plot3(outX',outY',outZ','.-'); grid on; hold off
axis equal
xlim([dx,dx*size(holo,2)]); ylim([dx,dx*size(holo,1)]);
zlim([opts.propRange(1),opts.propRange(2)]); xlabel('x [\mum]')
ylabel('y [\mum]'); zlabel('z [\mum]');
title('Microbeads paths of movement')
nms{1} = 'ground truth positions';
nms{2} = 'microbead size';
for ss = 3:size(outX,1)+2; nms{ss} = ['microbead ',num2str(ss-2)]; end
legend(nms);
set(gca, 'YDir','reverse')
%% 4. Using DarkTrack to recover real data objects (human sperm) 4D positions
clear
clc
close all

% Load holograms
mov = VideoReader('./Data/Data from article/RealData/HumanSpermHolo.avi');
t = 0;
while hasFrame(mov)
    t = t+1;
    tmp = readFrame(mov);
    holo(:,:,t) = tmp(:,:,3);
end

% DarkTrack 
% Load opts. When reconstructing your data, most of the opts are easy to
% set basing on the knowledge of your system. Only opts.dist and
% opts.propRange may be challenging - to set them properly we recommend to
% use DarkFocus autofocusing algorithm: https://github.com/MRogalski96/DarkFocus
load('./Data/Data from article/RealData/HumanSpermParameters.mat');

adv.NoF = 50; % reconstruct only first 50 frames
adv.showTmpRes = 2; % additionally show the EDOF results during reconstruction

[outX,outY,outZ,EDOF,CR] = DarkTrack(holo,opts,adv);

% Displaying results - EDOF vs CR
dx = opts.pixSize/opts.mag;
xx = (1:size(holo,2))*dx; yy = (1:size(holo,1))*dx;
miEDOF = min(EDOF(:)); maEDOF = max(EDOF(:));
miCR = min(CR(:)); maCR = max(CR(:));
figure('units','normalized','outerposition',[0 0 1 1]);
NoH = size(EDOF,3);
for tt = 1:NoH
    subplot(1,2,1); imagesc(xx,yy,EDOF(:,:,tt),[miEDOF,maEDOF]); axis image; xlabel \mum;
    ylabel \mum; colormap gray
    title(['EDF reconstruction ', num2str(tt),'/',num2str(NoH)])
    subplot(1,2,2); imagesc(xx,yy,CR(:,:,tt),[miCR,maCR]); axis image; xlabel \mum;
    ylabel \mum; colormap gray
    title(['CR reconstruction at middle of propRange ', num2str(tt),'/',num2str(NoH)])
    pause(0.1)
end

% Displaying results - Found X,Y,Z,T positions
m = size(outX,1); if m>5; m=5; end
for tt = 0:m
    figure;
    plot3(outX',outY',outZ','.-'); grid on;
    axis equal
    xlim([dx,dx*size(holo,2)]); ylim([dx,dx*size(holo,1)]);
    zlim([opts.propRange(1),opts.propRange(2)]); xlabel('x [\mum]')
    ylabel('y [\mum]'); zlabel('z [\mum]');
    title('Human sperm paths of movement')
    set(gca, 'YDir','reverse')
    if tt > 0
        xlim([min(outX(tt,:)),max(outX(tt,:))])
        ylim([min(outY(tt,:)),max(outY(tt,:))])
        zlim([min(outZ(tt,:)),max(outZ(tt,:))])
        title('Single spermatozoid movement path')
    end
    
end



%% Auxiliary functions
function uout = propagate_AS_2d(uin,z,n0,lambda,dx)
% Angular spectrum propagation method

[Ny,Nx] = size(uin);
k = 2*pi/lambda;

dfx = 1/Nx/dx; fx = -Nx/2*dfx : dfx : (Nx/2-1)*dfx;
dfy = 1/Ny/dx; fy = -Ny/2*dfy : dfy : (Ny/2-1)*dfy;

if  z<0
    kernel = exp(-1i*k*z*sqrt(n0^2 - lambda^2*(ones(Ny,1)*(fx.^2)+(fy'.^2)*ones(1,Nx))));
    ftu = kernel.*fftshift(fft2(fftshift(conj(uin))));
    uout = conj(fftshift(ifft2(ifftshift(ftu))));
else
    kernel = exp(1i*k*z*sqrt(n0^2 - lambda^2*(ones(Ny,1)*(fx.^2)+(fy'.^2)*ones(1,Nx))));
    ftu = kernel.*fftshift(fft2(fftshift(uin)));
    uout = fftshift((ifft2(ifftshift(ftu))));
end
end