function [holoSet,u0,obj2D,dz] = MSHoloSim(X,Y,Z,D,opts)
% Engine for simulating a set of lensless digital in-line holographic 
% microscopy (DIHM) holograms (holoSet) for a group of microbeads with
% given 4D (space + time) locations and diameters with respect to multiple 
% scattering
%
% Input parameters:
%   X(m,f),Y(m,f),Z(m,f) - microbeads locations in 3D volume
%       m - microbead number
%       f - frame number
%       e.g., X = [1:3;3:-1:1]; Y = [1,1,1;2,2,2]; Z = [1,1,1;1:3]; will
%       generate a set of 3 holograms, each containing 2 beads: first that
%       travels from left to right and second that travels from right to
%       left and increases its hight in each next hologram
%       Beads will only be generated if their X and Y values will be in
%       range dx:dx*imSize, where dx = pixSize/mag. There is no limit in
%       possible Z positions (can be negative), however, high range of Z 
%       values may cause to significantly slower generation time (see 
%       opts.red)
%   D(m) - vectors containing microbeads diameters (um)
%   opts - generation options
%       opts.n = [n1,n0];
%           n1 - object refractive index (default 1.002)
%           n0 - background refractive index (default 1)
%       opts.A - absorbtion level (%) - how many % of the light is absorbed
%           by the microbead of highest possible diameter (default 50)
%       opts.mag - system magnification (default 13)
%       opts.pixSize - camera pixel size (um) (default 5.5)
%       opts.dist - propagation distance (mm) - distance between object and
%           hologram central plane (where Z(m,f) = 0) (default 0.63)
%       opts.lambda - lightsource wavelength (um) (default 0.66)
%       opts.SNR - hologram signal-to-noise ratio (default 20)
%       opts.imSize - hologram size ([x] or [y,x]; pix) (default 500)
%       opts.b_sig - sigma parameter of the hologram Gaussian background
%           (the higher is the sigma, the more darker are hologram edges).
%           (default 0.4)
%       opts.red - generated object volume may be very large. When z
%           dimension is larger than max(y,x)/2 reduce it to be equal
%           max(y,x)/opts.ZF ?  red = 'y' - yes, red = 'n' - no, red = 'a'
%           - ask me when this will happen (default opts.red = 'y')
%       opts.ZF - sampling of simulated volume in Z dim (object volume has 
%           dimensions: imSize x imSize x imSize/ZF) (default opts.ZF = 2)
%
% Auxiliary drawings:
%                          dist
%    _________ <------------------------->
%   |  object |                           |
%   |* volume |                           |
%   |       * |                           | detector/hologram
%   |  *      |                           | plane
%   |     *   |                           |
%   |_________|                           |
%  Z = 10   Z = 0                      dist = 0
%
%                           dist
% _____________________ <---------
%                      |
%  object   beads --> / \
%  volume     |       \ /
%             v        |
%                      |
%            / \       |
%            \ /       |
%            <->     Z = 0
%             D
%
% Output parameters:
%   holoSet - generated set of holograms
%   u0 - set of ground truth optical fields at the hologram plane (padded
%       to 3*imSize size)
%   obj2D - object optical field 2D representation
%   dz = object sampling in z direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland
%
% Last modified - 09.06.2022
% 
% See the https://github.com/MRogalski96/DarkTrack for more info
% 
% Cite as:
% [1] Mikołaj Rogalski, Jose Angel Picazo-Bueno, Julianna Winnik, Piotr 
% Zdańkowski, Vicente Micó, Maciej Trusiak. "DarkTrack: a path across the 
% dark-field for holographic 4D particle tracking under Gabor regime." 
% 2021. Submitted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with the input
if ~isfield(opts,'n'); n = [1.002, 1]; else; n = opts.n; end
if ~isfield(opts,'A'); A = 50; else; A = opts.A; end
if ~isfield(opts,'mag'); mag = 13; else; mag = opts.mag; end
if ~isfield(opts,'pixSize'); pixSize = 5.5; else; pixSize = opts.pixSize; end
if ~isfield(opts,'dist'); dist = 0.63; else; dist = opts.dist; end
if ~isfield(opts,'lambda'); lambda = 0.66; else; lambda = opts.lambda; end
if ~isfield(opts,'SNR'); SNR = 20; else; SNR = opts.SNR; end
if ~isfield(opts,'imSize'); imSize = 500; else; imSize = opts.imSize; end
if ~isfield(opts,'b_sig'); b_sig = 0.4; else; b_sig = opts.b_sig; end
if ~isfield(opts,'red'); red = 'y'; else; red = opts.red; end
if ~isfield(opts,'ZF'); ZF = 2; else; ZF = opts.ZF; end

if length(imSize) == 1; imSize = [imSize,imSize]; end

%% Initialization
% Sampling parameters
dx = pixSize/mag;
dy = dx;
dz = dx/5;

% Number of objects and number of frames
[NoO,NoF] = size(X);

% Amount of light that is transported through single microbead layer
Ab = (1-A/100)^(1/(max(D)/dz-1));

% Object volume Z dimenzion
Z_dim = max(Z(:))-min(Z(:))+max(D(:));


% Ensure that the sampling in z direction is not to dense to save the
% computational time
if Z_dim/dz > max(imSize)/ZF && (red == 'y' || red == 'a')
    if red == 'a'
        warning(['The size of object volume is very large: ', ...
            num2str(imSize(2)),'x',num2str(imSize(1)),'x',...
            num2str(round(Z_dim/dz)),'. Reduce it to: ',...
            num2str(imSize(2)),'x',num2str(imSize(1)),'x',...
            num2str(round(max(imSize)/ZF)),' ?'])
        red = input(' 1 - Yes / 0 - No = ');
        if red == 1; red = 'y'; end
    end
    if red == 'y'
        dz = (Z_dim+1)/max(imSize)*ZF;
    end
end

% Initialize outputs
holoSet = zeros([imSize,NoF]);
u0 = zeros([3*imSize,NoF]);
obj2D = zeros([imSize,NoF]);

minZ = min(Z(:)); maxZ = max(Z(:));
xo=dx:dx:imSize(2)*dx;
yo=dx:dx:imSize(1)*dx;
zo=minZ-max(D)/2:dz:maxZ+max(D)/2;
% Object volume coordinates (in um)
[XX,YY,ZZ] = meshgrid(xo,yo,zo);

for ff = 1:NoF
    %% 3D obj generation
    obj3D = zeros(size(XX));
    for tt = 1:NoO
        R = sqrt((XX-X(tt,ff)).^2+(YY-Y(tt,ff)).^2+(ZZ-Z(tt,ff)).^2);
        obj3D_TMP = R<D(tt)/2;
        obj3D = obj3D + obj3D_TMP;
    end
    obj3D = flip(obj3D,3);
    obj3D(obj3D>1) = 1;
    
    % object 3D refractive index distribution
    obj3Dn = obj3D;
    obj3Dn(obj3D==1) = n(1);
    obj3Dn(obj3D==0) = n(2);
    
    % object 3D amplitude distribution
    if Ab ~=0
        obj3Da = obj3D*Ab;
        obj3Da(obj3Da==0) = 1;
    else
        obj3Da = 1-obj3D;
    end
    
    % object 2D optical field representation
    tmp = sum(obj3D,3);
    tmp2 = ones(imSize);
    for tt = 1:size(obj3Da,3)
        tmp2 = tmp2.*obj3Da(:,:,tt);
    end
    obj2D(:,:,ff) = tmp2.*exp(-1i*tmp*(n(1)-n(2))*2*pi/lambda(1)*dz);
    
    %% Hologram generation
    % Input optical field
    uill = ones(3*imSize);
    % Pad obj3Dn and obj3Da to avoid near image edges errors
    obj3DnP = padarray(obj3Dn,imSize,n(2));
    obj3DaP = padarray(obj3Da,imSize,1);
    
    % Propagate the optical fiel through the object volume with BPM method
    u2dB = propagate_BPM_2d(obj3DnP,obj3DaP,uill,lambda,dx,dy,dz,n(2));

    % Final propagation to the detector/hologram plane
    u2d = propagate_AS_2d(u2dB,dist*1000+zo(1),n(2),lambda,dx);
    
    % Crop previously padded regions
    u2df = u2d(imSize(1)+1:2*imSize(1),imSize(2)+1:2*imSize(2));
    
    % Generate the hologram from the optical field
    holo0 = abs(u2df).^2;
    AA = (max(holo0(:))-min(holo0(:))); % Signal amplitude
    N = AA/SNR; % Amount of noise
    noise = randn(imSize)/0.798*N; % Noise component
    bckr = gausswin(imSize(1),b_sig)*gausswin(imSize(2),b_sig)'; % Background component
    
    % Outputs
    holoSet(:,:,ff) = bckr.*(1 + 1*holo0) + noise; % Generated hologram
    u0(:,:,ff) = u2d; % Ground truth optical field at hologram plane
    
    disp(['Simulated: ',num2str(ff),'/',num2str(NoF),' frames'])
end
% memory
end

%% Auxiliary functions
function u = propagate_BPM_2d(n3d, obj3Da, ui, lambda, dx, dy, dz, n0)
% Beam propagation method for 3D amplitude and refractive index
% distribution
% 
% Input parameters:
%   n3d - 3D refractive index distribution;
%   obj3Da - 3D amplitude distribution
%   z - coordinate in the direction of propagation
%   ui -  illumination beam
%   dx, dy, dz - sampling interval
%   n0 - background refractive index
% 
% Output parameters:
%   u - the optical field distribution after propagation through the sample
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Julianna Winnik, 09.02.2021
% adapted to 2D: Mikołaj Rogalski, 17.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0=2*pi/lambda;
[Ny, Nx, Nz] = size(n3d);
uo=zeros(Ny,Nx,Nz);

n3d=n3d-n0;

dfx=1/Nx/dx;
dfy=1/Ny/dy;
fx = (-Nx/2:Nx/2-1)*dfx;
fy = (-Ny/2:Ny/2-1)*dfy;
[Fx,Fy] = meshgrid(fx,fy);
F = sqrt(Fx.^2+Fy.^2);

%calculation kernels for bpm method
kernelfq = ifftshift(exp(-1i*2*dz*pi*F.^2./( sqrt((n0/lambda)^2 - F.^2) + n0/lambda ) ));

u=ui;
for iz=1:Nz
    uo(:,:,iz)=u;
    dn=n3d(:,:,iz);
    ab = obj3Da(:,:,iz);
    
    kernelsp = exp(1i*dz*k0*dn);

    %P: first FFT
    u = ifft2(fft2(u.*ab).*kernelfq); 
 
    %Q: lens 
    u = u.*kernelsp; 
 
end
% memory
end


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