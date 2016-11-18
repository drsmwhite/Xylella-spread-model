function Xylella_Spread_Model_With_Control
% XYLELLA_SPREAD_MODEL_WITH_CONTROL Model code for simulating the spread of
% Xylella fastidiosa in Apulia with a dual buffer zone. Note that this code
% does not use the actual olive cover map as we do have a licence to share
% this. Instead, we have supplied a randomly generated map that gives the
% essence of the simulation. See manuscript for more details.
%
% Code created by Dr Steven White (18/11/16). Please email
% (smwhit@ceh.ac.uk) for further information.

bufferwidth1=2; %Length of first buffer
bufferwidth2=13; %Lenth of second buffer
eff=0.9; %Control efficiency in the buffer zone
a=0; %Relative carrying capacity in non olive plants
timerun=8; %How many years does the model run for?

% Run the model
spreadmap=spread_model(bufferwidth1,bufferwidth2,eff,a,timerun);

% Plot the results
load('olivegrowthprop.mat') % 1km^2 gridded landscape and olive cover
figure
for i=1:timerun
    pause(1)
    pic=spreadmap(:,:,i); % Set the population matrix for plotting
    pic(olivegrowthprop<0)=-0.01;
    imagesc(pic(144:260,224:310)); % Plot the zoomed map
    axis image; axis off;
    title(['{\it{Xylella fastidiosa}} spread, ',num2str(i),' years since initial outbreak'])
    map=[(0:0.01:1)', (1:-0.01:0)', zeros(101,1)];
    map(1,:)=[0,0,1];
    colormap(map);
    caxis([-0.01,1]);
    c=colorbar;
    ylabel(c,'Disease incidence')
    
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
end

end

function output = spread_model(bufferwidth,bufferwidth2,eff,a,timerun)
%The main spread model with buffer zone

load('olivegrowthprop.mat') % 1km^2 gridded landscape and olive cover

[N,M]=size(olivegrowthprop); % Obtain the map size

olivebuffer=-9999*ones(3*N,3*M); % create a buffer map
olivebuffer((N+1):2*N,(M+1):2*M) = olivegrowthprop; %put olive map in the middle

% Create a buffer zone by removing all trees in the zone
xp1=253; yp1=212; xp2=278; yp2=183;
A1=(yp1-yp2)/(xp1-xp2);
B1=yp1-A1*xp1;
B2_erad=B1-bufferwidth*sqrt(A1^2+1);
B2_rog=B1-(bufferwidth+bufferwidth2)*sqrt(A1^2+1);

for i=1:N
    for j=1:M
        if j<A1*i+B1 && j>A1*i+B2_erad && olivegrowthprop(i,j)>0
            olivegrowthprop(i,j)=0;
        end
    end
end

kp=3; B=14.069; % Recursive Gompertz model parameters
beta=0.1; % Mean dispersal distance
CC=olivegrowthprop+(1-olivegrowthprop)*a; % Carrying capacity definition
CC(olivegrowthprop<0)=0;
tol=10^-8; % Numerical noise tolerence


% We use convolution theory and FFTs to simulate dispersal for the sake of
% efficiency.  It is necessary to create a padded domain for these
% calculations.  The dispersal kernel is defined over the padded domain.

[x2, y2] = meshgrid(-(M-1):(M-1),-(N-1):(N-1)); % Padded domain meshgrid
k = exp(-(x2.^2+y2.^2).^(1/2)/beta); % Define a 2D exponential kernel
kernel = zeros(3*N,3*M); % Padded landscape for the kernel
IN=(floor(N/2)+2):(floor(N/2)+2*N); % Locations of kernel in padded kernel
IM=(floor(M/2)+2):(floor(M/2)+2*M);
for i=1:length(IN)
    for j=1:length(IM)
        kernel( IN(i) , IM(j) )=k(i,j);
    end
end
fftkernel = fft2(fftshift(kernel)); % Take the fft of the kernel

% Define the parameters for the stratified dispersal. This dispersal works
% by weighted (by infection level) random patches being chosen for random
% dispersal events. The distance is chosen from random drawers from a 2D
% Gaussian kernel. The infection jump is seeded by the Gompertz initial
% infection. The random drawers are performed by the code pinky.m.
strattol=0.2; % Stratified dispersal tolerance (p)
D=20; % Dispersal standard deviation
Dist=exp(-(y2.^2+x2.^2)/(2*D^2)); % Gaussian distribution

output=zeros(N,M,timerun); % Initialise the output storage

u = zeros(N,M); % Define the population maxtrix
u(235,266) = exp(-B); % The initial outbreak  near Gallipoli

for i=1:timerun
    % Calculate the temporal dynamics
    u(olivegrowthprop<0)=0; % Set the population in the sea to zero
    u=CC.^(1-exp(-kp)).*u.^exp(-kp); % Implement the population model everywhere
    % Perform local deterministic dispersal
    U = zeros(3*N,3*M); % Padded landscape
    U((N+1):2*N,(M+1):2*M) = u; % Put the original square landscape in the middle
    U = ifft2(fft2(U).*fftkernel); % The convolution (i.e. the dispersal)
    % Perform non-local stochastic dispersal
    Rdisp=U.*rand(size(U)); % Assign the probability of random dispersers
    Rdisp(Rdisp<strattol)=0; % Find the patches that randomly disperse
    [indi,indj]=find(Rdisp); % Find thier indices
    vals=zeros(2,length(indi)); % Matrix for storing dispersal distances
    Rdispdest=zeros(size(U)); % Matrix for storing map of random dispersers
    for dispind=1:length(indi)
        for randind=1:randi(5)
            exitstat=0;
            while exitstat==0
                [vals(1,dispind),vals(2,dispind)]=pinky(-(M-1):(M-1),-(N-1):(N-1),Dist);
                if olivebuffer(indi(dispind)+vals(1,dispind),indj(dispind)+vals(2,dispind))>=0
                    Rdispdest(indi(dispind)+vals(1,dispind),indj(dispind)+vals(2,dispind))=exp(-B);
                    exitstat=1;
                end
            end
        end
    end
    U=U+Rdispdest; % Add in the new random dispersal patches.
    % Book keeping
    u = U((N+1):2*N,(M+1):2*M); % retain only non-padded region
    u(olivegrowthprop<=0)=0; % Set the population in the sea and non-olives to zero
    u(u<tol)=0; % Remove very small populations occuring from numerical noise
    % Remove infected olives
    for i2=1:N
        for j2=1:M
            if j2<A1*i2+B2_erad && j2>A1*i2+B2_rog && olivegrowthprop(i2,j2)>0 && u(i2,j2)>0
                if rand<eff
                    % If the infection is detected then remove all infected
                    % olive trees and don't replace them.
                    olivegrowthprop(i2,j2)=max([olivegrowthprop(i2,j2)-u(i2,j2),0]);
                    u(i2,j2)=0;
                    CC=olivegrowthprop+(1-olivegrowthprop)*a; % Update arrying capacity
                    CC(olivegrowthprop<0)=0;
                end
            end
        end
    end
    output(:,:,i)=u;
end


end



