clear
close all
clc

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% User Defined Range and Velocity of target
% define the target's initial position and velocity. 
% Note : Velocity remains contant

R = 100; % anything under max range
v = -20; % between +/- max velocity


%% FMCW Waveform Generation

% *%TODO* :
% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

c = 3e8; % speed of light
rangeRes = 1;
Rmax = 200;
B = c/2/rangeRes;
Tchirp = 5.5*2*Rmax/c; % For FMCW radar sweep time atleast 5-6 times max travel time
slope = B/Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
                                                         
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
  
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R + v*t(i);
    td(i) = 2*r_t(i)/c;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi* (fc*t(i) + slope*t(i)^2/2) );
    Rx(i) = cos(2*pi* (fc*(t(i)-td(i)) + slope*(t(i)-td(i))^2/2) );
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix_resh=reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
L = Nr;
Y = fft(Mix_resh);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);


 % *%TODO* :
% Take the absolute value of FFT output

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.


%plotting the range
figure ('Name','Range from First FFT')

 % *%TODO* :
 % plot FFT output 
subplot(2,1,1)
plot(P1)
title('1D Range FFT')
xlabel('Range Measurement')
ylabel('|P1(f)|')
axis ([0 200 0 1]);

subplot(2,1,2)
plot(t(1:Nr),Mix_resh(:,1));
ylabel('f')
xlabel('t(s)')

saveas(gcf,'range_FFT.png')


%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);
figure,surf(Mix,'linestyle', 'none');
colorbar;
xlabel('chirp #')
ylabel('sample #')
title('IF Signal');

saveas(gcf,'IF_signal.png')


% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions

doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM,'linestyle', 'none');
colorbar;
xlabel('velocity estimate (m/s)')
ylabel('range estimate (m)')
title('Range-Doppler FFT Output');

saveas(gcf,'RD_FFT.png')

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.

Tr = 10;
Td = 8;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% offset the threshold by SNR value in dB
offset = 6;

%Create a vector to store noise_level for each iteration on training cells
noise_level = 0;


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

%Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR

RDM_pow = 10.^(RDM/10);

for i = Tr + Gr + 1:(Nr/2) - (Tr + Gr)
    for j = Td + Gd + 1:Nd - (Td + Gd)
        noise_level = 0;
        for p = i - (Tr + Gr):i + (Tr + Gr)
            for q = j - (Td + Gd):j + (Td + Gd)
                 if((abs(i - p) > Gr) || (abs(j - q) > Gd))
                        noise_level = noise_level + RDM_pow(p, q);   
                 end
            end
        end

        avg = noise_level / ((2 * (Td + Gd + 1) * 2 * (Tr + Gr + 1) - (Gr * Gd) - 1));
        threshold = pow2db(avg) + offset;
        CUT = RDM(i, j);
        if(CUT > threshold)
            RDM(i, j) = 1;
        else
            RDM(i, j) = 0;
        end
    end
end

% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 

for i = 1: Tr + Gr % rows
    RDM(i, :) = 0;
    RDM(Nr / 2 - i - 1:Nr / 2, :) = 0;
end
for i = 1:Td + Gd % cols
    RDM(:, i) = 0;
    RDM(:, Nd - i - 1:Nd) = 0;
end

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name','CA-CFAR Filtered RDM');
surf(doppler_axis, range_axis, RDM, 'linestyle', 'none');

colorbar;
view(0,-90)
xlabel('velocity estimate (m/s)')
ylabel('range estimate (m)')
title('Range-Doppler CA-CFAR Output');

saveas(gcf,'CACFAR.png')

function ydB = pow2db(y)
%POW2DB   Power to dB conversion
%   YDB = POW2DB(Y) convert the data Y into its corresponding dB value YDB
%
%   % Example:
%   %   Calculate ratio of 2000W to 2W in decibels
%
%   y1 = pow2db(2000/2)     % Answer in db

%#codegen
cond = all(y(:)>=0);
if ~cond
    coder.internal.assert(cond,'signal:pow2db:InvalidInput');
end

% We want to guarantee that the result is an integer
% if y is a negative power of 10.  To do so, we force
% some rounding of precision by adding 300-300.

ydB = (10.*log10(y)+300)-300;
end

