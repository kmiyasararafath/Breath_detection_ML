function [cepst]=cepstrogram(frame,Fs)
%CEPSTROGRAM to find the cepstrogram of a matrix
% INPUT
% Frame   - Frame whose cepstrogram needs to be found
% Fs      - Sampling frequency
% OUTPUT
% cepst   - The Ceptrogram matrix

emph_frame=filter([1 -0.95],1,frame);
Frame_time=10*10^(-3); 
Overlap_time=5*10^(-3);
Frame_length=ceil(Frame_time*Fs);
Overlap_length=ceil(Overlap_time*Fs);
No_of_sub_frames=floor((length(emph_frame)-Frame_length)/Overlap_length)+1;
fft_length=2^ceil(log(Frame_length)/log(2));
no_mel_fil_banks=26;
cepst=zeros(15,No_of_sub_frames);
for sf_count=1:No_of_sub_frames
    frame=emph_frame((sf_count-1)*Overlap_length+1:(sf_count-1)*...
    Overlap_length+Frame_length);
    cepst(:,sf_count)=find_mfcc(frame,Fs,no_mel_fil_banks,fft_length);
    mesb=mean(cepst(:,sf_count));
    cepst(:,sf_count)=cepst(:,sf_count)-mesb;
end
end

function [ cc ] = find_mfcc(speech_frame,fs,p,n)
%FIND_MFCC Find the Mel frequency cepstrum coefficients of a speech frame
% 1. Take the Fourier transform of (a windowed excerpt of) a signal.
% 2. Map the powers of the spectrum obtained above onto the mel scale, 
%       using triangular overlapping windows.
% 3. Take the logs of the powers at each of the mel frequencies.
% 4. Take the discrete cosine transform of the list of mel log powers, 
%       as if it were a signal.
% 
% The MFCCs are the amplitudes of the resulting spectrum.
% 
% Calcuate the Mel-frequency Cepstral Coefficients of Voice box help of
% melbankm

f=rfft(speech_frame,n);        % rfft() returns only 1+floor(n/2) coefficients
x=melbankm(p,n,fs);            % n is the fft length, p is the number of filters wanted
z=log(x*abs(f).^2);            % multiply x by the power spectrum
cc1=dct(z);                    % take the DCT
cc=cc1(2:16);                  % Selecting 13 of the 26 parameters

end








