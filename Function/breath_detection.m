function [Breathy] = breath_detection(Speech,Fs)
%breath_detection Finds the breath segments in the speech signal
% This function finds the breath sounds in the speech signals
%%%%% INPUT %%%%%%%%%%%%%%%%%
% Speech -  Speech signal for analysis
% Fs -      Sampling Frequency of the speech signal
%%%%%%%% OUTPUT %%%%%%%%%%%%%%
% Breathy is the signal showing the breath segments detected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Adding Voice box, Functions and voiceactivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('voicebox');
addpath('Function');
addpath('voiceactivity');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load svmmodel_ts_8p_25_6_19.mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(Speech,2)>1
    Speech=Speech(:,1);
end
% Covarep 
[Outs_Final,~,~,~,t] = VAD_Drugman(Speech,Fs,0);
t_s=[1:length(Speech)]/Fs;
vad_p=kron(Outs_Final(1),ones(ceil(Fs*0.0151),1));
vad_p=[vad_p' kron(Outs_Final(2:end),ones(1,ceil(Fs*0.01)))];
% Detecting the non-VAD part
Breathy_p=zeros(size(vad_p));
index=find(vad_p<0.08);
Breathy_p(index)=1;


diffe=diff(Breathy_p);
bindex=find(abs(diffe)==1);
if Breathy_p(1)==1
    bindex=[1  bindex];
end
if Breathy_p(end)==1
    bindex=[bindex  length(Breathy_p)];
end
plot(0:1/Fs:(length(Speech)-1)/Fs,Speech,'b',...
    0:1/Fs:(length(Breathy_p)-1)/Fs,0.75*Breathy_p,'r','LineWidth',2);
    axis([0 length(Speech)/Fs -1 1]);
    title(['Speech with VAD marked ']);
    xlabel('Time in seconds');
    set(gca,'FontSize',20);
    legend('Speech signal','VAD-detected');
    
vad_bframe=[];
Breathy=zeros(size(Breathy_p));
%% SVM for detected non-VAD
for i=1:length(bindex)/2
    b_frame=Speech(bindex(2*i-1):bindex(2*i));
    Overlap_time=10*10^(-3);
    Frame_length=Global_frame_length;
    Overlap_length=ceil(Overlap_time*Fs);
    No_frames=floor((length(b_frame)-Frame_length)/Overlap_length)+1;
    b_check=zeros(No_frames-1,1);
    for f_count=1:No_frames-1
%             range=bindex(2*i-1)+[(f_count-1)*Overlap_length+1:(f_count)*Overlap_length];
        frame=b_frame((f_count-1)*Overlap_length+1:(f_count-1)*...
        Overlap_length+Frame_length);
%********Feature extraction
        cepst=cepstrogram(frame,Fs);
        [m, n]=size(cepst);
        sf=[reshape(cepst,[1 m*n])];
        newf=(sf-feature_mean)./feature_std;
        for p_count=1:5
            [P_label(p_count),score]=predict(SVMModel,newf);
        end
        b_check(f_count)=mode(P_label);
    end
    if (sum(b_check)<2)
        Breathy_p(bindex(2*i-1):bindex(2*i))=0;
    else
        c_index=find(b_check==1);
         range=bindex(2*i-1)+[(min(c_index)-1)*Overlap_length+1 : ...
                (max(c_index))*Overlap_length];
            Breathy(range)=1;
        b_check(min(c_index):max(c_index))=1;
    end
    vad_bframe=[vad_bframe sum(b_check)/No_frames];
end
%*******Finding the breath edges    
diffe=diff(Breathy);
bindex=find(abs(diffe)==1);
if Breathy(1)==1
    bindex=[1  bindex];
end
if Breathy(end)==1
    bindex=[bindex length(Breathy)];
end

% Joining breathes which are less than 2 times minimum breath duration
for i=1:length(bindex)/2-1
    bb_length=bindex(2*i+1)-bindex(2*i);
    if bb_length < 2*Global_frame_length
        Breathy(bindex(2*i):bindex(2*i+1))=1;
    end
end 
diffe=diff(Breathy);
bindex=find(abs(diffe)==1);
if Breathy(1)==1
    bindex=[1  bindex];
end
if Breathy(end)==1
    bindex=[bindex length(Breathy)];
end

% % Removing breathes of smaller duration than minimum breath duration
for i=1:length(bindex)/2
    b_length=bindex(2*i)-bindex(2*i-1);
    if b_length < 100*10^(-3)*Fs
        Breathy(bindex(2*i-1):bindex(2*i))=0;
    end
end
diffe=diff(Breathy);
bindex=find(abs(diffe)==1);
if Breathy(1)==1
    bindex=[1  bindex];
end
if Breathy(end)==1
    bindex=[bindex length(Breathy)];
end
for i=1:length(bindex)/2
    bindex(2*i)=bindex(2*i)+(Frame_length-Overlap_length);
    Breathy(bindex(2*i-1):bindex(2*i))=1;
end 
diffe=diff(Breathy);
bindex=find(abs(diffe)==1);
if Breathy(1)==1
    bindex=[1  bindex];
end
if Breathy(end)==1
    bindex=[bindex length(Breathy)];
end
end

