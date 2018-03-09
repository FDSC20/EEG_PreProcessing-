clear
EEG=pop_biosig();%Open GDF
%EEG=pop_loadset(); unleashes lots of glitches
EEG.data=rmbase(EEG.data(:,:));%RemoveBaseline
EEG=pop_chanedit(EEG);
ALLEEG=[];
ALLCOM=[];
CURRENTSET=[];
xx=input('Quieres eliminar EMG? Y/N:','s');
if xx=='Y'
pop_eegplot( EEG, 1, 1, 1);%manual EMG component removal dont forget to hit have  dataset

pause();
EEG=EEGTMP;
clear EEGTMP;
end
EEG	=BandPass(EEG,.1,30); %filter 
clear xx
ALLEEG=EEG;

a=EEG;
a=BandPass(a,1,30); %filter 
CURRENTSET=2;
ALLEEG(double(CURRENTSET))=a;
EEGw=a;
aWICA=EEGw.data(:,:);
a=pop_runica(a,'verbose','off');% ica
%%
CURRENTSET=[3];
ALLEEG(double(CURRENTSET))=a;
pause();
EEGf=pop_editset(EEG);% change ica W from 1hz to .5 hz structure
pause();
EEGf=pop_selectcomps(EEGf); %select components
pause();
EEGf=pop_subcomp(EEGf); %eliminate components
pause();
pop_eegplot(EEGf,1,1,1); %display final signal
pause();
CURRENTSET=4;
%ALLEEG(double(CURRENTSET))=EEGf;
%%

[icasig, A, W] = fastica (aWICA);
z=zeros(size(icasig,1),2^ceil(log2(size(icasig,2)))-size(icasig,2));
icasig=[icasig z];

for i = 1:size(icasig,1)

	[SWA SWD]= swt(icasig(i,:),5,'coif5');

	k=(median(abs(SWA(:)))/.6745)*sqrt(2*log(length(SWA)));
	SWA(SWA>k)=0;

	icasig(i,:) = iswt(SWA(end,:),SWD,'coif5');
end
Data_wICA=single(A*icasig);
EEGw.data(:,:)=Data_wICA(:,1:length(aWICA));

pop_eegplot(EEGw,1,0)

EEGfepf=pop_epoch(EEGf);%Epoch split ICA
EEGfepf=pop_rmbase(EEGfepf,[-500 0]);%Remove baseline
EEGfepf.data=mean(EEGfepf.data(:,:,:),3);
EEGfepw=pop_epoch(EEGw);%Epoch split wICA
EEGfepw=pop_rmbase(EEGfepw,[-500 0]);%Remove baseline
EEGfepw.data=mean(EEGfepw.data(:,:,:),3);
for n =1:5
	figure
	subplot(211)
	plot(([1:length(EEGfepf.data)]/256)-.5,EEGfepf.data(n,:))
	hold on 
	plot([-.5 1],[0 0],'k')
	plot([0 0],[min(EEGfepf.data(n,:)),max(EEGfepf.data(n,:))],'k')
	hold off
	title('ICA','Fontsize',20)

	subplot(212)
	plot(([1:length(EEGfepw.data)]/256)-.5,EEGfepw.data(n,:))
	hold on 
	plot([-.5 1],[0 0],'k')
	plot([0 0],[min(EEGfepw.data(n,:)),max(EEGfepw.data(n,:))],'k')
	hold off
	title('wICA','Fontsize',20)
	pause()
end


%j





