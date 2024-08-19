%% Generate Koopman Traffic Modes
function [Psi]=GenerateKoopmanModes(data,mode1,mode2,save)
%% Load Data
clc; close all;
disp('Loading Data Set...')
tic
if strcmp(data,'Poi_Day_mean')
Data=dlmread('Daily PM2.5 in 13 cities.txt');
delay=3; dtype='Mean'; delt=1; delx=1;
hwy='day'; hwylength=76; xpath='x76.txt'; ypath='y76.txt';
end
toc
%% Compute KMD and Sort Modes
disp('Computing KMD via Hankel-DMD...')
tic
Avg=mean(Data,2);% Compute and Store Time Average
[eigval,Modes1,bo] = H_DMD(Data-repmat(Avg,1,size(Data,2)),delay); 
% [eigval,Modes1,bo] = H_DMD(Data-Avg,delay); % Compute HDMD
toc
disp('Sorting Modes...')
tic
% Sampling Frequency of PeMs/NGSIM Data is 5 Minutes/Seconds.
scatter(real(diag(eigval)),imag(diag(eigval))) 
omega=log(diag(eigval))./delt;% Compute Cont. Time Eigenvalues
Freal=imag(omega)./(2*pi);% Compute Frequency
[T,Im]=sort((1./Freal),'descend');% Sort Frequencies
omega=omega(Im); Modes1=Modes1(:,Im); bo=bo(Im); % Sort Modes
toc

%% Compute and Plot Modes 
disp('Computing and Plotting Modes...')
tic
[nbx,nbt]=size(Data); % Get Data Size
time=(0:nbt-1)*delt;% Specify Time Interval
Psi=zeros(nbx,nbt,mode2-mode1+1);
res=[]
for i=mode1:mode2 % Loop Through all Modes to Plot.
psi=zeros(1,nbt);% Preallocate Time Evolution of Mode.
omeganow=omega(i);% Get Current Eigenvalue.
bnow=bo(i);% Get Current Amplitude Coefficient.
parfor t=1:length(time) 
psi(:,t)=exp(omeganow*time(t))*bnow; % Evolve for Time Length.
end
psi=Modes1(1:nbx,i)*psi;% Compute Mode.
Psi(:,:,i)=psi;% Store & Output Modes
k=[xpath,ypath]
n=real(psi)
n1=mean(n,1)
n2=mean(n,2)
m=abs(psi)

% -------------------------------------------------------------------------%
FONTSIZE = 35;
TICKSIZE = 28;
if strcmp(hwy,'day')
[X,Y]=meshgrid(time./1,linspace(0,hwylength,nbx));
h=figure
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	
jFrame = get(h,'JavaFrame');	
pause(0.3);					
set(jFrame,'Maximized',1);	
pause(0.5);					
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		

s1=surfc(X,Y,real(psi));
set(s1,'LineStyle','none')

set(gca,'position',[0.1,0.15,0.60,0.78],'TickLabelInterpreter','latex','linewidth',2.5,'FontSize',30)
title(strcat('Mode #',num2str(i)),... 
                     'Interpreter','Latex','FontSize',30)
xlabel('Time (day)','Interpreter','tex','FontSize',30,'rotation',13); 
h=colorbar;
ylabel('Monitoring station','rotation',-20,'position',[-23 -10],'FontSize',30);%day

if strcmp(dtype,'Mean')
set(get(h,'title'),'string',{'¦Ìg/m^{3} per day'},'FontSize',30);
elseif strcmp(dtype,'Mean')
set(get(h,'title'),'string', {'¦Ìg/m^{3} per hour'});
end

end %End Modes to Plot Loop
%toc
disp('All Done')

end %End function
