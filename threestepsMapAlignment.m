'Genotype can be: Wildtype, ephrinA3KIplus, ephrinA3KIKI, EphA3KIKI, EphA3KIplus'
n = input('Genotype: ');
X=linspace(0,100,100);

switch n
    case 'Wildtype'
        ephrinA3KI=0;
        Genotype='WT'
        EphA3=0;
        EphA4=1.05;
        EphA5=0.14*exp(0.018*X);
        EphA6=0.09*exp(0.029*X);
        
    case 'ephrinA3KIplus'    
        ephrinA3KI=0.22;
        Genotype='ephrinA3KIplus'
        EphA3=0;
        EphA4=1.05;
        EphA5=0.14*exp(0.018*X);
        EphA6=0.09*exp(0.029*X);
         
    case 'ephrinA3KIKI'
        ephrinA3KI=0.44
        Genotype='ephrinA3KIKI'
        EphA3=0
        EphA4=1.05
        EphA5=0.14*exp(0.018*X)
        EphA6=0.09*exp(0.029*X)     
        
    case 'EphA3KIKI'
        ephrinA3KI=0;
        EphA3=1.86;
        Genotype='EphA^{KI/KI}'
        threshold=8.68;
        EphA4=1.05;
        EphA5=0.14*exp(0.018*X);
        EphA6=0.09*exp(0.029*X);
          
    case 'EphA3KIplus'
        ephrinA3KI=0;
        EphA3=0.93;
        Genotype='Epha3^{KI/+}'
        threshold=8.68;
        EphA4=1.05;
        EphA5=0.14*exp(0.018*X);
        EphA6=0.09*exp(0.029*X);
        
    case 'EphA3KIxephrinA3'
        ephrinA3KI=0.22;
        EphA3=0.93;
        Genotype='EphA3KIxephrinA3'
        threshold=8.68;
        EphA4=1.05;
        EphA5=0.14*exp(0.018*X);
        EphA6=0.09*exp(0.029*X);
        
    otherwise
        disp('Invalid genotype')
end

%%
% 1D simulation of the development of retinocollicular maps in wild-type species
% From Koulakov and Tsigankov 2004

%Intial parameters

NIterations =1000000  %number of measurement of the axon distribution
N = 100; %number of axons 
alpha = 200; %simulation parameter, the transition probability coefficient in rostral-caudal direction
gamma = 1; % the overall strength of the activity-dependent component            
%Constants for the activity-dependent portion
d=3;%Interaction distance in the superior colliculus
R=0.11*N;%Interaction distance in the retina

%% Creating gradients

%Retinal EphA
X=linspace(0,100,100);
RA= (EphA5 + EphA6 + EphA4); %Equation from Reber et al., 2004
RA(RA<0)=0;
ISL2=randi(2,100,1)-1;
ISL2=ISL2';
A3=ISL2.*EphA3;
[Isl2posIdx]=find(A3);
RA=plus(RA,A3);

%Collicular ephrinA
X=linspace(0,100,100);
LA = exp((X-N)/N)-exp((-X-N)/N); %Equation from Tsigankov and Koulakov, 2010
LA(LA<0)=0;

%Generating the cortical EphA gradient
X = 1:N;
RASC =(exp(-X/N)-exp((X-2*N)/N))+1; %Equation from Tsigankov and Koulakov, 2010 plus constant
RASC=fliplr(RASC);
RASC(RASC<0)=0;

%Generating the retinal ephrinA gradient: Exponential fit of ISH data (unpublished data, Reber)
X=linspace(0,10,100);
gA2=0.5462*exp(0.08403*X); %ephrinA2
gA3=0.006833*X+0.44;%ephrinA3
gA5=0.5612*exp(0.1447*X);%ephrinA5
gA5=gA5/1.2; % 20% used for axon-axon interaction 
gEphrins=gA2+gA5+gA3; %summed ephrinAs
N3=ISL2.*ephrinA3KI;
gEphrins=plus(gEphrins,N3);

%% Simulation of development process for the retinocollicular map
%setting initial configuration
 
RofL = randperm(N); 
LofR = zeros(N,1);
for i=1:N     
    LofR(i) = find(RofL == i);
end
 
%simulation of development process for the retinocollicular map

for t = 1:NIterations%(iterator)
    %Specify the coordinates of the axons in the SC.
    iL = ceil(rand*(N)); % ceil rounds up towards an integer.  
    iL2 = ceil(rand*(N)); % Picks second axon NT index.
   
    % Calculate Echem
    EA=alpha* ( RA(RofL(iL2)) - RA(RofL(iL)) ) * (LA(iL2)-LA(iL)); 
    Echem=EA;
    
    % Source position
    RetiL=RofL(iL);
	RetiL2=RofL(iL2);
    
    %Calculate C with RetiL and RetiL2
    Rret=abs(RetiL-RetiL2);%
    C=exp(-Rret/R); 
    
    %Calculate U with iL and iL2 
    Rsc_sq=(iL-iL2)^2;
    U=exp(-Rsc_sq/(2*d^2));
    EactA= U.*C; 
    
    %This is the change in Eact. 
    Eact=-(gamma/2)*(sum(sum(EactA)));

    % Calculate Etot
    Etot= Eact+Echem;
    
    p=1/(1+exp(4*Etot));
    if (rand < p) %Switches the two axons
        tmp = LofR(RofL(iL)); %Extracting the position of the first axon
        LofR(RofL(iL)) = LofR(RofL(iL2)); %Replacing it with the position of the second axon
        LofR(RofL(iL2)) = tmp; %Putting the first axon into the former
        %position of the second axon
        tmp = RofL(iL); % Extracts the value of the first axon
        RofL(iL) = RofL(iL2); %Replacing it with the value of the
        % second axon
        RofL(iL2) = tmp; %Putting the first axon's value into the second
        %axon
    end   
end


 
 MapLR = zeros(N,N);
        for j=1:N
        MapLR(j,N+1-RofL(N+1-j)) = 1;
        end
        
%Extracting indexes
inds = find(MapLR==1);
[row, col] = ind2sub(size(MapLR),inds); 

%% Extracting indexes
inds = find(MapLR==1);
[row, col] = ind2sub(size(MapLR),inds); 
    
%% Transposing retinal ephrinAs to the SC

[RCaxis,col2] = find(MapLR>0); %Extracting indexes
LASC(RCaxis)=gEphrins; %Attributing ephrinA values to the positions in the retinocollicular map

%% Generation of the corticocollicular map
 
RofL2 = randperm(N); 
LofR2 = zeros(N,1);
for i=1:N     
    LofR2(i) = find(RofL2 == i);
end

 
%simulation of development process for the retinocollicular map

for t = 1:NIterations%(iterator)
    %Specify the coordinates of the axons in the SC.
    iL = ceil(rand*(N)); % ceil rounds up towards an integer.  
    iL2 = ceil(rand*(N)); % Picks second axon NT index.
   
    % Calculate Echem
    EA=alpha* ( RASC(RofL2(iL2)) - RASC(RofL2(iL)) ) * (LASC(iL2)-LASC(iL)); 
    Echem=EA;
    
    % Source position
    RetiL=RofL2(iL);
	RetiL2=RofL2(iL2);
    
    %Calculate C with RetiL and RetiL2
    Rret=abs(RetiL-RetiL2);%
    C=exp(-Rret/R); 
    
    %Calculate U with iL and iL2
    Rsc_sq=(iL-iL2)^2;
    U=exp(-Rsc_sq/(2*(d^2)));
    EactA= U.*C; 
    
    %This is the change in Eact. 
    Eact=-(gamma/2)*(sum(sum(EactA)));

    % Calculate Etot
    Etot= Eact+Echem;
    
    p=1/(1+exp(4*Etot));
    if (rand < p) %Switches the two axons
        tmp = LofR2(RofL2(iL)); %Extracting the position of the first axon
        LofR2(RofL2(iL)) = LofR2(RofL2(iL2)); %Replacing it with the position of the second axon
        LofR2(RofL2(iL2)) = tmp; %Putting the first axon into the former
        %position of the second axon
        tmp = RofL2(iL); % Extracts the value of the first axon
        RofL2(iL) = RofL2(iL2); %Replacing it with the value of the
        % second axon
        RofL2(iL2) = tmp; %Putting the first axon's value into the second
        %axon
    end   
end


 
 MapLRSC = zeros(N,N);
        for j=1:N
        MapLRSC(j,N+1-RofL2(N+1-j)) = 1;
        end
        
%Extracting indexes
inds = find(MapLRSC==1);
[row2, col2] = ind2sub(size(MapLRSC),inds); 




%% Plotting all figures
figure
%Plotting retinal EphA gradient

subplot(2,4,1)
a=area(RA);
a.FaceColor=[.5 .5 .5];
a.EdgeColor=[.5 .5 .5];
set(gca,'XTick',[0,20,40,60,80,100]);
set(gca,'XtickLabel',{'N';'20'; '40'; '60';'80';'T'});
ylim([0 5])
xlabel('Retinal axis (%)');
ylabel('Measured Relative [Ephas]');
set(gca,'FontSize',12)
box off
title([Genotype ' retinal Ephas gradients']);

%Plotting collicular ephrinA gradient
subplot(2,4,2)
b=area(LA);
b.FaceColor=[.8 .8 .8];
b.EdgeColor=[.8 .8 .8];
set(gca,'XTick',[0,20,40,60,80,100])
set(gca,'XtickLabel',{'R';'20'; '40'; '60';'80';'C'})
axis([0 100 0 1])
ylim([0 5])
xlabel('Collicular axis (%)')
ylabel('Estimated Relative [Efnas]')
set(gca,'FontSize',12)
box off
title([Genotype ' collicular Efnas gradients'])

% Plotting the retinocollicular map 
subplot(2,4,3)
scatter(col,row,'filled');
s=scatter(col,row,'filled');
s.MarkerFaceColor = [0.5 0.5 0.5];
set(gca,'XDir','Reverse');
set(gca,'XTick',[0,20,40,60,80,100]);
set(gca,'XtickLabel',{'T';'80'; '60'; '40';'20';'N'})
set(gca,'YTick',[0,20,40,60,80,100]);
set(gca,'YtickLabel',{'R';'20'; '40'; '60';'80';'C'})
xlabel('Retinal axis (%)')
ylabel('Collicular axis (%)')
set(gca,'FontSize',12)
title([Genotype ' retinocollicular map']);

%Plotting the collicular ephrinA gradient
subplot(2,4,4)
b=area(LASC);
b.FaceColor=[.8 .8 .8];
b.EdgeColor=[.8 .8 .8];
set(gca,'XTick',[0,20,40,60,80,100])
set(gca,'XtickLabel',{'R';'20'; '40'; '60';'80';'C'})
axis([0 100 0 5])
xlabel('Collicular position (%)')
set(gca,'FontSize',12)
ylabel('Relative [Efnas]')
box off
title(['Transposed retinal Efnas in ' Genotype])

%Plotting the cortical EphA gradient
subplot(2,4,5)
a=area(RASC);
a.FaceColor=[.5 .5 .5];
a.EdgeColor=[.5 .5 .5];
set(gca,'XTick',[0,20,40,60,80,100]);
set(gca,'XtickLabel',{'M';'20'; '40'; '60';'80';'L'});
ylim([0 5]);
set(gca,'FontSize',12)
xlabel('Cortical axis (%)');
box off
ylabel('Estimated Relative [Ephas]');
title([Genotype ' cortical Ephas gradients']);

% Plotting the corticocollicular map 
subplot(2,4,6)
s=scatter(col2,row2,'filled');
s.MarkerFaceColor = [0.5 0.5 0.5];
set(gca,'XTick',[0,20,40,60,80,100]);
set(gca,'YtickLabel',{'C';'80'; '60'; '40';'20';'R'})
set(gca,'YTick',[0,20,40,60,80,100]);
set(gca,'YDir','Reverse');
set(gca,'XtickLabel',{'L';'80'; '60'; '40';'20';'M'})
xlabel('Cortical axis (%)');
ylabel('Collicular axis (%)');
title([Genotype ' corticocollicular map']);
set(gca,'FontSize',12)

%Plotting both retino and corticocollicular map
subplot(2,4,7)
rowbis=flipud(row);
s=scatter(col,rowbis,'filled');
s.MarkerFaceColor = [0 0 0];
ylabel('Collicular axis (%)')
title([Genotype ' retino and cortico-collicular map']);
xlabel('Retinal axis | Cortical axis (%)')
hold on

row2bis=(((ones(100,1))-row2))+100;
s=scatter(col,row2bis,'filled');
s.MarkerEdgeColor=[0 0 0];
s.MarkerFaceColor = [1 1 1];
set(gca,'XTick',[0,20,40,60,80,100]);
set(gca,'YTick',[0,20,40,60,80,100]);
set(gca,'FontSize',12)
