function [O] = Dance(pData,option)
%% INSTRUCTION
% 1) set path of interest:
%    1) Drag path into Command Window
%    2) Type in [p=cd] and press ENTER
% 2) Drag in MWT_Programs folder (located under your MultiWormTrackerPortal
% folder)
% 3) Type in Dance(p,cd)
%
% Variables:
%   p = path of interest, different for each analysis
%   pF = path to function folder, you can cd to your
% Coding History: 
%   created from AnalysisOne(pExp,pFun,ID)
%   created by Conny Lin 2013

%% STEP1: CLEAN UP
clearvars -except option pData;
%% STEP1: PATHS [Need flexibility] 
a = regexp(userpath,':','split'); cd(a{:,1}); PathCommonList; clearvars a;
addpath(genpath(paths.pFunMWT));
pMaestro = '/Volumes/AWLIGHT/MWT/Archived Analysis/Maestro';
pRosePick = '/Volumes/Rose/MWT_Analysis_20130811';
pRose = '/Volumes/Rose/MWT_Analysis_20131020_Liquid';
pSet = '/Users/connylinlin/Documents/MATLAB/MATLAB MWT/MatSetfiles';
pSum = '/Volumes/Rose/MultiWormTrackerPortal/Summary';
pSave = paths.pSave;


%% STEP2: INTERPRET INPUT AND DEAL OUTPUT
HomePath = pData;
O = []; % [CODE] Define output O % temperary assign O = nothing
% STEP2: OPTION INTERPRETATION
% STEP2A: CREATE OPTION LIST 
ID = {'ShaneSpark';'LadyGaGa';'DrunkPosture'};
% STEP2B: CREATE OPTION NOT NEEDING TIME INPUTS
OptionNoTimeInput = {'ShaneSpark'};
% STEP2C: OPTION SELECTION 
i = celltakeout(regexp(ID,option),'logical');
if sum(i)~=1;
    display 'option entered not found'l
    display 'please select from the following list:';
    [show] = makedisplay(ID,'bracket'); disp(show);
% ask for analysis ID
display 'Enter analysis number,press [ENTER] to abort';
a = input(': ');
if isempty(a) ==1; return; else option = ID{a}; end
end
% STEP2D: TIME INTERVAL
optioni = celltakeout(regexp(OptionNoTimeInput,option),'logical');
if optioni ==0; timeintervalstatus = 'input'; 
else timeintervalstatus = 'noneed'; end
% STEP2E: DEFINE CHOROUTPUTS AND CHOR CODE
chorcode = option;
switch option
    case 'ShaneSpark'; fnvalidate = {'*.trv';'*shanespark.dat'}; 
    case 'DrunkPosture'; fnvalidate = {'*drunkposture.dat'};
    case 'LadyGaGa'; fnvalidate = {'*evan.dat';'*.sprevs'};
end

%% STEP3: DEFINE EXP OBJECTIVES
%% STEP3A: GET BASIC EXP INFO
display 'Checking if MWTExpDataBase is already loaded';
if exist('pExpfD','var')==0||exist('ExpfnD','var')==0||...
        exist('pMWTfD','var')==0|| exist('MWTfnD','var')==0;
    display 'Loading MWTExpDataBase...';
    [A] = MWTDataBaseMaster(HomePath,'GetStdMWTDataBase');
    pExpfD = A.pExpfD; ExpfnD = A.ExpfnD; GfnD = A.GfnD;
    pGfD = A.pGfD; pMWTf = A.pMWTf; MWTfn = A.MWTfn;
end
%% STEP4A: SELECT TARGET EXPERIMENTS
% search for target paths
display 'Search for...'
display 'MWT files [M], Exp files [E], Group name [G] or no search [A]?';
searchclass = input(': ','s');
switch searchclass
    case 'M' % search for MWT
        display 'option under construction'
    case 'G' % search for group name
        display 'option under construction'
        
    case 'E' % search for Experiment folders
        display 'Enter search term:';
        searchterm = input(': ','s');
        % find matches
        k = regexp(ExpfnD,searchterm,'once');
        searchindex = logical(celltakeout(k,'singlenumber'));
        pExpfS = pExpfD(searchindex);
        ExpfnS = ExpfnD(searchindex);
        disp(ExpfnS);
        moresearch = input('Narrow down search (y=1,n=0)?: ');
        while moresearch==1;
            display 'Enter search term:';
            searchterm = input(': ','s');
            % find matches
            k = regexp(ExpfnS,searchterm,'once');
            searchindex = logical(celltakeout(k,'singlenumber'));
            pExpfS = pExpfS(searchindex);
            ExpfnS = ExpfnS(searchindex);
            disp(ExpfnS);
            moresearch = input('Narrow down search (y=1,n=0)?: ');
        end
        pExpfT = pExpfS;
        ExpfnT = ExpfnS;
        display 'Target experiments:';
        disp(ExpfnT); 
        O.ExpfnT = ExpfnT; % export
    case 'A'
        pExpfT = pExpfD;
        ExpfnT = ExpfnD;
end  
if isempty(ExpfnT)==1; display 'No target experiments'; return; end
%% STEP4B. GET EXP AND GROUP INFO FROM TARGET EXP
[A] = MWTDataBaseMaster(pExpfT,'GetExpTargetInfo');
Gfn = A.GfnT; pGf = A.pGfT; pMWTfT = A.pMWTfT; MWTfnT = A.MWTfnT;
%% STEP2C: CHOOSE GROUP TO ANALYZE
gnameU = unique(Gfn);
a = num2str((1:numel(gnameU))');
[b] = char(cellfunexpr(gnameU,')'));
c = char(gnameU);
show = [a,b,c];
disp(show);
display 'Choose group(s) to analyze separated by [SPACE]';
display 'enter [ALL] to analyze all groups';
i = input(': ','s');
if strcmp(i,'ALL'); gnamechose = gnameU;
else k = cellfun(@str2num,(regexp(i,'\s','split')'));
    gnamechose = gnameU(k); 
end
% STEP2D: SORT GROUP DISPLAY SEQUENCE
[show] = makedisplay(gnamechose,'bracket');
disp(show);
display 'is this the sequence to be appeared on graphs (y=1 n=0)';
q2 = input(': ');
while q2 ==0;
    display 'Enter the correct sequence separated by [SPACE]';
    s = str2num(input(': ','s'));
    gnamechose = gnamechose(s,1);
    [show] = makedisplay(gnamechose,'bracket');
    disp(show);
    q2 = input('is this correct(y=1 n=0): ');
end
% STEP2B: GET MWT FILES UNDER SELECTED GROUP
% get experiment count for each group
MWTfG = [];
for g = 1:numel(gnamechose);
    str = 'Searching MWT files for group [%s]';
    display(sprintf(str,gnamechose{g}));
    % find MWT files under group folder
    search = ['\<',gnamechose{g},'\>'];
    i = logical(celltakeout(regexp(Gfn,search),'singlenumber'));
    pGfT = pGf(i);
    [MWTfn,pMWTf] = cellfun(@dircontentmwt,pGfT,'UniformOutput',0);
    MWTfn = celltakeout(MWTfn,'multirow');
    pMWTf = celltakeout(pMWTf,'multirow');
    a = celltakeout(regexp(MWTfn,'^(\d{8})','match'),'match');
    expcount = numel(unique(a));
    str = 'Got %d MWT files from %d experiments';
    display(sprintf(str,numel(MWTfn),expcount));
    % create output file
    MWTfG.(gnamechose{g})(:,1:2) = [MWTfn,pMWTf];
end

%% STEP4B: USER INPUT-TIME INTERVALS
% works only for the same run conditions

switch timeintervalstatus
    case 'input'
        display 'Enter analysis time periods: ';
        tinput = input('Enter start time, press [Enter] to use MinTracked: ');
        intinput = input('Enter interval, press [Enter] to use default (10s): ');
        tfinput = input('Enter end time, press [Enter] to use MaxTracked: ');
        display 'Enter duration after each time point to analyze';
% (optional) survey duration after specifoied target time point
        durinput = input('press [Enter] to analyze all time bewteen intervals: '); 
    case 'auto'
        tinput = []; int = []; tf = []; dur = [];
    case 'noneed';
end
%% STEP4C: CREATE OUTPUT FOLDER
display 'Name your analysis output folder';
a = clock;
time = [num2str(a(4)),num2str(a(5))];
name = [input(': ','s'),'_',datestr(now,'yyyymmdd'),time];
mkdir(pSave,name);
pSaveA = [pSave,'/',name];


%% STEP5: CHOR 
% STEP4B: CHECK IF CHOR HAD BEEN DONE
display ' '; display 'Checking chor outputs...'
% prepare pMWTf input for validation
A = celltakeout(struct2cell(MWTfG),'multirow');
pMWTfT = A(:,2); MWTfnT = A(:,1); pMWTf = pMWTfT; MWTfn = MWTfnT;
% check chor ouptputs
pMWTfC = {}; MWTfnC = {}; 
for v = 1:size(fnvalidate,1);
    [valname] = cellfunexpr(pMWTf,fnvalidate{v});
    [fn,path] = cellfun(@dircontentext,pMWTf,valname,'UniformOutput',0);
    novalfn = cellfun(@isempty,fn);
    if sum(novalfn)~=0;
        pMWTfnoval = pMWTf(novalfn); MWTfnoval = MWTfn(novalfn);
        pMWTfC = [pMWTfC;pMWTfnoval]; MWTfnC = [MWTfnC;MWTfnoval];
    end
end
pMWTfC = unique(pMWTfC); MWTfnC = unique(MWTfnC);
% reporting
if isempty(pMWTfC)==0;
    str = 'Need to Chore %d MWT files';
    display(sprintf(str,numel(pMWTfC)));
else display 'All files have required Chor outputs';
end
% STEP4C: RUN CHORE
if isempty(pMWTfC)==0; [~] = chormaster(pMWTfC,chorcode);
else display 'No Chor needs to be ran.'; end

%% STEP6: ANALYSIS SWITCH BOARD
display ' '; display 'Importing data generated by Choreography...';
switch option
    case 'ShaneSpark'
        %% STEP4: STATS AND GRAPHING
        %% STEP3E: EXCLUDE CHOR PROBLEM MWT FILES FROM ANALYSIS
        % RE-CHECK CHOR OUTPUTS
        display ' '; display 'Double checking chor outputs...'
        % prepare pMWTf input for validation
        A = celltakeout(struct2cell(MWTfG),'multirow');
        pMWTfT = A(:,2); MWTfnT = A(:,1);
        pMWTf = pMWTfT; MWTfn = MWTfnT;
        % check chor ouptputs
        pMWTfC = {}; MWTfnC = {}; 
        for v = 1:size(fnvalidate,1);
            [valname] = cellfunexpr(pMWTf,fnvalidate{v});
            [fn,path] = cellfun(@dircontentext,pMWTf,valname,'UniformOutput',0);
            novalfn = cellfun(@isempty,fn);
            if sum(novalfn)~=0;
                pMWTfnoval = pMWTf(novalfn); MWTfnoval = MWTfn(novalfn);
                pMWTfC = [pMWTfC;pMWTfnoval]; MWTfnC = [MWTfnC;MWTfnoval];
            end
        end
        pMWTfC = unique(pMWTfC); MWTfnC = unique(MWTfnC);
        % reporting
        if isempty(pMWTfC)==1; display 'All files have required Chor outputs';
        elseif isempty(pMWTfC)==0;
            str = 'Chore unsuccessful in %d MWT files';
            display(sprintf(str,numel(pMWTfC)));disp(MWTfnC);
            % STEP3F: EXCLUDE PROBLEM MWT FILES
            display 'Excluding problem MWT from analysis';
            gname = fieldnames(MWTfG);
            for x = 1:numel(MWTfnC) % each problem folders
                for g = 1:numel(gname)
                    A = MWTfG.(gname{g})(:,1);
                    i = logical(celltakeout(regexp(A,MWTfnC{x}),'singlenumber'));
                    if sum(i)>0; 
                        str = '>removing [%s]';
                        display(sprintf(str,MWTfnC{x}));
                        MWTfG.(gname{g})(i,:)=[]; % remove that from analysis
                    end
                end
            end
        end
        %% STEP4A: IMPORT .TRV 
        % get MWT list
        A = celltakeout(struct2cell(MWTfG),'multirow');
        pMWTfT = A(:,2); MWTfnT = A(:,1);
        pMWTf = pMWTfT; MWTfn = MWTfnT;
        % creating MWTftrv legend
        % ALegend = {1,'MWTfile name';3,'RawData'};
        % trvL = {1,'time';2,'N?';3,'Noresponse';4,'NReversed';5,'RevDist'};
        % import 
        A = MWTfn;
        for m = 1:size(pMWTf,1);
            [fn,~] = dircontentext(pMWTf{m},'*.trv'); 
            a = dlmread(fn{1},' ',0,0);
            i = [1,3:5,8:10,12:16,19:21,23:27]; % index to none zeros
            A{m,2} = a(:,i); % remove zeros
        end
        MWTfnImport = A;
        %% STEP4X: CHECK TAP CONSISTENCY
        [r,c] = cellfun(@size,MWTfnImport(:,2),'UniformOutput',0);
        rn = celltakeout(r,'singlenumber');
        rfreq = tabulate(rn);
        rcommon = rfreq(rfreq(:,2) == max(rfreq(:,2)),1);
        str = 'Most common tap number = %d';
        display(sprintf(str,rcommon));
        rproblem = rn ~= rcommon;
        if sum(rproblem)~=0;
            MWTfnP = MWTfn(rproblem); 
            str = 'The following MWT did not have the same tap(=%d)';
            display(sprintf(str,rcommon)); disp(MWTfnP);
            display 'Removing from analysis...'; gname = fieldnames(MWTfG);
            for x = 1:numel(MWTfnP) % each problem folders
                for g = 1:numel(gname)
                    A = MWTfG.(gname{g})(:,1);
                    i = logical(celltakeout(regexp(A,MWTfnP{x}),'singlenumber'));
                    if sum(i)>0; 
                        str = '>removing [%s] from group file';
                        display(sprintf(str,MWTfnP{x}));
                        MWTfG.(gname{g})(i,:)=[]; % remove that from analysis
                    end
                    k = logical(celltakeout(regexp(MWTfnImport(:,1),MWTfnP{x}),'singlenumber'));
                    if sum(k)>0; 
                        str = '>removing [%s] from imported file';
                        display(sprintf(str,MWTfnP{x}));
                        MWTfnImport(k,:)=[]; % remove that from analysis
                    end
                end
            end
        end
        %% STEP4B: MAKING SENSE OF TRV 
        % reload MWT files
        A = celltakeout(struct2cell(MWTfG),'multirow');
        pMWTf = A(:,2); MWTfn = A(:,1);
        % calculate
        B = [];
        B.MWTfn = MWTfn;
        A = MWTfnImport;
        for m = 1:size(pMWTf,1);
            B.X.TapTime(:,m) = A{m,2}(:,1); 
            B.N.NoResponse(:,m) = A{m,2}(:,3);
            B.N.Reversed(:,m) = A{m,2}(:,4); 
            B.N.TotalN(:,m) = B.N.Reversed(:,m)+B.N.NoResponse(:,m);
            B.Y.RevFreq(:,m) = B.N.Reversed(:,m)./B.N.TotalN(:,m);
            B.Y.RevDist(:,m) = A{m,2}(:,5); 
        %     B.Y.SumRevDist(:,m) = B.Y.RevDist(:,m).*B.N.Reversed(:,m); 
            B.Y.RevDistStd(:,m) = B.Y.RevDist(:,m)/B.Y.RevDist(1,m);
            B.Y.RevFreqStd(:,m) = B.Y.RevFreq(:,m)/B.Y.RevFreq(1,m); % freqStd
        end
        Stats = B;
        O.Stats = Stats;
        %% STEP4B: PREPARE GRAPH SUMMARY FOR SUBPLOTS
        Graph = [];
        MWTfnimport = Stats.MWTfn;
        M = fieldnames(Stats.Y);
        gnameL = gnamechose;
        for m = 1:numel(M);% for each measure
            Y = []; X = []; E = [];
            for g = 1:numel(gnameL); % for each group
                gname = gnameL{g};
                %pMWTf = MWTfG.(gname)(:,2); 
                MWTfn = MWTfG.(gname)(:,1);
                N = size(MWTfG.(gname),1);
                [~,i,~] = intersect(MWTfnimport,MWTfn);
                %Graph.Y.(M{m})(:,g) = mean(Stats.Y.(M{m})(:,i),2);
                %Graph.E.(M{m})(:,g) = std(Stats.Y.(M{m})(:,i)')'./sqrt(N);
                Graph.Raw.(M{m}).(gname) = Stats.Y.(M{m})(:,i);
                Graph.Y.(M{m})(:,g) = nanmean(Stats.Y.(M{m})(:,i),2);
                Graph.E.(M{m})(:,g) = nanstd(Stats.Y.(M{m})(:,i)')'./sqrt(N);
                Graph.X.(M{m})(:,g) = (1:size(Stats.X.TapTime,1));
            end
        end 
        O.Graph = Graph;
        %% STEP4C.SUBPLOTS
        % source code: LadyGaGaSubPlot(MWTftrvG,pExp,SavePrefix)
        % define universal settings
        % switch graphing sequence
        i = [2,3,1,4];
        k = fieldnames(Stats.Y)';
        M = k(i);
        groupname = gnamechose';
        groupsize = numel(gnamechose);
        gnshow = regexprep(groupname,'_',' ');
        xmax = size(Graph.X.(M{m}),1)+1;
        subplotposition(1,1:4) = [0.065 0.55 0.4 0.4];
        subplotposition(2,1:4) = [0.55 0.55 0.4 0.4];
        subplotposition(3,1:4) = [0.065 0.11 0.4 0.4];
        subplotposition(4,1:4) = [0.55 0.11 0.4 0.4];
        legendposition = 2;
        % preset color codes
        color(1,:) = [0,0,0];
        color(2,:) = [0.30 0.50 0.92]; %[0.04 0.14 0.42];
        color(3,:) = [0.168 0.505 0.337];
        color(4,:) = [0.847 0.16 0];
        % Create figure
        figure1 = figure('Color',[1 1 1]); 
        for m = 1:numel(M);
            axes1 = axes('Parent',figure1,'FontName','Arial',...
                'Position',subplotposition(m,1:4));
            % 'XTickLabel','', (remove setting it off
            ylim(axes1,[0 1.1]); xlim(axes1,[0 xmax]); hold(axes1,'all');
            errorbar1 = errorbar(Graph.X.(M{m}),Graph.Y.(M{m}),...
                Graph.E.(M{m}),...
                'Marker','.','LineWidth',2);
            ylabel(M{m},'FontName','Arial'); % Create ylabel
            if numel(groupname) <=4
                for g = 1:numel(groupname);
                    set(errorbar1(g),'DisplayName',gnshow{g},...
                            'LineWidth',2,'Color',color(g,1:3),...
                            'MarkerSize',20,'Marker','.'); 
                end
            elseif numel(groupname) >=5
                for g = 1:4;
                    set(errorbar1(g),'DisplayName',gnshow{g},...
                            'LineWidth',2,'Color',color(g,1:3),...
                            'MarkerSize',20,'Marker','.'); 
                end
                for g = 5:numel(groupname);
                    set(errorbar1(g),'DisplayName',gnshow{g},...
                            'LineWidth',2,...
                            'MarkerSize',20,'Marker','.'); 
                end  
            end
                
                
            if m ==legendposition; % if making righttop cornor
                for g = 1:numel(groupname);
                    %set(errorbar1(g),'DisplayName',gnshow{g},...
                        %'LineWidth',2);
    
                    legend1 = legend(axes1,'show');
                    set(legend1,'EdgeColor',[1 1 1],'YColor',[1 1 1],...
                        'XColor',[1 1 1],'TickDir','in',...
                        'LineWidth',1);
                end
            end
        end
        
        % create textbox for N
        for g = 1:numel(groupname); gN(g) = size(MWTfG.(groupname{g}),1); end
        n = num2str(gN'); b = cellfunexpr(groupname',' ');
        a = char(cell2mat(cellstr([n,char(b)])'));
        v = a(1,1:end); 
        t = 'N = '; 
        N = [t,v];
        annotation(figure1,'textbox',[0.003 0.02 0.20 0.05],'String',{N},...
            'FontName','Arial','FitBoxToText','off','EdgeColor','none');
        
        % save figure 
%         h = (gcf);
        titlename = ['ShaneSpark_CombineGraph']; % set name of the figure
        savefigjpeg(titlename,pSaveA);
%         set(h,'PaperPositionMode','auto'); % set to save as appeared on screen
%         cd(pSaveA);
%         print (h,'-dtiff', '-r0', titlename); saveas(h,titlename,'fig');
        % finish up
%         display 'Graph done.';
%         close;
        %% STEP6C: SAVE MATLAB
        cd(pSaveA); save('matlab.mat');  
        O.pSaveA = pSaveA;
    
    case 'DrunkPosture';
        %% EXCLUDE CHOR PROBLEM MWT FILES FROM ANALYSIS
        % RE-CHECK CHOR OUTPUTS
        display ' '; display 'Double checking chor outputs...'
        % prepare pMWTf input for validation
        A = celltakeout(struct2cell(MWTfG),'multirow');
        pMWTfT = A(:,2); MWTfnT = A(:,1);
        pMWTf = pMWTfT; MWTfn = MWTfnT;
        % check chor ouptputs
        pMWTfC = {}; MWTfnC = {}; 
        for v = 1:size(fnvalidate,1);
            [valname] = cellfunexpr(pMWTf,fnvalidate{v});
            [fn,path] = cellfun(@dircontentext,pMWTf,valname,'UniformOutput',0);
            novalfn = cellfun(@isempty,fn);
            if sum(novalfn)~=0;
                pMWTfnoval = pMWTf(novalfn); MWTfnoval = MWTfn(novalfn);
                pMWTfC = [pMWTfC;pMWTfnoval]; MWTfnC = [MWTfnC;MWTfnoval];
            end
        end
        pMWTfC = unique(pMWTfC); MWTfnC = unique(MWTfnC);
        % reporting
        if isempty(pMWTfC)==1; display 'All files have required Chor outputs';
        elseif isempty(pMWTfC)==0;
            str = 'Chore unsuccessful in %d MWT files';
            display(sprintf(str,numel(pMWTfC)));disp(MWTfnC);
            % STEP3F: EXCLUDE PROBLEM MWT FILES
            display 'Excluding problem MWT from analysis';
            gname = fieldnames(MWTfG);
            for x = 1:numel(MWTfnC) % each problem folders
                for g = 1:numel(gname)
                    A = MWTfG.(gname{g})(:,1);
                    i = logical(celltakeout(regexp(A,MWTfnC{x}),'singlenumber'));
                    if sum(i)>0; 
                        str = '>removing [%s]';
                        display(sprintf(str,MWTfnC{x}));
                        MWTfG.(gname{g})(i,:)=[]; % remove that from analysis
                    end
                end
            end
        end
        %% IMPORT AND PROCESS CHOR DATA
        display 'Importing drunkposture.dat';
        A = celltakeout(struct2cell(MWTfG),'multirow');
        pMWTfT = A(:,2); MWTfnT = A(:,1);
        pMWTf = pMWTfT; MWTfn = MWTfnT;
        % tnNslwakb
        % drunkposturedatL = {1,'time';2,'number';3,'goodnumber';4,'speed';5,'length';...
        %     6,'width';7,'aspect';8,'kink';9,'bias'};
        for p = 1 : numel(pMWTfT);
            str = 'Importing from [%s]';
            display(sprintf(str,MWTfnT{p}));
            [~,datevanimport] = dircontentext(pMWTfT{p},'*drunkposture.dat');  
            A(p,1) = MWTfnT(p);
            A(p,2) = {dlmread(datevanimport{1})};
        end
        MWTfdrunkposturedat = A;
        display 'Got all the drunkposture.dat.';
        %% PREPARE TIME POINTS
        % prepare universal timepoints limits
        % timepoints
        % find smallest starting time and smallest ending time
        MWTfdrunkposturedatL = {1,'MWTfname';2,'Data';3,'time(min)';4,'time(max)'};
        Raw = MWTfdrunkposturedat;    
        for p = 1:numel(MWTfn)
        Raw{p,3} = min(Raw{p,2}(:,1)); 
        Raw{p,4} = max(Raw{p,2}(:,1));
        end
        valstartime = max(cell2mat(Raw(:,3)));
        valendtime = min(cell2mat(Raw(:,4)));
        str = 'Earliest time tracked (MinTracked): %0.1f';
        display(sprintf(str,valstartime));
        str = 'Max time tracked  (MaxTracked): %0.1f';
        display(sprintf(str,valendtime));
        % processing inputs
        if tinput ==0; ti = valstartime; 
        elseif isempty(tinput)==1; ti = valstartime; tinput = 0; 
        elseif tinput>valstartime; ti = tinput; 
        end
        if isempty(intinput)==1; int = 10; else int = intinput; end
        if isempty(tfinput)==1; tf = valendtime; else tf = tfinput; end
        if isempty(durinput)==0; duration = 'restricted'; else duration = 'all'; end

        % reporting
        str = 'Starting time: %0.0fs';
        display(sprintf(str,ti));
        switch duration
        case 'all'
        timepoints = [0,ti+int:int:tf];
        str = 'Time points: %0.0f ';
        timeN = numel(timepoints);
        display(sprintf(str,timeN));
        case 'restricted'
        display 'Under construction';% need coding
        end       
        %% STATS
        % Stats.MWTfdrunkposturedat
        % drunkposturedatL = {1,'time';2,'number';3,'goodnumber';4,'Speed';5,'Length';...
        %     6,'Width';7,'Aspect';8,'Kink';9,'Bias'};
        Raw = MWTfdrunkposturedat;
        Graph = [];
        for p = 1:numel(MWTfn);
        Graph.X = timepoints; Graph.MWTfn = MWTfn';
            % summary 
            for t = 1:numel(timepoints)-1; % for each stim
                % get timeframe
                k = Raw{p,2}(:,1)>timepoints(t) & Raw{p,2}(:,1)<timepoints(t+1); 
                dataVal = Raw{p,2}(k,:);
                % create Graph.N
                Nrev = size(dataVal(:,2),1);
                Graph.N.Ndatapoints(t,p) = Nrev;
                Graph.N.NsumN(t,p) = sum(dataVal(:,2));
                Graph.N.NsumNVal(t,p) = sum(dataVal(:,3));
                Graph.Y.Speed(t,p) = mean(dataVal(:,4));
                Graph.E.Speed(t,p) = std(dataVal(:,4))./sqrt(Nrev);
                Graph.Y.Length(t,p) = mean(dataVal(:,5));
                Graph.E.Length(t,p) = std(dataVal(:,5))./sqrt(Nrev);
                Graph.Y.Width(t,p) = mean(dataVal(:,6));
                Graph.E.Width(t,p) = std(dataVal(:,6))./sqrt(Nrev);
                Graph.Y.Aspect(t,p) = mean(dataVal(:,7));
                Graph.E.Aspect(t,p) = std(dataVal(:,7))./sqrt(Nrev);        
                Graph.Y.Kink(t,p) = mean(dataVal(:,8));
                Graph.E.Kink(t,p) = std(dataVal(:,8))./sqrt(Nrev);         
                Graph.Y.Bias(t,p) = mean(dataVal(:,8));
                Graph.E.Bias(t,p) = std(dataVal(:,8))./sqrt(Nrev);
            end
        end
        Graph.YLegend = fieldnames(Graph.Y);
        clearvars Y X E;
        MWTfnimport = (Graph.MWTfn');
        M = Graph.YLegend;
        gnameL = gnamechose;       
        %% GRAPHING
        for m = 1:numel(M);% for each measure
            for g = 1:numel(gnameL);
                gname = gnameL{g};
                pMWTf = MWTfG.(gname)(:,2); 
                MWTfn1 = MWTfG.(gname)(:,1);
                A.MWTfn = MWTfn1;
                [~,i,~] = intersect(MWTfnimport(:,1),MWTfn1);
                Y(:,g) = mean(Graph.Y.(M{m})(:,i),2);
                E(:,g) = std(Graph.Y.(M{m})(:,i)')'./sqrt(sum(i));
                X(:,g) = Graph.X(2:end);
            end        
            % Create figure
            figure1 = figure;
            axes1 = axes('Parent',figure1);
            box(axes1,'on');
            hold(axes1,'all');
            errorbar1 = errorbar(X,Y,E);
            % create legend        
            gnshow = regexprep(gnameL,'_',' ');
            for g = 1:numel(gnshow)
                set(errorbar1(g),'DisplayName',gnshow{g});
            end
            legend1 = legend(axes1,'show');
            set(legend1,'Location','NorthEastOutside');            
            ylabel(M{m},'FontName','Arial','FontSize',30); % Create ylabel
            figname = [M{m},'[',num2str(ti,'%.0f'),':',num2str(int,'%.0f'),':',num2str(tf,'%.0f'),']'];
            savefig(figname,pSaveA);
        end  
        % STEP6C: SAVE MATLAB
        cd(pSaveA); save('matlab.mat');
        %% STATS & GRAPH: Standardize to 0mM
% Stats.MWTfdrunkposturedat
drunkposturedatL = {1,'time';2,'number';3,'goodnumber';4,'Speed';5,'Length';...
                    6,'Width';7,'Aspect';8,'Kink';9,'Bias'};
Raw = MWTfdrunkposturedat;
Graph = [];
Sum = [];
% for each plate
for p = 1:numel(MWTfn);
    Graph.X = timepoints; Graph.MWTfn = MWTfn';
    % for each stim
    for t = 1:numel(timepoints)-1; % for each stim
        % get timeframe
        k = Raw{p,2}(:,1)>timepoints(t) & Raw{p,2}(:,1)<=timepoints(t+1); 
        dataVal = Raw{p,2}(k,:);        
        % create summary data
        Sum{p,t} = dataVal;
    end
end


% get legend
iT = 1; iN = 2; igoodN = 3; iSpeed = 4; iLength = 5; iWidth = 6;
iAspect = 7; iKink =8; iBias = 9;
datL = {'time';'number';'goodnumber';'Speed';'Length';...
                    'Width';'Aspect';'Kink';'Bias'};

% create controls pair
%i = celltakeout(regexp(gnameL,'_'),'singlenumber');
%ctrlN = gnameL(i==0); 
ctrlN = gnameL(1); 
%ctrl_N2 = celltakeout(regexp(gnameL,'\<N2\>'),'singlenumber'); % N2 control
%expN = gnameL(i~=0);
expN = gnameL(2:end);

% find pairs
    Exp = []; Graph = []; 
for j = 1:numel(ctrlN)
                      
    % calculate control mean
    ctrl = ctrlN{j};
    m = MWTfG.(ctrl)(:,1);
    i = ismember(MWTfn,m);
    A = Sum(i,:);
    Control = []; B = [];
    for a = 1:size(datL,1) % for each analysis
        for t = 1:numel(timepoints)-1; % for each stim
            for p = 1:size(A,1)   
                dataVal = A{p,t};  % get timeframe data
                B.(datL{a})(t,p) = mean(dataVal(:,a));
            end
        end
        Control.(datL{a}) = mean(B.(datL{a}),2); 
    end

    % standardize to control
    a = celltakeout(regexp(expN,'_400mM','split'),'split');
    i = ismember(a(:,1),ctrlN{j});
    exp = expN{i};
    m = MWTfG.(exp)(:,1);
    i = ismember(MWTfn,m);
    A = Sum(i,:);

    for a = 1:size(datL,1) % for each analysis
        for t = 1:numel(timepoints)-1; % for each stim
            tc = Control.(datL{a})(t); % get control mean
            for p = 1:size(A,1)   
                b = (A{p,t}(:,a))./tc.*100;
                Exp.(datL{a})(t,p) = mean(b); % get data
                %Graph.N.(datL{a})(t,p) = size(A{p,t}(:,iN),1);
            end
        end
        d = Exp.(datL{a});
        n = sum(i);
        Graph.Y.(datL{a})(:,j) = mean(d,2);
        Graph.E.(datL{a})(:,j) = (std(d')/sqrt(n))';
        Graph.X.(datL{a})(:,j) = timepoints(2:end)';
    end

end

% graphing 
X = []; Y = []; E = [];
for m = 1:numel(M);% for each measure
    g = numel(expN);
    Y = Graph.Y.(M{m});
    E = Graph.E.(M{m});
    X = Graph.X.(M{m});
    % Create figure
    figure1 = figure;
    axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Calibri');
    hold(axes1,'all');
    
    % create errorbar
    errorbar1 = errorbar(X,Y,E);       
    gnshow = regexprep(expN,'_',' ');  
    % color codes
    color(1,:) = [0,0,0];
    color(2,:) = [0.04 0.14 0.42];
    color(3,:) = [0.847 0.16 0];
    color(4,:) = [0.168 0.505 0.337];
    for g = 1:numel(expN)
        if strcmp(gnshow{g},'N2 400mM')==1;
             set(errorbar1(g),'DisplayName',gnshow{g},...
            'MarkerSize',30,'Marker','.','LineWidth',2.5,'Color',[0 0 0]);
        elseif g == 2||3||4;
            set(errorbar1(g),'DisplayName',gnshow{g},...
            'MarkerSize',30,'Marker','.','LineWidth',2.5,'Color',color(g,1:3));
        end
    end
    
    
    
    % Create legend
    legend1 = legend(axes1,'show');
    set(legend1,'EdgeColor',[1 1 1],'Location','NorthEastOutside',...
    'YColor',[1 1 1],...
    'XColor',[1 1 1],...
    'FontSize',14);  
    
    % Create xlabel
    xlabel('Time (s)','FontSize',16,'FontName','Calibri');
    
    % Create ylabel
    ylabel(M{m},'FontName','Arial','FontSize',30); 
    
    % create 100% line
    Xmax = max(X); Xmax = Xmax(1);
    plot(repmat(100,1,round(Xmax)+1),'Parent',axes1,'LineWidth',3,'LineStyle',':',...
    'DisplayName','100','Color',[0 0 0]);
    
% save figure
    figname = [M{m},'_std[',num2str(ti,'%.0f'),':',num2str(int,'%.0f'),':',num2str(tf,'%.0f'),']'];
    savefigjpeg(figname,pSaveA);
end
        % STEP6C: SAVE MATLAB
        cd(pSaveA); save('matlab.mat');


    case 'LadyGaGa'
        %% STEP5: IMPORT CHOR OUTPUTS
        % STEP5A: Import and process Chor generated data
        % MWTfevandat
        % note: this take a long time
        display 'Importing evan.dat';
        % evandat legends (%tnNss*b12)
        evandatL = {1,'time';2,'number';3,'goodnumber';4,'speed';5,'speedStd';...
        6,'bias';7,'tap';8,'puff'};
        MWTfevandatL = {1,'MWT filename';2,'Data';3,'mintime';4,'maxtime';...
        5,'mean number tracked';6,'mean good number';7,'mean percent N valid';...
        8,'mean speed';9,'min speed';10,'max speed';11,'mean bias'};
        % import
        MWTfevandat = {}; 
        for p = 1:numel(pMWTf); % for each mwt
        [~,datevanimport] = dircontentext(pMWTf{p},'*evan.dat');  
        MWTfevandat(p,1) = MWTfn(p);
        MWTfevandat(p,2) = {dlmread(datevanimport{1})};
        end
        % summarize
        for x = 1:size(MWTfevandat,1)
        MWTfevandat{x,3} = min(MWTfevandat{x,2}(:,1)); % min time
        MWTfevandat{x,4} = max(MWTfevandat{x,2}(:,1)); % max time
        MWTfevandat{x,5} = mean(MWTfevandat{x,2}(:,2)); % mean number tracked
        MWTfevandat{x,6} = mean(MWTfevandat{x,2}(:,3)); 
        MWTfevandat{x,7} = MWTfevandat{x,6}/MWTfevandat{x,5}*100; 
        MWTfevandat{x,8} = mean(MWTfevandat{x,2}(:,4)); 
        MWTfevandat{x,9} = min(MWTfevandat{x,2}(:,4)); 
        MWTfevandat{x,10} = max(MWTfevandat{x,2}(:,4)); 
        MWTfevandat{x,11} = mean(MWTfevandat{x,2}(:,6)); 
        end
        display 'done.';
        % sprevs
        display 'Importing *.sprevs, this will take a while...';
        % sprevs legends
        MWTfsprevsL = {1,'MWT filename'; 2,'Data';3,'mintime';4,'maxtime';...
        5,'objects tracked';6,'mean reversal';7,'minreversal'; ...
        8,'max reversal'};
        sprevsL = {1,'object ID'; 2,'time'; 3,'reversal'; 4,'trackduration?'};
        % import data
        MWTfsprevs = {};
        for p = 1 : numel(pMWTf); 
        [~,ps] = dircontentext(pMWTf{p},'*.sprevs');
        A = []; for k = 1:size(ps,1); a = dlmread(ps{k}); A = [A;a]; end
        MWTfsprevs(p,1) = MWTfn(p); MWTfsprevs(p,2) = {A};
        end 
        % summarize sprevs
        for x = 1:size(MWTfsprevs,1)
        MWTfsprevs{x,3} = min(MWTfsprevs{x,2}(:,2));
        MWTfsprevs{x,4} = max(MWTfsprevs{x,2}(:,2)); 
        MWTfsprevs{x,5} = numel(unique(MWTfsprevs{x,2}(:,1))); 
        MWTfsprevs{x,6} = mean(MWTfsprevs{x,2}(:,3)); 
        MWTfsprevs{x,7} = min(MWTfsprevs{x,2}(:,3)); 
        MWTfsprevs{x,8} = max(MWTfsprevs{x,2}(:,3)); 
        end
        display 'done';
        % summary Data
        Data.MWTfsprevs = MWTfsprevs;
        Data.MWTfsprevsL= MWTfsprevsL;
        Data.MWTfevandat = MWTfevandat;
        Data.MWTfevandatL = MWTfevandatL;
        Data.evandatL = evandatL;
        Data.sprevsL = sprevsL;
        %% STEP6: DESCRIPTIVE STATS
        % STEP6A: PREPARE TIME POINTS
        % prepare groups
        gnameL = fieldnames(MWTfG);
        % prepare universal timepoints limits
        % find smallest starting time and smallest ending time
        t = max(cell2mat(Data.MWTfsprevs(:,3)));
        t = str2num(num2str(t,'%.0f'));
        MWTfsprevStartTime = t;
        t = min(cell2mat(Data.MWTfsprevs(:,4)));
        t = str2num(num2str(t,'%.0f'));
        MWTfsprevEndTime = t;
        str = 'Earliest time tracked (MinTracked): %0.1f';
        display(sprintf(str,MWTfsprevStartTime));
        str = 'Max time tracked  (MaxTracked): %0.1f';
        display(sprintf(str,MWTfsprevEndTime));
        % processing inputs
        if tinput ==0; ti = MWTfsprevStartTime; 
        elseif isempty(tinput)==1; ti = MWTfsprevStartTime; tinput = 0; 
        elseif tinput>MWTfsprevStartTime; ti=tinput; end
        if isempty(intinput)==1; int = 10; else int = intinput; end
        if isempty(tfinput)==1; tf = MWTfsprevEndTime; else tf = tfinput; end
        if isempty(durinput)==0; duration = 'restricted'; else duration = 'all'; end
        % reporting
        str = 'Starting time: %0.0fs';
        display(sprintf(str,ti));
        switch duration
        case 'all'
        timepoints = [ti:int:tf];
        str = 'Time points: %0.0f ';
        timeN = numel(timepoints);
        display(sprintf(str,timeN));
        case 'restricted'
        % need coding
        end

        % STEP6B: DESCRIPTIVE STATS BY GROUPS
        % summarize group data
        G = [];
        for g = 1:numel(gnameL);
        % get group info
        A = [];
        gname = gnameL{g};    
        pMWTf = MWTfG.(gname)(:,2);
        MWTfn = MWTfG.(gname)(:,1);
        A.MWTfn = MWTfn;
        [~,i,~] = intersect(Data.MWTfsprevs(:,1),MWTfn);
        MWTfsprevs = Data.MWTfsprevs(i,:);
        [~,i,~] = intersect(Data.MWTfevandat(:,1),MWTfn);
        MWTfevandat = Data.MWTfevandat(i,:);

        %% MWTfsprevsum summary
        for p = 1:numel(MWTfn); % for each MWT plate
        A.X(:,p) = timepoints';
        for t = 1:numel(timepoints)-1; % for each stim
            % get data durint time frame 
            k = MWTfsprevs{p,2}(:,2)>timepoints(t) ....
                & MWTfsprevs{p,2}(:,2)<timepoints(t+1);
            %%
            sprevsValid = MWTfsprevs{p,2}(k,:);
            m = MWTfevandat{p,2}(:,1)>timepoints(t) ...
                & MWTfevandat{p,1}(:,1)<timepoints(t+1);
            evandatValid = MWTfevandat{p,2}(m,:);
            %% get N.Min
            if isempty(evandatValid)==1; % if no valid times
                A.N.Minimum(t,p) = 0; 
            else A.N.Minimum(t,p) = min(evandatValid(:,3)); 
            end
            % get N.Max
            if isempty(evandatValid)==1; 
                A.N.Maximum(t,p) = 0; 
            else A.N.Maximum(t,p) = max(evandatValid(:,3)); 
            end
            Nrev = size(sprevsValid(:,2),1);
            A.Y.RevIncidents(t,p) = Nrev;
            A.Y.RevWorm(t,p) = size(unique(sprevsValid(:,1)),1);
            A.Y.RevDist(t,p) = mean(sprevsValid(:,3));
            A.E.RevDist(t,p) = std(sprevsValid(:,3))./sqrt(Nrev);
            A.Y.RevDur(t,p) = mean(sprevsValid(:,4));
            A.E.RevDur(t,p) = std(sprevsValid(:,4))./sqrt(Nrev);
            % calculate total duration reversed
            RevEndTime = sprevsValid(:,2)+sprevsValid(:,4); % time start+time end
            overTimei = RevEndTime(:,1)>timepoints(t+1); % time tracked till after the end of timepoint
            RevEndTime(overTimei,1) = timepoints(t+1);
            RevDur = RevEndTime(:,1)-sprevsValid(:,2);
            A.Y.TotalTimeRev(t,p) = sum(RevDur);
        end % end of loop timepoints
        end % end of loop MWT
        % summarize group data
        Ytype = fieldnames(A.Y);
        B.X(:,g) = mean(A.X,2);
        for y = 1:numel(Ytype)
        Stat = Ytype{y};
        B.Y.(Stat)(:,g) = mean(A.Y.(Stat),2);
        B.E.(Stat)(:,g) = (std(A.Y.(Stat)')')./sqrt(size(A.Y.(Stat),2));
        end
        B.MWTf.(gname) = MWTfn;
        end % end of loop for group
        Graph = B;
        Graph.GroupName = gnameL;
        %% STEP7: GRAPH GROUPED DATA
        MeasureType = fieldnames(Graph.Y);
        for a = 1:numel(MeasureType);
        X = Graph.X(2:end,:);
        Y = Graph.Y.(MeasureType{a});
        E = Graph.E.(MeasureType{a});
        errorbar(X,Y,E);
        figname = [(MeasureType{a}),'[',num2str(ti,'%.0f'),':',num2str(int,'%.0f'),':',num2str(tf,'%.0f'),']'];
        savefig(figname,pSaveA);
        end  
        %% STEP6C: SAVE MATLAB
        cd(pSaveA); save('matlab.mat');
    otherwise
end
end



%% SUBFUNCTIONS


    
    
    
  

%% OLD CODE
% switch ID
%     case 'ShaneSpark'
%         %% standard habituation curve anaylysis
%         %% set paths
%         pExp = pInput;
%         ShaneSpark(pExp,pFun);
% 
%     case 'LadyGaGa'
%         %% analyze reversal responses without taps
%         choice = 'suspend'; % coding status
%         graph = 'combined'; % graph choice
%         % ask if have done chor analysis
%         display('Was LadyGaGa choreography analysis completed previously?');
%         Overwrite = input('[y=1,n=0]: ');
%         display('Zip current file as backup?');
%         zipfile = input('[y=1,n=0]: ');
%         display('New run condition?');
%         sprevsQ = input('[y=1,n=0]: ');
%         % set path
%         pExp = pInput;
%         % run program
%         LadyGaGa(pExp,pFun,choice,graph,Overwrite,zipfile,sprevsQ);
% 
% 
%     case 'Igor'        
%         %% gives analysis of raw chor output
%         %% add chor analysis to IgoreRaw2
%         pExp = pInput;
%         Igor(pExp,pFun); 
% 
%     case 'combineExp'
%         %% [unfinished script] combine experiment
%         combineexp;
%         
%         
%     case 'Org1'
%     %% Step1: organize files, standard
%     MWTorgfiles2(pFun,pExp); 
% 
% 
%     case 'Hab1'
%     IgoreRaw3(pExp,pFun);% ungrouped data, standard habitutaion curve
% 
%     case 'Raw1' % grouped data, MWTfdata analysis only (to be removed soon)
%     Igor3Raw1v2(pExp,pFun);
%     [Time,Graph,Stats] = igoreGraphing2(pFun,pExp,'time','import*.mat'); 
% 
% 
%     case 'Raw2' % ungrouped data, MWTfdata analysis only
%     % [BUG] fixing no renaming on IgoreRaw2v2
%     IgoreRaw2(pExp,pFun); 
% 
%     case 'G1' % graphing, ungrouped data, MWTfdata analysis only
%     IgoreGraph1(pExp,pFun);
%     case 'G2' % graphing with standard habituation curve
%     IgoreGraph2(pExp,pFun,'import*.mat');
%     %case 'G3'
%     %   igoreGraphing3(pFun,pExp,Xn,rawfilename,intmin,int,intmax,m);
%     case 'C1' % combine data already analyzed by Raw
%     [MWTfdatGG,ExpGA] = getMWTfdatcombined(pBet);
%     case 'GC1'; % make all graphs and save it in pBet
%     [~,~,~,intmin,int,intmax,Xn] = selectmeasure2(pFun,...
%         'analysisetting.mat',2);
%     for m = 1:15;
%         [Stats,Graph,Data] = igoreGraphCombineM(pFun,pBet,MWTfdatGG,intmin,int,intmax,m);        
%     end
% 
%     case 'GdG1'; % get graph data within groups, use data analyzed by C1
%     [GraphData,Time,Xaxis,Yn,titlename] = IgoreGraphDataexpwithingroup(MWTfdatGG,g,intmin,int,intmax,m,pFun);
% 
%     case 'Gexp1'
%     makefigwithinexp(GraphData,g,setfilename,titlename,pFun,Yn,pBet);
% 
%     otherwise
%     error('Invalid analysis code "%s"',ID);
% end
% end

% %% STATS & GRAPH: Standardize to 0mM
% % Stats.MWTfdrunkposturedat
% drunkposturedatL = {1,'time';2,'number';3,'goodnumber';4,'Speed';5,'Length';...
%     6,'Width';7,'Aspect';8,'Kink';9,'Bias'};
% Raw = MWTfdrunkposturedat;
% Graph = [];
% Sum = [];
% % for each plate
% for p = 1:numel(MWTfn);
% Graph.X = timepoints; Graph.MWTfn = MWTfn';
% % for each stim
% for t = 1:numel(timepoints)-1; % for each stim
% % get timeframe
% k = Raw{p,2}(:,1)>timepoints(t) & Raw{p,2}(:,1)<=timepoints(t+1); 
% dataVal = Raw{p,2}(k,:);        
% % create summary data
% Sum{p,t} = dataVal;
% end
% end
% 
% 
% % get legend
% iT = 1; iN = 2; igoodN = 3; iSpeed = 4; iLength = 5; iWidth = 6;
% iAspect = 7; iKink =8; iBias = 9;
% datL = {'time';'number';'goodnumber';'Speed';'Length';...
%     'Width';'Aspect';'Kink';'Bias'};
% 
% % create controls pair
% i = celltakeout(regexp(gnameL,'_'),'singlenumber');
% ctrlN = gnameL(i==0); % 0mM control = no _
% ctrl_N2 = celltakeout(regexp(gnameL,'\<N2\>'),'singlenumber'); % N2 control
% expN = gnameL(i~=0);
% % find pairs
% Exp = []; Graph = []; 
% for j = 1:numel(ctrlN)
% 
% % calculate control mean
% ctrl = ctrlN{j};
% m = MWTfG.(ctrl)(:,1);
% i = ismember(MWTfn,m);
% A = Sum(i,:);
% Control = []; B = [];
% for a = 1:size(datL,1) % for each analysis
% for t = 1:numel(timepoints)-1; % for each stim
% for p = 1:size(A,1)   
% dataVal = A{p,t};  % get timeframe data
% B.(datL{a})(t,p) = mean(dataVal(:,a));
% end
% end
% Control.(datL{a}) = mean(B.(datL{a}),2); 
% end
% 
% % standardize to control
% a = celltakeout(regexp(expN,'_400mM','split'),'split');
% i = ismember(a(:,1),ctrlN{j});
% exp = expN{i};
% m = MWTfG.(exp)(:,1);
% i = ismember(MWTfn,m);
% A = Sum(i,:);
% 
% for a = 1:size(datL,1) % for each analysis
% for t = 1:numel(timepoints)-1; % for each stim
% tc = Control.(datL{a})(t); % get control mean
% for p = 1:size(A,1)   
% b = (A{p,t}(:,a))./tc.*100;
% Exp.(datL{a})(t,p) = mean(b); % get data
% %Graph.N.(datL{a})(t,p) = size(A{p,t}(:,iN),1);
% end
% end
% d = Exp.(datL{a});
% n = sum(i);
% Graph.Y.(datL{a})(:,j) = mean(d,2);
% Graph.E.(datL{a})(:,j) = (std(d')/sqrt(n))';
% Graph.X.(datL{a})(:,j) = timepoints(2:end)';
% end
% 
% end
% 
% % graphing 
% X = []; Y = []; E = [];
% for m = 1:numel(M);% for each measure
% g = numel(expN);
% Y = Graph.Y.(M{m});
% E = Graph.E.(M{m});
% X = Graph.X.(M{m});
% % Create figure
% figure1 = figure;
% axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Calibri');
% hold(axes1,'all');
% 
% % create errorbar
% errorbar1 = errorbar(X,Y,E);       
% gnshow = regexprep(expN,'_',' ');  
% % color codes
% color(1,:) = [0,0,0];
% color(2,:) = [0.04 0.14 0.42];
% color(3,:) = [0.847 0.16 0];
% color(4,:) = [0.168 0.505 0.337];
% for g = 1:numel(expN)
% if strcmp(gnshow{g},'N2 400mM')==1;
% set(errorbar1(g),'DisplayName',gnshow{g},...
% 'MarkerSize',30,'Marker','.','LineWidth',2.5,'Color',[0 0 0]);
% elseif g == 2||3||4;
% set(errorbar1(g),'DisplayName',gnshow{g},...
% 'MarkerSize',30,'Marker','.','LineWidth',2.5,'Color',color(g,1:3));
% end
% end
% 
% 
% 
% % Create legend
% legend1 = legend(axes1,'show');
% set(legend1,'EdgeColor',[1 1 1],'Location','NorthEastOutside',...
% 'YColor',[1 1 1],...
% 'XColor',[1 1 1],...
% 'FontSize',14);  
% 
% % Create xlabel
% xlabel('Time (s)','FontSize',16,'FontName','Calibri');
% 
% % Create ylabel
% ylabel(M{m},'FontName','Arial','FontSize',30); 
% 
% % create 100% line
% Xmax = max(X); Xmax = Xmax(1);
% plot(repmat(100,1,round(Xmax)+1),'Parent',axes1,'LineWidth',3,'LineStyle',':',...
% 'DisplayName','100','Color',[0 0 0]);
% 
% % save figure
% figname = [M{m},'_std[',num2str(ti,'%.0f'),':',num2str(int,'%.0f'),':',num2str(tf,'%.0f'),']'];
% savefigjpeg(figname,pSaveA);
% end
% 
% %% STEP6C: SAVE MATLAB
% cd(pSaveA); save('matlab.mat');
% 










