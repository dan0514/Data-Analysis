function [O] = Dance_20140204r
%% INSTRUCTION
% modified from Dance to suit Flame
% Lined modified by Daniel: 27, 142, 143



%% STEP1: GET PATHS 
restoredefaultpath;
% add program path (Mac only)
pProgram = cd;
p = pProgram;
a = {}; % create cell array for output
a = dir(cd); % list content
a = {a.name}'; % extract folder names only
a(ismember(a,{'.','..','.DS_Store'})) = []; 
for x = 1:size(a,1); % for all files 
    p1 = [p,'/',a{x,1}]; % make path for files
    if isdir(p1) ==1; % if a path is a folder
        addpath(p1);
    end
end


%% Get experiment home folders
display ' ';
display 'Select experiment home folder';
[fn,p] = dircontent(fileparts(p));
% get rid of .xxx folders
i = cellfun(@isempty,regexp(fn,'[.]|(AnalysisProgram)'));
fn = fn(i); p = p(i);
display(makedisplay(fn));
display ' ';
i = input('folder: ');
pData = p{i};

% get save folder path
savefoldername = 'AnalysisResults';
pSave = [pData,'/',savefoldername ];
if exist(pSave,'dir') ~=7
    mkdir(pData,savefoldername);
    display 'made analysis folder';
end



% Get data paths
% display '1) drag data folder into Commmand window';
% display '2) copy and paste the path portion and press enter';
% display 'i.e. cd('Volumnes/Conny/20120223C_CL_100s30x10s10s')';
% display 'path portion is Volumnes/Conny/20120223C_CL_100s30x10s10s';
% display ' ';
% display 'Getting data folder paths';
% display 'drag experiment folder and press [Enter]';
% i = input(': ','s');
% p = cd

% a = regexp(userpath,':','split'); cd(a{:,1}); PathCommonList; clearvars a;
% addpath(genpath(paths.pFunMWT));
% pMaestro = '/Volumes/AWLIGHT/MWT/Archived Analysis/Maestro';
% pRosePick = '/Volumes/Rose/MWT_Analysis_20130811';
% pRose = '/Volumes/Rose/MWT_Analysis_20131020_Liquid';
% pSet = '/Users/connylinlin/Documents/MATLAB/MATLAB MWT/MatSetfiles';
% pSum = '/Volumes/Rose/MultiWormTrackerPortal/Summary';
% pSave = paths.pSave;
% 
% 


%% STEP2: INTERPRET INPUT AND DEAL OUTPUT
HomePath = pData;
O = []; % [CODE] Define output O % temperary assign O = nothing
% STEP2: OPTION INTERPRETATION
%% STEP2A: CREATE OPTION LIST 
display ' ';
display 'Select analysis option...';
ID = {'ShaneSpark';'LadyGaGa';'DrunkPosture'};
display(makedisplay(ID))
display ' ';
option = ID{input('option: ')};


%% STEP2B: CREATE OPTION NOT NEEDING TIME INPUTS
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
% [release only by experiment search] search for target paths
display ' ';
% display 'Select data by...'
% a = {'Experiment name';'MWT name';'Group name'};
% b = {'E'; 'M'; 'G'};
% display(makedisplay(a));
% display ' ';
% i = input('search method: ');
% searchclass = b{i};
searchclass = 'E';
display 'Enter search term for experiment folder name';


% 
% display 'MWT files [M], Exp files [E], Group name [G] or no search [A]?';
% searchclass = input(': ','s');
switch searchclass
    case 'M' % search for MWT
        display 'option under construction'
    case 'G' % search for group name
        display 'option under construction'
        
    case 'E' % search for Experiment folders
        display 'Choose a search term: (from folder names)';
        display 'eg. 20140221B_DH_100s30x10s10s or DH ';
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
display ' ';
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
display ' ';
display 'Name your analysis output folder';
a = clock;
time = [num2str(a(4)),num2str(a(5))];
name = [input(': ','s'),'_',option,'_',datestr(now,'yyyymmdd'),time];
mkdir(pSave,name);
pSaveA = [pSave,'/',name];



%% STEP5: CHOR 
% STEP4B: CHECK IF CHOR HAD BEEN DONE
% display ' '; display 'Checking chor outputs...'
% % prepare pMWTf input for validation
A = celltakeout(struct2cell(MWTfG),'multirow');
pMWTfT = A(:,2); MWTfnT = A(:,1); pMWTf = pMWTfT; MWTfn = MWTfnT;
% % check chor ouptputs
% pMWTfC = {}; MWTfnC = {}; 
% for v = 1:size(fnvalidate,1);
%     [valname] = cellfunexpr(pMWTf,fnvalidate{v});
%     [fn,path] = cellfun(@dircontentext,pMWTf,valname,'UniformOutput',0);
%     novalfn = cellfun(@isempty,fn);
%     if sum(novalfn)~=0;
%         pMWTfnoval = pMWTf(novalfn); MWTfnoval = MWTfn(novalfn);
%         pMWTfC = [pMWTfC;pMWTfnoval]; MWTfnC = [MWTfnC;MWTfnoval];
%     end
% end
% pMWTfC = unique(pMWTfC); MWTfnC = unique(MWTfnC);
% % reporting
% if isempty(pMWTfC)==0;
%     str = 'Need to Chore %d MWT files';
%     display(sprintf(str,numel(pMWTfC)));
% else display 'All files have required Chor outputs';
% end
% STEP4C: RUN CHORE
% if isempty(pMWTfC)==0; 
    [~] = chormaster_20140204(pMWTf,pProgram,chorcode);
% else
%     display 'No Chor needs to be ran.';
% end




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
        display 'Importing drunkposture.dat...';
        A = celltakeout(struct2cell(MWTfG),'multirow');
        pMWTfT = A(:,2); MWTfnT = A(:,1);
        pMWTf = pMWTfT; MWTfn = MWTfnT;
        % tnNslwakb
        % drunkposturedatL = {1,'time';2,'number';3,'goodnumber';4,'speed';5,'length';...
        %     6,'width';7,'aspect';8,'kink';9,'bias'};
        for p = 1 : numel(pMWTfT);
            %str = 'Importing from [%s]';
            %display(sprintf(str,MWTfnT{p}));
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
            savefigjpeg(figname,pSaveA);
        end  
        % STEP6C: SAVE MATLAB
        cd(pSaveA); save('matlab.mat');
        
%         %% STATS & GRAPH: Standardize to 0mM
% % Stats.MWTfdrunkposturedat
% drunkposturedatL = {1,'time';2,'number';3,'goodnumber';4,'Speed';5,'Length';...
%                     6,'Width';7,'Aspect';8,'Kink';9,'Bias'};
% Raw = MWTfdrunkposturedat;
% Graph = [];
% Sum = [];
% % for each plate
% for p = 1:numel(MWTfn);
%     Graph.X = timepoints; Graph.MWTfn = MWTfn';
%     % for each stim
%     for t = 1:numel(timepoints)-1; % for each stim
%         % get timeframe
%         k = Raw{p,2}(:,1)>timepoints(t) & Raw{p,2}(:,1)<=timepoints(t+1); 
%         dataVal = Raw{p,2}(k,:);        
%         % create summary data
%         Sum{p,t} = dataVal;
%     end
% end
% 
% 
% % get legend
% iT = 1; iN = 2; igoodN = 3; iSpeed = 4; iLength = 5; iWidth = 6;
% iAspect = 7; iKink =8; iBias = 9;
% datL = {'time';'number';'goodnumber';'Speed';'Length';...
%                     'Width';'Aspect';'Kink';'Bias'};
% 
% % create controls pair
% %i = celltakeout(regexp(gnameL,'_'),'singlenumber');
% %ctrlN = gnameL(i==0); 
% ctrlN = gnameL(1); 
% %ctrl_N2 = celltakeout(regexp(gnameL,'\<N2\>'),'singlenumber'); % N2 control
% %expN = gnameL(i~=0);
% expN = gnameL(2:end);
% 
% % find pairs
%     Exp = []; Graph = []; 
% for j = 1:numel(ctrlN)
%                       
%     % calculate control mean
%     ctrl = ctrlN{j};
%     m = MWTfG.(ctrl)(:,1);
%     i = ismember(MWTfn,m);
%     A = Sum(i,:);
%     Control = []; B = [];
%     for a = 1:size(datL,1) % for each analysis
%         for t = 1:numel(timepoints)-1; % for each stim
%             for p = 1:size(A,1)   
%                 dataVal = A{p,t};  % get timeframe data
%                 B.(datL{a})(t,p) = mean(dataVal(:,a));
%             end
%         end
%         Control.(datL{a}) = mean(B.(datL{a}),2); 
%     end
% 
%     % standardize to control
%     a = celltakeout(regexp(expN,'_400mM','split'),'split');
%     i = ismember(a(:,1),ctrlN{j});
%     exp = expN{i};
%     m = MWTfG.(exp)(:,1);
%     i = ismember(MWTfn,m);
%     A = Sum(i,:);
% 
%     for a = 1:size(datL,1) % for each analysis
%         for t = 1:numel(timepoints)-1; % for each stim
%             tc = Control.(datL{a})(t); % get control mean
%             for p = 1:size(A,1)   
%                 b = (A{p,t}(:,a))./tc.*100;
%                 Exp.(datL{a})(t,p) = mean(b); % get data
%                 %Graph.N.(datL{a})(t,p) = size(A{p,t}(:,iN),1);
%             end
%         end
%         d = Exp.(datL{a});
%         n = sum(i);
%         Graph.Y.(datL{a})(:,j) = mean(d,2);
%         Graph.E.(datL{a})(:,j) = (std(d')/sqrt(n))';
%         Graph.X.(datL{a})(:,j) = timepoints(2:end)';
%     end
% 
% end
% 
% % graphing 
% X = []; Y = []; E = [];
% for m = 1:numel(M);% for each measure
%     g = numel(expN);
%     Y = Graph.Y.(M{m});
%     E = Graph.E.(M{m});
%     X = Graph.X.(M{m});
%     % Create figure
%     figure1 = figure;
%     axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Calibri');
%     hold(axes1,'all');
%     
%     % create errorbar
%     errorbar1 = errorbar(X,Y,E);       
%     gnshow = regexprep(expN,'_',' ');  
%     % color codes
%     color(1,:) = [0,0,0];
%     color(2,:) = [0.04 0.14 0.42];
%     color(3,:) = [0.847 0.16 0];
%     color(4,:) = [0.168 0.505 0.337];
%     for g = 1:numel(expN)
%         if strcmp(gnshow{g},'N2 400mM')==1;
%              set(errorbar1(g),'DisplayName',gnshow{g},...
%             'MarkerSize',30,'Marker','.','LineWidth',2.5,'Color',[0 0 0]);
%         elseif g == 2||3||4;
%             set(errorbar1(g),'DisplayName',gnshow{g},...
%             'MarkerSize',30,'Marker','.','LineWidth',2.5,'Color',color(g,1:3));
%         end
%     end
%     
%     
%     
%     % Create legend
%     legend1 = legend(axes1,'show');
%     set(legend1,'EdgeColor',[1 1 1],'Location','NorthEastOutside',...
%     'YColor',[1 1 1],...
%     'XColor',[1 1 1],...
%     'FontSize',14);  
%     
%     % Create xlabel
%     xlabel('Time (s)','FontSize',16,'FontName','Calibri');
%     
%     % Create ylabel
%     ylabel(M{m},'FontName','Arial','FontSize',30); 
%     
%     % create 100% line
%     Xmax = max(X); Xmax = Xmax(1);
%     plot(repmat(100,1,round(Xmax)+1),'Parent',axes1,'LineWidth',3,'LineStyle',':',...
%     'DisplayName','100','Color',[0 0 0]);
%     
% % save figure
%     figname = [M{m},'_std[',num2str(ti,'%.0f'),':',num2str(int,'%.0f'),':',num2str(tf,'%.0f'),']'];
%     savefigjpeg(figname,pSaveA);
% end
%         % STEP6C: SAVE MATLAB
%         cd(pSaveA); save('matlab.mat');


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
        f1 = figure('Visible','on'); 
        hold on;
        errorbar(X,Y,E);
        figname = [(MeasureType{a}),'[',num2str(ti,'%.0f'),':',num2str(int,'%.0f'),':',num2str(tf,'%.0f'),']'];
        %saveas(gcf,figname,'fig'); % save as matlab figure 
        savefigjpeg(figname,pSaveA);
        end  
        %% STEP6C: SAVE MATLAB
        cd(pSaveA); save('matlab.mat');
    otherwise
end


%% STEP7: REPORT AND END
display 'Analysis completed';

end










%% SUBFUNCTION
function [a,b,c,d] = dircontent(p)
% a = full dir, can be folder or files, b = path to all files, 
% c = only folders, d = paths to only folders
cd(p); % go to directory
a = {}; % create cell array for output
a = dir; % list content
a = {a.name}'; % extract folder names only
a(ismember(a,{'.','..','.DS_Store'})) = []; 
b = {};
c = {};
d = {};
for x = 1:size(a,1); % for all files 
    b(x,1) = {strcat(p,'/',a{x,1})}; % make path for files
    if isdir(b{x,1}) ==1; % if a path is a folder
        c(end+1,1) = a(x,1); % record name in cell array b
        d(end+1,1) = b(x,1); % create paths for folders
    else
    end
end
end



function [show] = makedisplay(varargin)
% option = 'bracket', 1) list
% option = string, i.e. '-', then 1-list; 2-list
%%
% A = varargin{1}; 
% A

list = varargin{1};

if numel(varargin) >=2
    A = varargin{2};
else
    A = 'none';    
end

switch A
    case 'bracket'
        numb = num2str((1:numel(list))');
        [b] = cellfunexpr(list,')');
        show = [numb,char(b),char(list)];
    case 'none'
        numb = num2str((1:numel(list))');
        [b] = cellfunexpr(list,')');
        show = [numb,char(b),char(list)];
end
    
end


function [Output] = MWTDataBaseMaster(HomePath,option)
switch option
    case 'GetStdMWTDataBase'
        display 'Searching for Experiment folders...';
        [~,~,fn,p] = dircontent(HomePath);
        expnameidentifider = '^\d{8}[A-Z][_]([A-Z]{2})[_](.){1,}';
        i = not(cellfun(@isempty,regexp(fn,expnameidentifider,'match')));
        pNonExpf = p(not(i)); NonExpfn = fn(not(i));
        pExpfD = p(i); ExpfnD = fn(i);
        str = '%d [%d standardized] experiment folders';
        display(sprintf(str,numel(fn),sum(i)));
        Output.pExpfD = pExpfD; Output.ExpfnD = ExpfnD;
        Output.NonExpfn = NonExpfn; Output.pNonExpf = pNonExpf;


        % pGfD and GfnD
        display 'Searching for Group folders...';
        [~,~,fn,p] = cellfun(@dircontent,pExpfD,'UniformOutput',0);
        fn = celltakeout(fn,'multirow');
        i = not(celltakeout(regexp(fn,'MatlabAnalysis'),'singlenumber'));
        Gfn = fn(i);
        p = celltakeout(p,'multirow');
        pGf = p(i);
        [fn,p,~,~] = cellfun(@dircontent,pGf,'UniformOutput',0);
        empty = cellfun(@isempty,fn); % see which group folder is empty
        pGfD = pGf(not(empty));
        GfnD = Gfn(not(empty));
        if sum(empty)>1; 
            pGfproblem = pGf(empty); 
            display '> The following folders are empty:';
            disp(Gfn(empty));
        end
        str = '> %d folders found under Exp folders';
        display(sprintf(str,numel(Gfn)));
        str = '> %d [%d unique] Group folders';
        display(sprintf(str,numel(GfnD),numel(unique(GfnD))));
        Output.GfnD = GfnD; Output.pGfD = pGfD;
        
        % pMWTfD & MWTfnD
        display 'Searching for MWT folders...';
        fn = celltakeout(fn(not(empty)),'multirow');
        p = celltakeout(p(not(empty)),'multirow');
        mwt = logical(celltakeout(regexp(fn,'\<(\d{8})[_](\d{6})\>'),'singlenumber'));
        pMWTf = p(mwt); MWTfn = fn(mwt);
        str = '> %d [%d unique] MWT folders';
        display(sprintf(str,numel(MWTfn),numel(unique(MWTfn))));
        Output.pMWTf = pMWTf; Output.MWTfn = MWTfn;
        
        % Zip files?
        display 'Searching for zipped files...';
        zip = logical(celltakeout(regexp(fn,'.zip'),'singlenumber'));
        if sum(zip)>1; display 'zipped files found'; 
            pZipf = p(zip); Zipfn = fn(zip);
            disp(Zipfn); 
            str = '> %d [%d unique] zip files';
            display(sprintf(str,numel(Zipfn),numel(unique(Zipfn))));
            Output.pZipf = pZipf; Output.Zipfn = Zipfn;
        else
            display '> No zip files found.';
        end
    case 'GetExpTargetInfo'
        pExpfD = HomePath;
        % pGfD and GfnD
        display 'Searching for Group folders';
        [~,~,fn,p] = cellfun(@dircontent,pExpfD,'UniformOutput',0);
        fn = celltakeout(fn,'multirow');
        i = not(celltakeout(regexp(fn,'MatlabAnalysis'),'singlenumber'));
        Gfn = fn(i);
        p = celltakeout(p,'multirow');
        pGf = p(i);
        [fn,p,~,~] = cellfun(@dircontent,pGf,'UniformOutput',0);
        empty = cellfun(@isempty,fn); % see which group folder is empty
        pGfD = pGf(not(empty));
        GfnD = Gfn(not(empty));
        if sum(empty)>1; 
            pGfproblem = pGf(empty); 
            display 'the following folders are empty:';
            disp(Gfn(empty));
        end
        str = '%d folders found under Exp folders';
        display(sprintf(str,numel(Gfn)));
        str = '%d [%d unique] Group folders';
        display(sprintf(str,numel(GfnD),numel(unique(GfnD))));
        Output.GfnT = GfnD; Output.pGfT = pGfD;
        
        % pMWTfD & MWTfnD
        display 'Searching for MWT folders';
        fn = celltakeout(fn(not(empty)),'multirow');
        p = celltakeout(p(not(empty)),'multirow');
        mwt = logical(celltakeout(regexp(fn,'\<(\d{8})[_](\d{6})\>'),'singlenumber'));
        pMWTf = p(mwt); MWTfn = fn(mwt);
        str = '%d [%d unique] MWT folders';
        display(sprintf(str,numel(MWTfn),numel(unique(MWTfn))));
        Output.pMWTfT = pMWTf; Output.MWTfnT = MWTfn;
        
        % Zip files?
        display 'Searching for zipped files';
        zip = logical(celltakeout(regexp(fn,'.zip'),'singlenumber'));
        if sum(zip)>1; display 'zipped files found'; 
            pZipf = p(zip); Zipfn = fn(zip);
            disp(Zipfn); 
            str = '%d [%d unique] zip files';
            display(sprintf(str,numel(Zipfn),numel(unique(Zipfn))));
            Output.pZipf = pZipf; Output.Zipfn = Zipfn;
        else
            display 'No zip files found.';
        end
    case 'GetSingleExpInfo'
        pExpfD = HomePath;
        % pGfD and GfnD
        display 'Searching for Group folders';
        [~,~,fn,p] = dircontent(pExpfD);
        fn = celltakeout(fn,'multirow');
        i = not(celltakeout(regexp(fn,'MatlabAnalysis'),'singlenumber'));
        Gfn = fn(i);
        p = celltakeout(p,'multirow');
        pGf = p(i);
        [fn,p,~,~] = cellfun(@dircontent,pGf,'UniformOutput',0);
        empty = cellfun(@isempty,fn); % see which group folder is empty
        pGfD = pGf(not(empty));
        GfnD = Gfn(not(empty));
        if sum(empty)>1; 
            pGfproblem = pGf(empty); 
            display 'the following folders are empty:';
            disp(Gfn(empty));
        end
        str = '%d folders found under Exp folders';
        display(sprintf(str,numel(Gfn)));
        str = '%d [%d unique] Group folders';
        display(sprintf(str,numel(GfnD),numel(unique(GfnD))));
        Output.GfnT = GfnD; Output.pGfT = pGfD;
        
        % pMWTfD & MWTfnD
        display 'Searching for MWT folders';
        fn = celltakeout(fn(not(empty)),'multirow');
        p = celltakeout(p(not(empty)),'multirow');
        mwt = logical(celltakeout(regexp(fn,'\<(\d{8})[_](\d{6})\>'),'singlenumber'));
        pMWTf = p(mwt); MWTfn = fn(mwt);
        str = '%d [%d unique] MWT folders';
        display(sprintf(str,numel(MWTfn),numel(unique(MWTfn))));
        Output.pMWTfT = pMWTf; Output.MWTfnT = MWTfn;
        
        % Zip files?
        display 'Searching for zipped files';
        zip = logical(celltakeout(regexp(fn,'.zip'),'singlenumber'));
        if sum(zip)>1; display 'zipped files found'; 
            pZipf = p(zip); Zipfn = fn(zip);
            disp(Zipfn); 
            str = '%d [%d unique] zip files';
            display(sprintf(str,numel(Zipfn),numel(unique(Zipfn))));
            Output.pZipf = pZipf; Output.Zipfn = Zipfn;
        else
            display 'No zip files found.';
        end

    case 'FindAllMWT'
        % code based on
            % [b] = getalldir(home); 
            % [Output] = dircontentmwtall(HomePath);
        %display 'Searching for MWT in all drives, this will take a while...';         
        a = cellfun(@genpath,HomePath,'UniformOutput',0);
        paths = regexpcellout(a,pathsep,'split');
        paths(cellfun(@length,paths)<1) = []; % get rid of cell has zero lengh
        paths = paths';
        [~,fn] = cellfun(@fileparts,paths, 'UniformOutput',0);
        search = '(\<(\d{8})[_](\d{6})\>)';
        k = regexpcellout(fn,search);
        Output.pMWTf = paths(k);
        Output.MWTfn = fn(k);
  
        
        
        
    otherwise
end
end


function [exprname] = cellfunexpr(cell2match,expr)
    exprname(1:size(cell2match,1),1) = {expr};
end



function [A] = celltakeout(cell,option)

switch option
    case 'split'
        A = {};
        for x = 1:numel(cell); col = size(cell{x},2); 
            A(x,1:col) = cell{x}; end
    case 'multirow'
        A = {}; for x = 1:numel(cell); A = [A;cell{x}]; end
    case 'singlerow'
        A = {};
        for x = 1:numel(cell);
            if isempty(cell{x})==0; A(x,1) = cell{x,1};
            else A(x,1) = {''}; end
        end
    case'match'
        A = {};
        for x = 1:numel(cell);
            if isempty(cell{x})==0; A(x,1) = cell{x,1};
            else A(x,1) = {''}; end
        end
    case 'singlenumber'
        A = [];
        for x = 1:numel(cell)
            if isempty(cell{x})==0; A(x,1) = cell{x,1};
            else A(x,1) = 0; end
        end
    case 'logical'
         A = [];
        for x = 1:numel(cell)
            if isempty(cell{x})==0; A(x,1) = cell{x,1};
            else A(x,1) = 0; end
        end
        A = logical(A);
        
    otherwise
end

%

end

function [MWTf,pMWTf] = dircontentmwt(pExp)

%% main function
% You might like this function I made called 'dircontentmwt'! 
% [MWTf,pMWTf] = dircontentmwt(pExp)
% - pExp = character variable of the path of your experiment folder
% - MWTf = a cell array containing the name of MWT folders
% - pMWTf = a cell array containing the paths of MWT folders
% 
% It gets you only MWT folders name and path from a folder path 'pExp'. In the function file attached you can see how I validate the name of the dir..
% Use it like any regular function in Matlab, like [path,folder] = fileparts(x). fileparts is the function, and x, path,folder are all variables you can define.
% For example, if your experiment folder is in current path you can just use 'cd' like this:  [MWTf,pMWTf] = dircontentmwt(cd)
% If you don't want to have both output, for example, if you just want pMWTf you can go [~,pMWTf] =  dircontentmwt(cd)
% The '~' ignores the MWTf output.
% 
% 
% The pMWTf is useful for loop cd into each MWT folders, like this:
% 
% for x = 1:numl(pMWTf); % for each validated MWT folders
% 	cd(pMWTf{x});  % path to MWT folder on cell array row x
% end



[~,~,f,p] = dircontent(pExp);
if isempty(f)==1;
%    MWTf = []; pMWTf = [];
else
    % check if there is mwt under group folder
    for x = 1:numel(p)
        [~,~,f1,p1] = dircontent(p{x});
        k(1:numel(p1),x) = cellfun(@isdir,p1);    
    end
end
if sum(sum(k)) > 0;
    Status = 'group folder';
else
    Status = 'MWT folder'; 
end


MWTf = []; pMWTf = [];
switch Status
    case 'group folder';
        for x = 1:numel(p)
            [~,~,f1,p1] = dircontent(p{x});      
            filestruct = cellfunexpr(f1,'\<(\d{8})[_](\d{6})\>');
            match = cellfunexpr(f1,'match');
            t = cellfun(@regexp,f1,filestruct,'UniformOutput',0);
            i = celltakeout(t,'logical'); 
            a = f1(i,1);
            b = p1(i,1); 
            MWTf = [MWTf;a];
            pMWTf = [pMWTf;b];

        end
        
    case 'MWT folder';
        % get mwt folders
        filestruct(1:size(f,1),1) = {'\<(\d{8})[_](\d{6})\>'};
        t = cellfun(@regexp,f,filestruct,'UniformOutput',0);
        i = find(cell2mat(t))';
        MWTf = f(i,1);
        pMWTf = p(i,1); 
end


end


function [a,b] = dircontentext(p,ext)
% a = full dir, can be folder or files
% b = path to all files
% if want all files under a certain foler called ext, then do not add * mark

cd(p); % go to directory
a = {}; % create cell array for output
a = dir(ext); % list content
a = {a.name}'; % extract folder names only
a(ismember(a,{'.','..','.DS_Store'})) = []; 
b = {};
for x = 1:size(a,1); % for all files 
    b(x,1) = {strcat(p,'/',a{x,1})}; % make path for files
end
end


function savefigjpeg(titlename,pSave)
% save figures 
% titlename = 'CombinedGraph';
cd(pSave);
h = (gcf);
set(h,'PaperPositionMode','auto'); % set to save as appeared on screen
print (h,'-djpeg', '-r150', titlename); % save as tiff
saveas(h,titlename,'fig'); % save as matlab figure 
close;
end



function [A] = regexpcellout(C,searchterm,varargin)
% function [A] = regexpcellout(C,searchterm,option)
% upated 20140204
optionlist = {'logical'};
reglist = {'split','match'};


% make sense of inputs
switch nargin
    case 2
        option = 'logical';
        % validate first input as cell
        %if iscell(varargin{1}), C = varargin{1}; end
        %if ischar(varargin{2}), searchterm = varargin{2}; end
        B = regexp(C,searchterm);
    
    case 3
        %if iscell(varargin{1}), C = varargin{1}; end
        %if ischar(varargin{2}), searchterm = varargin{2}; end
        a = varargin{1};
        if ischar(a); option = a; end
        i = strcmp(reglist,option);
        if sum(i)==1
            B = regexp(C,searchterm,option);
        else
            B = regexp(C,searchterm);
        end
        
    otherwise
        error 'incorect number of inputs';
end

  
switch option
    case 'split'
        A = {};
        for x = 1:numel(B); col = size(B{x},2); 
            A(x,1:col) = B{x}; end
    case'match'
        col = cell2mat(cellfun(@size,B,'UniformOutput',0));
        A = cell(numel(B),max(col(:,2)));
        for x = 1:numel(B);
            if isempty(B{x})==0; 
                col = size(B{x,1},2);
                A(x,1:col) = B{x};
            else
                A(x,1) = {''}; 
            end
        end
    case 'logical'
         A = [];
        for x = 1:numel(B)
            if isempty(B{x})==0; A(x,1) = B{x,1};
            else A(x,1) = 0; end
        end
        A = logical(A); 
%     case 'multirow'
%         A = {}; for x = 1:numel(B); A = [A;B{x}]; end
%     case 'singlerow'
%         A = {};
%         for x = 1:numel(B);
%             if isempty(B{x})==0; A(x,1) = B{x,1};
%             else A(x,1) = {''}; end
%         end
%     case 'singlenumber'
%         A = [];
%         for x = 1:numel(B)
%             if isempty(B{x})==0; A(x,1) = B{x,1};
%             else A(x,1) = 0; end
%         end       
    otherwise
end

end



function [chorscript] = chormaster_20140204(pMWT,pProgram,varargin)

% function [chorscript] = chormaster2(paths,option,pData)

% CHOREGRAPHY
% input must be a rx1 cell array containing rows of paths to MWT folders
% source code: MWT003_chormaster



Var = varargin;


% PATH
% % pData & pBadFiles
% i = cellfun(@isdir,Var); % see if pData path is provided
% 
% if sum(i) ==1;  % if pData is provided
%     pData = Var(i); 
% 
%     
% % if provided with pData and pBadFiles
% elseif sum(i) ==2; 
%     i = find(i);
%     pData = i(1);
%     pBadFiles = i(2);
% 
%     
% % if no path provided
% elseif sum(i) ==0;
%     
%     PathCommonList;
%     
%     % check if rose is a dir
%     p = paths.pRose;
%     if isdir(p) ==1; 
%         pData = paths.MWT.pRoseData; 
%         pBadFiles = paths.MWT.pBadFiles;
% 
%     else
%         error 'Rose is not attached';
%     end
%     
% end




% VALIDATE pMWT INTPUTS
% get all MWT folders under input paths
p = [];
for x = 1:numel(pMWT)
    if exist(pMWT{x},'dir')==7
        p2 = genpath(pMWT{x});
        p = [p,':',p2];
    end
end

a = regexp(p,':','split');
b = regexp(a,'\<(\d{8})[_](\d{6})\>');
i = celltakeout(b','singlenumber');
i = i~=0;
pMWT = a(i)';
if isempty(pMWT)
    error 'path contains no MWT folder';
end


% CHECK INTEGRITY OF MWT FILES
% prepare pMWTf input for validation
display 'Validating MWT folder contents...';
% check for files
fname = '*.blobs';
a = cellfun(@numel,(cellfun(@dircontentext,pMWT,...
    cellfunexpr(pMWT,fname),'UniformOutput',0))); % get numer of files
fname = '*.summary';
b = cellfun(@numel,(cellfun(@dircontentext,pMWT,...
    cellfunexpr(pMWT,fname),'UniformOutput',0))); % get numer of files
fname = '*.set';
c = cellfun(@numel,(cellfun(@dircontentext,pMWT,...
    cellfunexpr(pMWT,fname),'UniformOutput',0))); % get numer of files
% fname = '*.png';
% d = cellfun(@numel,(cellfun(@dircontentext,paths,...
%     cellfunexpr(paths,fname),'UniformOutput',0))); % get numer of files


% MOVE BAD MWT FILES TO BAD FOLDER 

[r,c] = find([a,b,c]==0);
if isempty(r)==0;
    % make pBadFile folder
    p = fileparts(fileparts(fileparts(fileparts(pMWT{1}))));
    filename = 'MWT_BadFiles';
    pBadFiles = [p,'/',filename];
    if exist(pBadFiles,'dir')~=7
        mkdir(p,filename)
    end
    pData = p;
    
    p = pMWT(r); % % find files with missing MWT files
    a = cellfun(@strrep,p,cellfunexpr(p,pData),cellfunexpr(p,pBadFiles),...
        'UniformOutput',0); % replace pData path with pBadFiles paths
    %makemisspath(cellfun(@fileparts,a,'UniformOutput',0)); % make folders
    cellfun(@movefile,p,a); % move files
    str = '[%d] bad MWT files moved out of Analysis folder';
    display(sprintf(str,numel(a)));
    pMWT = pMWT(~r); % update path
else
    display 'All MWT files contain .blob, .summary & .set';
end


% GET CHOR OPTION
AList = {'LadyGaGa';
        'ShaneSpark';
        'Beethoven';
        'DrunkPosture';
        'AnnaPavlova';
        'Rastor';
        'SwanLake'};


% get option from varargin
i = cellfun(@isstr,Var) & ~cellfun(@isdir,Var);
option = AList{regexpcellout(AList,Var(i))};


% %% chor options list
% %  ------[update as you add more options] ---- 
% AnaList = '(LadyGaGa)|(ShaneSpark)|(Beethoven)|(DrunkPosture)|(AnnaPavlova)|(Rastor)|(SwanLake)';
% %  ------[update as you add more options] ----
% 
% Chor = char(regexp(option,AnaList,'match'));
% if isempty(Chor)==1; % if no match, select options
%     % create display list
%     select = sortrows(AnaList);
%     [b1] = cellfunexpr(AnaList,'['); [b2] = cellfunexpr(AnaList,']');
%     num = cellstr(num2str((1:numel(AnaList))'));
%     display(char(cellfun(@strcat,b1,num,b2,select,'UniformOutput',0)));
%     % ask for analysis id
%     display ' '; answer = input('Select analysis or [Enter] to abort: ');
%     if isempty(answer)==1; return; end
%     Chor = select{answer};
% end

%% JAVA ARGUMENTS--------------------------------------------------
% path to java programs
%javapath = [strrep(userpath,pathsep,''),'/MATLAB MWT/SubFun_Java'];
javapath = [pProgram,'/Java'];

b = blanks(1); % blank
% call java 
javacall = 'java -jar'; javaRAM = '-Xmx8G'; javaRAM7G = '-Xmx7G';
beethoven = ['''',javapath,'/Beethoven_v2.jar','''']; % call beethoven 
chor = ['''',javapath,'/Chore_1.3.0.r1035.jar','''']; % call chor 
% chor calls 
map = '--map';
% settings 
pixelsize = '-p 0.027'; speed = '-s 0.1'; 
mintime = '-t 20'; minmove = '-M 2'; shape = '--shadowless -S';
% plugins 
preoutline = '--plugin Reoutline::exp';  
prespine = '--plugin Respine';
% plugins (reversals) 
revbeethoven_trv = '--plugin MeasureReversal::tap::dt=1::collect=0.5::postfix=trv';
revignortap_sprevs = '--plugin MeasureReversal::postfix=sprevs';
rev_ssr = '--plugin MeasureReversal::tap::collect=0.5::postfix=ssr';

% dat output 
odrunkposture = '-O drunkposture -o nNslwakb';
oconny = '-o 1nee#e*ss#s*SS#S*ll#l*LL#L*ww#w*aa#a*mm#m*MM#M*kk#k*bb#b*pp#p*dd#d'; % Conny's 
obeethoven = '-o nNss*b12M'; % standard for Beethoven
oshanespark = '-O shanespark -o nNss*b12M'; % standard for Beethoven
oevan = '-O evan -o nNss*b12'; % Evan's dat output
oevanall = '-O evanall -N all -o nNss*b12';
oswanlakeall = '-O swanlakeall -N all -o tnNemMawlkcspbd1';
oswanlake = '-O swanlake -o tnNemMawlkcspbd1e#m#M#a#w#l#k#c#s#p#b#d#e-m-M-a-w-l-k-c-s-p-b-d-e*m*M*a*w*l*kvc*s*p*b*d*';




%% CREATE JAVA SYNTAX (chorescript) ---------------------------------------
chorscript = {};
switch option
    case 'LadyGaGa'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oevan,b,preoutline,b,prespine,b,...
            revignortap_sprevs,b]; 
        fval = {'*evan.dat';'*.sprevs'};
    case 'DrunkPosture'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oevan,b,odrunkposture,b,....
            preoutline,b,prespine,b]; 
        fval = {'*drunkposture.dat'};
    case 'ShaneSpark'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oshanespark,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b]; 
        fval = {'*.trv';'*shanespark.dat'};
    case 'Beethoven'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,obeethoven,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b]; 
        fval = {'*.trv'};
    case 'AnnaPavlova'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oshanespark,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b]; 
        fval = {'*.trv';'*shanespark.dat'};
    case 'Rastor'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oevanall,b,preoutline,b,prespine,b];
        fval = {'*evanall*'};
    case 'SwanLake'
        chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oswanlakeall,b,preoutline,b,...
            prespine,b,revbeethoven_trv,b];   
        chorscript{2} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,oswanlake,b,preoutline,b,...
            prespine,b]; 
        fval = {'*.trv';'*swanlake.dat'; '*swanlakeall*'};
    otherwise
        error('option does not exist in chore master');
end
% validate chorescript is cell
if iscell(chorscript) ==0; error 'chorescript must be in cell array'; end



%% CHECK IF CHOR HAD BEEN DONE
% DEFINE CHOROUTPUTS AND CHOR CODE
display ' '; display 'Checking chor outputs...'
% check chor ouptputs
val = ones(numel(pMWT),numel(fval));
for v = 1:size(fval,1);
    [fn,~] = cellfun(@dircontentext,pMWT,cellfunexpr(pMWT,fval{v})...
        ,'UniformOutput',0); % search for files
    val(:,v) = cellfun(@isempty,fn);
end

% get files that need chor
i = find(sum(val,2));
if isempty(i) == 0; 
    pMWT = unique(pMWT(i)); 
    display(sprintf('Need to Chore %d MWT files',numel(pMWT)));
else display 'All files contain required Chor outputs'; return;end


%% STEP4C: RUN CHORE
% flexibility with different paths inputs
str = 'Chor-ing MWTfile [%s]...';
for x = 1:numel(pMWT); 
    [~,fn] = fileparts(pMWT{x}); file = strcat('''',pMWT{x},''''); 
    display(sprintf(str,fn));
    for cs = 1:numel(chorscript) 
        system([chorscript{cs} file], '-echo'); 
    end  
end
display 'Chor Finished.';
end % function end










