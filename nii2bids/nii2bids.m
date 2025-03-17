%
clear all;
%%
baseDir = '/Users/dennis.jungchildmind.org/OneDrive - Child Mind Institute/bolt2025/bidsformat/dataset_chang';
bidsDir = fullfile(baseDir,'bidsit');%output directory where bids format will be applied 
dataFolders = {'anat','func'};%,'physio'};%folders for organization
taskName = strcat('_task-rest_');
bidsPostfix = {'T1w','bold'};%,'physio'};%
subfolder = 'raw';%if there is a subfolder, if not remove it
sessionPrefixRegex = 'mr_[0-9]*';%session name prefix is "mr" in this case
%create a top-level folder
if ~exist(bidsDir,'dir'), mkdir(bidsDir); end
%add dataset_description.json
if ~exist(fullfile(bidsDir,'dataset_description.json'),'file')
    fid = fopen(fullfile(bidsDir,'dataset_description.json'),'w');
    jsonInput = struct('BIDSVersion',"1.0.0","Name","Dataset_chang","DatasetDOI","Something");
    encodedJSON = jsonencode(jsonInput);
    fprintf(fid,encodedJSON);
    fclose('all');
end
%add participants.tsv
if ~exist(fullfile(bidsDir,'participants.tsv'),'file'); 
    %generate participants.tsv, if I get error.
end

%let's setup file names and sort what the heck is going on
allFiles = [];

for i = 1:length(dataFolders)
    %grab filenames from the specified folders
    temp = dir(fullfile(baseDir,dataFolders{i},'raw'));
    temp = {temp.name};
    fInd= logical(cell2mat(cellfun(@(x) length(regexp(x,'.*sub_.*')),temp,'uni',0)));
    allFiles{i} = temp(fInd);
end

%assume that there is sMRI for each session (I mean this is fMRI/EEG data
%so anatomical (sMRI) data should exist ('anat' folder)
anatInd = find(cell2mat(cellfun(@(x) regexp(x,'anat'), dataFolders,'uni',0)));
%extract subject ID from the folder names 
subjectID = cellfun(@(x) regexp(x,'sub_[0-9]*','match'),allFiles{anatInd});
subjectIDbids = cellfun(@(x) strrep(x,'_','-'),subjectID,'uni',0);
%check if there are any dupes
if sum(cell2mat(cellfun(@(x) cellfun(@(y) strcmpi(x,y), subjectID),subjectID,'uni',0)'),'all') ~= length(subjectID)
    error('Duplicate subjects exist');
end

% iterate each subject
for i = 1:length(subjectID)
    %new subject folder will be name in "sub-" with a hyphen
    subjectFolder = fullfile(bidsDir,subjectIDbids{i});
    %generate a session folder (in Bolt et al. 2024 data, it will be "mr") before data folders
    checkSessionPrefix = [];
    %%
    for k = 1:length(dataFolders)
        folders2check = fullfile(baseDir,dataFolders{k},subfolder);
        folders2check = dir(folders2check);
        files2check = {folders2check.name};
        %check (1) number of files associated with the subject in each data
        %folder (2) and availablity of the data file for the particular subject
        fileCheck = cellfun(@(x) length(regexp(x,subjectID{i})), files2check);
        
        if sum(fileCheck) > 0 
            %create session folder in the subjectFolder
            %(1) check the naming convention and session number
            availFiles = files2check(logical(fileCheck));
            checkSessionPrefix{k} = cellfun(@(x) regexp(x,sessionPrefixRegex,'match'), availFiles,'uni',0);
            %check if the session is the same, if so just reduce it to 1
            A = cellfun(@(x) regexp(x,sessionPrefixRegex,'match'), availFiles,'uni',0);
            B = [A{:}];
            checkMat = cell2mat(cellfun(@(x) cellfun(@(y) strcmpi(x,y), B)', B,'uni',0));

            if ~isempty(checkMat)
                if sum(checkMat,'all') == length(availFiles).^2
                    %they are all same sessions, although multiple data file exist
                    checkSessionPrefix{k} = checkSessionPrefix{k}(1);
                else
                    %more than 2 sessions
                    nSess = unique(length(availFiles)./sum(checkMat));
                    nSpacing = sum(checkMat,'all')./length(availFiles);
                    fprintf('%d sessions exist.\n',nSess);
                    checkSessionPrefix{k} = checkSessionPrefix{k}(1:nSpacing:length(checkSessionPrefix{k}));
                end
            else
                disp('wha')
            end
            %*should write elseif for the case (e.g., 3 files, 2 files from 1 session and 1 file from another session).
        end   
    end
    %% check if there are any files associated with anatomy file
    %if only anatomy file, the length of checksessionPrefix should be 1
    
    if length(checkSessionPrefix) == 1
        %if the length equals to 1, just go to next file
        continue;
    end

    %% then check if I should create the session folder here
    [~,I] = max(cellfun(@(y) sum(cellfun(@length,y)),checkSessionPrefix));%grab the one with the most session names (since we need to create based on that?)
    
    sessionNumbers = checkSessionPrefix{I};
    sessionNumbers = [sessionNumbers{:}];
    sessionNumbers= cellfun(@(x) regexp(x,'[0-9]*','match'), sessionNumbers);

    %make it BIDS format
    sessionFolders2create = cellfun(@(x) strcat('ses-',x),sessionNumbers,'uni',0);

    if isempty(sessionFolders2create)
        continue;
        %('Session folder for creation is empty...!');
    end
    
    %create session folders 
    for sf = 1:length(sessionFolders2create)
        if ~exist(subjectFolder,'dir'), mkdir(subjectFolder); end%create the subject folder first 
        %then add the session folders
        tempFolder = fullfile(subjectFolder,sessionFolders2create{sf});
        if ~exist(tempFolder,'dir'), mkdir(tempFolder); end%create it

        %for each session folder, generate data folders ('anat','func',...) based on the dataFolders
        try
            if exist(tempFolder,'dir')
                cellfun(@(x) mkdir(x), cellfun(@(x) fullfile(tempFolder,x), dataFolders, 'uni', 0));
            end
        catch
            warning('Data folders cannot be generated.');
            return;
        end

        %we not have folders organized based on the subjects and data types
        %we will move and reorganize files to the appropriate folders.
        for ii = 1:length(dataFolders)
            moveFrom = fullfile(baseDir,dataFolders{ii},subfolder);%where the files are at originally
            moveTo = fullfile(tempFolder,dataFolders{ii});%where the files should go
            fileList = dir(moveFrom);
            fileList = {fileList.name};
            fileInd = find(cellfun(@(x) length(regexp(x,subjectID{i})),fileList));
            if length(fileInd) > 1
                disp('stop ')
            end
            if isempty(fileInd), warning('No file associated with the subjectID exists (data type: %s, subject ID: %s)', dataFolders{ii}, subjectID{i}); continue; end
            
            %anatomical file (make exception to this for now since there
            %are one anat for other data files).

            files2move = fileList(fileInd);%file to move

            
            if ~strcmpi('anat',dataFolders{ii})
                ind2move = find(cellfun(@(x) length(regexp(x,sessionNumbers{sf})), files2move));
                files2move = files2move(ind2move);
            else
                files2move = files2move(1);%if "anat" just the first file (for this dataset specifically)
            end
            

            %just check echo number using lookbehind assertion (?<=echo)
            echoNum =  cellfun(@(x) regexp(x,'(?<=echo)[0-9]*','match'), files2move,'uni',0);
            fileExtension = cell2mat(regexp(files2move{1},'(\.mat|\.nii\.gz)$','match'));

            if logical(sum(cellfun(@isempty,echoNum)))
                newFileName = {strcat(subjectIDbids{i},'_',sessionFolders2create{sf},taskName,bidsPostfix{ii})};
            else
                %add echo number here
                newFileName = cellfun(@(x) strcat(subjectIDbids{i},'_',sessionFolders2create{sf},taskName,'echo-',num2str(x{1}),'_',bidsPostfix{ii}),...
                    echoNum,'uni',0);
            end
            
            dataFileName = cellfun(@(x) strcat(x,fileExtension), newFileName, 'uni',0);
            jsonFileName = cellfun(@(x) strcat(x,'.json'), newFileName, 'uni',0);
            
            %(1) save the data
            cellfun(@(x,y) copyfile(fullfile(moveFrom,x),fullfile(moveTo,y)), files2move,dataFileName, 'uni',0);%move the file from the original directory to destination
            %(2) then save JSON sidecar corresponding to the saved data
            cellfun(@(x) saveJsonFile(moveTo,x), jsonFileName);
        end
    
    end
    
end


%
function saveJsonFile(moveTo,jsonFileName)
if ~exist(fullfile(moveTo,jsonFileName),'file')
    fid = fopen(fullfile(moveTo,jsonFileName),'w');
    S = struct();%json file
    S.TaskName = cell2mat(regexp(jsonFileName,'(?<=task-).[a-zA-Z]*','match'));
    S.Manufacturer = "Siemens";
    S.ManufacturersModelName = "Avanto";
    S.MagneticFieldStrength = 3;
    S.RepetitionTime = 2.1;
    echoNum = cell2mat(regexp(jsonFileName,'(?<=echo-).[0-9]*','match'));
    
    if ~isempty(echoNum)
        echoNum = str2num(echoNum);
    end
        
    echoTime = 1;
    if echoNum == 1
        echoTime = 0.013;%second
    elseif echoNum == 2
        echoTime = 0.0294;%second
    elseif echoNum == 3
        echoTime = 0.0457;%second
    else
        disp("No echo time, likely not BOLD")
    end

    S.EchoTime = echoTime;
    S.FlipAngle = 75;
    S.PhaseEncodingDirection = "j-";
    S.liceEncodingDirection = "k";
    S.SliceTiming = [0.07, 0.21, 0.35, 0.49, 0.63, 0.77, 0.91, 1.05, 1.19,...
                     1.33, 1.47, 1.61, 1.75, 1.89, 2.03, 0, 0.14, 0.28, 0.42,...
                     0.56, 0.70, 0.84, 0.98, 1.12, 1.26, 1.4, 1.54, 1.68, 1.82, 1.96];

    encodedJSON = jsonencode(S);
    fprintf(fid,encodedJSON);
    fclose('all');
end
end