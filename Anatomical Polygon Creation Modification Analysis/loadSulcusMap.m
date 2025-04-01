function vertices = loadSulcusMap(SulcusMap_fileLocation,SessionRunList)

%% Load sulcus map for selected sessions and runs


sulcusMap_folderContents = dir(SulcusMap_fileLocation);
%Removing possible invalid entries returned by dir command
sulcusMap_folderContents = sulcusMap_folderContents(~cellfun('isempty', {sulcusMap_folderContents.date}));
%Selecting only files.
sulcusMap_folderContents = sulcusMap_folderContents([sulcusMap_folderContents.isdir]~=1);

%Shows variable structure. Not used for anything.
FileInfoLabels = {'Session Number' 'Run' 'Date' 'LatestMap'};

%parse filenames into session information
for file = 1:length(sulcusMap_folderContents)
    try
        sulcusMap_fileInfo(file,1) = str2num(sulcusMap_folderContents(file).name(strfind(sulcusMap_folderContents(file).name,'Sess')+4:strfind(sulcusMap_folderContents(file).name,'Run')-1));
        sulcusMap_fileInfo(file,2) = str2num(sulcusMap_folderContents(file).name(strfind(sulcusMap_folderContents(file).name,'Run')+3:strfind(sulcusMap_folderContents(file).name,'_Sulcus')-1));
        sulcusMap_fileInfo(file,3) = datenum(sulcusMap_folderContents(file).name(strfind(sulcusMap_folderContents(file).name,'Map')+4:end-4));
    catch
        warning([sulcusMap_folderContents(file).name ' does not follow designated naming structure.'])
    end
end

sulcusMap_fileInfo(:,4) = 0;

%Only keep latest sulcus map for each session
uniqueSessions = unique(sulcusMap_fileInfo(:,1));
for i = 1:length(uniqueSessions)
    uniqueRuns = unique(sulcusMap_fileInfo(sulcusMap_fileInfo(:,1)==uniqueSessions(i),2));
    for j = 1:length(uniqueRuns)
        latestROI = max(sulcusMap_fileInfo(sulcusMap_fileInfo(:,1)==uniqueSessions(i) & sulcusMap_fileInfo(:,2)==uniqueRuns(j),3));
        useSessionInd = sulcusMap_fileInfo(:,1)==uniqueSessions(i) & ...
            sulcusMap_fileInfo(:,2)==uniqueRuns(j) & ...
            sulcusMap_fileInfo(:,3)==latestROI;
        sulcusMap_fileInfo(useSessionInd,4)=1;
    end
end
%Load sulcus map
for i = 1:length(SessionRunList(:,1))
    loadInd = sulcusMap_fileInfo(:,1)==SessionRunList(i,1) & ...
        sulcusMap_fileInfo(:,2)==SessionRunList(i,2) & ...
        sulcusMap_fileInfo(:,4)==1;
    load(fullfile(SulcusMap_fileLocation, sulcusMap_folderContents(loadInd).name),'VerticePositions')
    vertices{i} = VerticePositions;
end

end