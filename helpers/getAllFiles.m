function fileList = getAllFiles(dirName, pattern, do_recursiv)
%getAllFiles list all filesnames matching a pattern recursively
% Example: filelistCSV = getAllFiles(RootPath,'\d+_\d+_\d+\.csv$');
% match 0123_20110101_20111231.csv

  dirData = dir(dirName);      %# Get the data for the current directory
  dirIndex = [dirData.isdir];  %# Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
    matchstart = regexp(fileList, pattern);
    fileList = fileList(~cellfun(@isempty, matchstart));
  end
  if exist('do_recursiv','var')
      if do_recursiv
          subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
          validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                                       %#   that are not '.' or '..'
          for iDir = find(validIndex)                  %# Loop over valid subdirectories
            nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
            fileList = [fileList; getAllFiles(nextDir, pattern, do_recursiv)];  %# Recursively call getAllFiles
          end
      end
  end

end
