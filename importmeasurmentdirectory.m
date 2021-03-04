function [time, timeMS, tempcell] = importmeasurmentdirectory(dirpath)
    temp = regexp(dir(dirpath),'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match');
    time = str2num(temp) * 0.001;
    timeMS = str2num(temp);

    fileList = dir(dirpath + '/*.txt');
    for i = 1:length(fileList)
        tempcell{i} = importfiledata(fileList(i).name);
    end
    time, timeMS, tempcell
end