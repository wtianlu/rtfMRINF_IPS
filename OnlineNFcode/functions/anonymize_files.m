function anonymize_files(mr8,ds,run)
fprintf('***anonymize_files.m***********************\n')
if mr8 
    disp('Move and rename DRIN files...')
    exportdrin = ds.dirs.exportdrin;
    drinfiles = dir([exportdrin filesep '*Dyn00*']);
    if isempty(drinfiles)
        drinfiles2 = dir([exportdrin filesep 'DRIN-00*']); % check SHORT format
        if ~isempty(drinfiles2)
            for ii = 1:length(drinfiles2)
                movefile([exportdrin filesep drinfiles2(ii).name],[ds.r(run).dir filesep 'tbv_watch' filesep ds.idstr drinfiles2(ii).name])
            end
            disp('Finished moving DRIN files')
        else
            disp('No DRIN files found!')
        end
    else
        idxcut = cellfun(@length,strsplit(drinfiles(1).name,'-'));
        for ii = 1:length(drinfiles)
            movefile([exportdrin filesep drinfiles(ii).name],[ds.r(run).dir filesep 'tbv_watch' filesep ds.idstr drinfiles(ii).name(idxcut(1)+1:end)])
        end
        disp('DRIN files saved!')
    end
else
    disp('Not at MR8, not moving DRIN files')
end

% --- Rename TBV settings file -> set format to SHORT
% fname = [ds.dirs.tbvs filesep ds.idstr '-' num2str(run) '.tbv'];
% fprintf('Renaming %s\n',fname)
% filecontents = fileread(fname);
% idxn =  strfind(filecontents,ds.exp.drinfirstfname);
% if isempty(idxn)
%     fprintf('%s not found in %s\n',ds.exp.drinfirstfname,fname);
% else
%     filecontents(idxn:idxn+length(ds.exp.drinfirstfname)-1) = [ds.idstr repmat('0',1,length(ds.exp.drinfirstfname)-length(ds.idstr))];
%     fid = fopen(fname,'w');
%     fprintf(fid,filecontents);
%     fclose(fid);  
%     disp('TBV settings file renamed!')
% end


fprintf('*******************************************\n')
