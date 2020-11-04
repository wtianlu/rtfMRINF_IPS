% Lulu 2020.06.05
% Extract coordinates from *.roi used during online data processing
% case when there are 2 ROIs

function coords = get_tbvroicoords(fnroi)
% Extract ROI coordinates in resolution 96x96x52
fid=fopen(fnroi,'r'); tline = fgetl(fid); tlines = cell(0,1);
while ischar(tline), tlines{end+1,1} = tline; tline = fgetl(fid); end
fclose(fid);
idxstart = find(contains(tlines,'NrOfVoxels'));
nvox = cellfun(@(x) {strsplit(x)},tlines(contains(tlines,'NrOfVoxels')));
nvox = cellfun(@(x) str2double(x(2)),nvox);
for h = 1:2
coords{h} = cellfun(@(x) {str2double(strsplit(x))},tlines(idxstart(h)+[1:nvox(h)]));
coords{h} = vertcat(coords{h}{:}); coords{h} = coords{h}(:,2:end);
end
[~,fn,ext] = fileparts(fnroi);
fprintf('Extracting ROI coordinates from %s%s\tNumber of voxels: %d - %d\n',...
    fn,ext,nvox(1),nvox(2))