% Check ds.mat online quickly. needs ds and run
% load ds, define run

fprintf('***check ds run %d********************\n',run)

if run == 1
    disp('ds');     disp(ds);
    disp('ds.exp'); disp(struct2table(ds.exp));
end

disp('ds.r');   disp(struct2table(ds.r));
fprintf('ds.r(%d).rtp\n',run); disp(ds.r(run).rtp);

fprintf('Experiment timings\n')
t = diff(ds.r(run).timings.rtp);
fprintf('RTP arrival M±SD: %.3f±%.3f, range %.3f-%.3f seconds\n',mean(t),std(t),min(t),max(t))
t = diff(ds.r(run).screenflips.timings(contains(ds.r(run).screenflips.type, {'Rest','Regulate'}),1));
fprintf('Screen flips M±SD: %.3f±%.3f, range %.3f-%.3f seconds\n',mean(t),std(t),min(t),max(t))

%% Show activity time course from TBV

figure('pos',[136         242        1169         563]); sp = [3 1];
cmap = colormap(lines);

subplot(sp(1),sp(2),1)
plot(ds.r(run).rtp.raw(:,1),'o-','Color',cmap(1,:));
hold on;
s1 = scatter(ds.exp.vols.regu(:), ds.r(run).rtp.raw(ds.exp.vols.regu(:),1),'filled','MarkerFaceColor',cmap(1,:));
plot(ds.r(run).rtp.raw(:,2),'o-','Color',cmap(3,:));
s2 = scatter(ds.exp.vols.regu(:), ds.r(run).rtp.raw(ds.exp.vols.regu(:),2),'filled','MarkerFaceColor',cmap(3,:));
s3 = plot(ds.exp.vols.rest([6 10],:),repmat(ds.r(run).rtp.allbl(:,1)',2,1),'Color',[1 0 0],'LineWidth',3);
plot(ds.exp.vols.rest([6 10],:),repmat(ds.r(run).rtp.allbl(:,2)',2,1),'Color',get(s3(1),'Color'),'LineWidth',3)
title('Raw BOLD (filled: regulate, empty: rest)')
legend([s1 s2 s3(1)],'Left IPS','Right IPS','Baseline');xlim([0 210]); xlabel('volume');ylabel('BOLD')

subplot(sp(1),sp(2),2)
s1 = plot(ds.exp.vols.regu,reshape(ds.r(run).rtp.pscraw(ds.exp.vols.regu(:),1),15,8),'o-','MarkerFaceColor',cmap(1,:),'Color',cmap(1,:)); hold on
s2 = plot(ds.exp.vols.regu,reshape(ds.r(run).rtp.pscraw(ds.exp.vols.regu(:),2),15,8),'o-','MarkerFaceColor',cmap(3,:),'Color',cmap(3,:));
s3 = plot(ds.exp.vols.regu,reshape(ds.r(run).rtp.dpsc(ds.exp.vols.regu(:)),15,8),'Color',[1 0 0],'LineWidth',2);
plot([0 210],[0 0],'Color',.5*ones(1,3));
title('PSC during regulate blocks'); legend([s1(1) s2(1) s3(1)],'Left IPS','Right IPS','dPSC','Location','southeast')
xlim([0 210]); xlabel('volume');ylabel('PSC')

subplot(sp(1),sp(2),3)
drawlvl = cell2table(vertcat(ds.r(run).rtp.drawlvl{:}),'VariableNames',strsplit('lvl displvl colour')); % disable for smooth transitions
plot(ds.exp.vols.regu, reshape(drawlvl.displvl(ds.exp.vols.regu(:))-7,15,8),'o-','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:)); hold on; 
for l = -6:6,plot([0 210],[l l],'Color',.5*ones(1,3));end
plot([0 210],[0 0],'Color',[1 0 0]);
title('Presented feedback level')
xlim([0 210]); xlabel('volume'); ylabel('Level'); ylim([-6 6])

%% Check motion parameters

fn = dir([ds.r(run).dir filesep 'tbv_target' filesep '*.log']);
fid = fopen([fn.folder filesep fn.name],'r');
if fid ~= -1
    tline = fgetl(fid); tlines={};
    while ischar(tline)
        tlines{end+1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid); clear fid tline;

    tmp = cellfun(@(x) strsplit(x),tlines(cellfun(@(x) contains(x,'deg'),tlines))','UniformOutput',false); tmp = vertcat(tmp{:});
    mp = str2double(tmp(:,[7:3:end]));

    figure
    subplot(2,1,1)
    plot(mp(:,1:3));
    xlabel('volumes');ylabel('mm'); xlim([0 210]);legend(strsplit('x y z'),'location','northwest')
    title(sprintf('Max translation x: %.4f y: %.4f z: %.4f mm',max(abs(mp(:,1:3)))));
    subplot(2,1,2)
    plot(mp(:,4:6));
    xlabel('volumes');ylabel('deg'); xlim([0 210]);legend(strsplit('rx ry rz'),'location','northwest')
    title(sprintf('Max rotation rx: %.4f ry: %.4f rz: %.4f degs',max(abs(mp(:,4:6)))));

    if any(max(abs(mp))>1)==0, disp('All motion parameters within bounds (1mm, 1deg)');end
else
    disp([ds.r(run).dir filesep 'tbv_target' filesep '*.log not found'])
end
%% Show display
pt = ds.r(run).pt;
figure('pos', pt.win.ws/3)
for i = 1:length(ds.r(run).screenflips.imgArray)
imagesc((ds.r(run).screenflips.imgArray{i}));axis('off');
title(i)
pause(0.2)
end
close all;

%% Check draw_NF function
pt = load_PTBsettings(bool_fullscreen,1);
Screen('Preference', 'SkipSyncTests', 1);
[pt.win.w, pt.win.windowRect] = PsychImaging('OpenWindow', mr8*pt.win.screenNumber, pt.stim.colors.dgrey);
pt = load_PTBsettings(bool_fullscreen,2,pt);
Screen('TextSize', pt.win.w, pt.win.screenXpixels/24);
pt.win.WindowSize=Screen('WindowSize', pt.win.w);pt.win.Resolution = Screen('Resolution', pt.win.w);

isrest = 0; istraining = 1; dpsc = 1;
drawlvl = draw_NF(pt,isrest,istraining, dpsc);
Screen('Flip', pt.win.w);

pause(2)
Screen('closeall');

% pt.stim.thwidth = round(tan(deg2rad(1))*viewDis/pixSize);





