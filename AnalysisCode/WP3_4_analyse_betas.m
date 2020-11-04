% PROJECT:      WP3 - rt-fMRI NF for self-regulation of interhemispheric IPS activity
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Get NF GLM betas 
% -------------------------------------------------------------------------
% 2020.06.06 Removed betas from ms; calculate GLMs to validate the TBV output - v2
% 2020.06.25 
function dirs = WP3_4_analyse_betas(dirs)
%% Initialisation

fprintf('\n\n***************** Start WP3_4_analyse_betas *****************\n\n')
load('visualisationsettings')

% Get paths to betas
dirbetas = cellfun(@(x) {[x 'beta_0001.nii']},dirs.results.func); % r, s, p

%% Start selection
clear bdatab
outfn = [dirs.data.main 'Group/WP3_bdatab_native_210.txt'];
if exist(outfn,'file')==0
betas = zeros(dirs.n.r,dirs.n.s,dirs.n.p);
for p = 1%:dirs.n.p
    fprintf('\n------------------------------\nParticipant %d\n\n',p)
    for s = 1:dirs.n.s
        for r = 1:dirs.n.r
            fnroi = sprintf('%ssub-0%d/ses-%d/online/r%d/tbv_target/NSL%dS%d_r%d.roi',...
                dirs.raw.main,p,s+1,r,p,s,r);                              % Get ROIs used during online training
            IPS_native = get_tbvroicoords(fnroi);
            
            Vbetas = spm_vol(dirbetas{r,s,p});                             % Get beta values
            fprintf('S%dr%d\t',s,r)
            for h = 1:2                                                    % Extract vals per IPS
                Dbetas = spm_get_data(Vbetas, IPS_native{h}');
                fprintf('ROI%d %d/%d - ',h,sum(~isnan(Dbetas)),(length(Dbetas)))
                
                % Get mean betas
                betas(r,s,p) = nanmean(Dbetas);
                fprintf('%.2f\t',betas(r,s,p))
                
                % Save data
                tmptab = struct2table(struct('pid',p,'ses',s,'runnr',r,...
                    'istestr',ismember(r,[1 5]),'hem',h,...
                    'istargetroi',vis.hgroup(p,1)==h,'betas',betas(r,s,p)));
                if exist('bdatab','var')==0,bdatab = tmptab;else, bdatab = [bdatab;tmptab];end
            end
            fprintf('\n')
        end
        fprintf('\n')
    end
end
    writetable(bdatab,outfn)
else
    bdatab = readtable(outfn);
    disp([outfn ' loaded!'])
end


%% Visualise raw data
close all;clc;
newplabsp = [2 1 3 4 5 6];
hemmarker = strsplit('T C'); font = 'Arial';

for r = 1:2 % 1:test, 2:training
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font)
nr = size(vis.runnr{r},1);
for p = 1%:dirs.n.p
    tmpbdatab = bdatab(and(bdatab.pid==p,bdatab.istestr==not(r-1)),:);
    
    % Visualise
    subplot(3,2,newplabsp(p));hold on;
    plot([.2 17.8],[0 0],'-','Color',ones(1,3)*.75)                        % zero line
    for h = 1:2
        y = tmpbdatab.betas; y = y(tmpbdatab.hem==h); % Observed
        if r==1,x = [1 3 5 7 9 11];else,x = [1:3 5:7 9:11];end
        for s = 1:3,plot(x([1:nr]+(s-1)*nr),y([1:nr]+(s-1)*nr),'k:');end
        text(x,y,hemmarker{vis.hgroup(p,h)},'FontWeight','bold','Color',...
            vis.cmap(setdiff(1:2,h),:),'HorizontalAlignment',...
            'center','VerticalAlignment','middle')
    end
    xticks([2 6 10]);xticklabels(1:3);xlabel('Session');ylabel('\beta')
    xlim([0 12]);%ylim([-1.8 2.5]);%title(['P' num2str(newplabsp(p))])
    if p == 1, title('Right-IPS group');elseif p == 2, title('Left-IPS group');end 
    text(6,max(ylim),['P' num2str(newplabsp(p))],'VerticalAlignment','top','HorizontalAlignment','center')
    
    if r == 1
        % paired t-test
        x1 = bdatab.betas(and(bdatab.pid==p,bdatab.hem==1),:);
        x2 = bdatab.betas(and(bdatab.pid==p,bdatab.hem==2),:);
        [H,P,CI,STATS] = ttest(x1,x2);
        fprintf('P%d mean L: %.2f, R:%.2f\tH = %d, p = %.4f, t(%d) = %.2f\n',...
            p,mean(x1),mean(x2),H,P,STATS.df,STATS.tstat)
    end
end
end

fprintf('\n\n***************** End WP3_4_analyse_betas *****************\n\n')







