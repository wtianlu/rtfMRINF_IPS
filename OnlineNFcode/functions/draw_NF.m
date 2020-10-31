% present feedback displays - circle or ther version with discrete steps

function drawlvl = draw_NF(pt,isrest,istraining, dpsc)

if nargin < 4,     dpsc = 0;
elseif nargin < 3,    istraining = false;
elseif nargin < 2,    isrest = true;
end

twidth = pt.stim.thwidth; %2*pt.stim.thwidth; 
thmin = pt.win.yCenter+pt.stim.thheight; 
thmax = pt.win.yCenter-pt.stim.thheight;
radlist = thmax:(thmin-thmax)/pt.stim.thnr:thmin; radlist = radlist(end:-1:1); % higher index~higher level
psc2lvl = -pt.stim.maxlvl:2*pt.stim.maxlvl/(length(radlist)-1):pt.stim.maxlvl;

[~,lvl]=min(abs(psc2lvl-dpsc));
if isrest == 1
    drawlvl = {ceil(length(radlist)/2),pt.stim.colors.blue}; 
elseif isrest == 0
    drawlvl = {ceil(length(radlist)/2),pt.stim.colors.red};
    if istraining, drawlvl = {lvl,pt.stim.colors.red}; end
end

% --- %% Draw thermometer %% --- %

% level ladder
for i = 1:length(radlist)
    if i == drawlvl{1}
        Screen('DrawLine',pt.win.w,drawlvl{2},...
            pt.win.xCenter-twidth,radlist(i),pt.win.xCenter+twidth,radlist(i),pt.stim.linewidth*2)
    else
        lcol = ones(1,3)*pt.stim.colors.dgrey*4/3;
        Screen('DrawLine',pt.win.w,lcol, pt.win.xCenter-twidth,radlist(i),pt.win.xCenter+twidth,radlist(i),pt.stim.linewidth)
    end
end
        
% --- frame outline
d = 2;
Screen('FrameRect',pt.win.w,ones(1,3),...
    [pt.win.xCenter-twidth thmax-pt.stim.linewidth*d pt.win.xCenter+twidth thmin+pt.stim.linewidth*d],pt.stim.linewidth)

% --- white fixation
Screen(pt.stim.fixation{1},pt.win.w,pt.stim.fixation{2:end}); % Screen('Flip', pt.win.w);

drawlvl = [lvl drawlvl];
fprintf(' - level: %d/%d\n\n',drawlvl{2},length(radlist));


