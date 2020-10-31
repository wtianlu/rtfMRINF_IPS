% script with Psychtoolbox settings saved, ready to load

function pt = load_PTBsettings(bool_fullscreen,whichstep,pt)
fprintf('***load_PTBsettings.m %d********************\n',whichstep)
switch whichstep
    case 1
        % -- Keyboard/input information
        KbName('UnifyKeyNames'); FlushEvents('KeyDown');
        pt.input.triggerkey = KbName('5%');
        pt.input.abortkey = KbName('ESCAPE');

        [id,name] = GetKeyboardIndices;
        disp('Available keyboards: ')
        for i=1:length(id), disp([num2str(id(i)) ': ' name{i}]); end
        pt.input.device = id(i); %selects externally connected keyboard
        disp(['Selected keyboard: ' num2str(id(i)) ' ' name{end}])

        %  -- PsychToolbox settings
        PsychDefaultSetup(2); 
        screens = Screen('Screens'); disp('Available screens: ');disp(screens);
        pt.win.screenNumber = max(screens); fprintf('Selected screen: %d\n',pt.win.screenNumber);
        sr = Screen('Resolution',pt.win.screenNumber);
        disp('Native screen settings:');disp(sr);

        if bool_fullscreen==0,ws = [0 0 1920 1080]./2; pt.win.fullscr = false;
        else, ws = [0 0 1920 1080]; pt.win.fullscr = true; end
        pt.win.ws = ws;
        fprintf('Screen size set to: %d %d %d %d\n',ws);

        % -- Define screen colors
        pt.stim.colors.white = WhiteIndex(pt.win.screenNumber);
        pt.stim.colors.black = BlackIndex(pt.win.screenNumber);
        pt.stim.colors.dgrey = pt.stim.colors.white * .5; 
        pt.stim.colors.red = [158 0 0]/255; % red: stop
%         pt.stim.colors.regu = [0 90 0]/255;% green: go
        pt.stim.colors.blue = [0 0 1];
        % color luminance: Evernote WP3 brain/programming 2019.09.24 
    case 2
%----------------------------------------------------------------------
%                       Image information
%----------------------------------------------------------------------
        ws = pt.win.windowRect;
        % -- Open window and get coordinates
        pt.win.screenXpixels = ws(3); pt.win.screenYpixels = ws(4);
        pt.win.xCenter = ws(3)/2; pt.win.yCenter = ws(4)/2;

        % -- Define stimulus sizes
        viewDis = 620; pixSize = (290^2+165^2)^.5/(1920^2+1080^2)^.5; % projector screen in mm, resolution in px
        pt.stim.crosswidth = round(tan(deg2rad(.25))*viewDis/pixSize); % full size
        xCoords = [-pt.stim.crosswidth pt.stim.crosswidth 0 0];
        yCoords = [0 0 -pt.stim.crosswidth pt.stim.crosswidth];
        allCoords = [xCoords; yCoords];
        pt.stim.fixation = {'DrawLines', allCoords,round(pt.stim.crosswidth*.52), ones(1,3), [pt.win.xCenter pt.win.yCenter]};

        pt.stim.linewidth = round(pt.stim.crosswidth*.25); % 3
        pt.stim.thwidth = round(tan(deg2rad(1.2))*viewDis/pixSize); % half size
        pt.stim.thheight = round(tan(deg2rad(6))*viewDis/pixSize); % half size
        pt.stim.thnr = 12; % number of steps excluding 0
        pt.stim.maxlvl = 1; 
        
        disp('pt.stim'); disp(pt.stim)
end
fprintf('*******************************************\n')

