% function for initialisation eyelink, start and stop recording, and
% transferring files
% Lulu 190918

function eye = NF_eyelink(which,edfFile,pt)
fprintf('***NF_eyelink.m****************************\n')
eye = true;

switch which
    case 'init'
        try
            res = EyelinkInit();
            if res==1, eye = true; disp('Eyelink connected');end
        catch
            eye = false; disp('Eyelink not connected')
        end
    case 'startr'
        try
            i = Eyelink('Openfile', edfFile);
            if i~=0, fprintf('Cannot create EDF file ''%s'' ', edfFile); return; end

            Eyelink('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox for NSL experiment''');
            Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, pt.win.screenXpixels-1, pt.win.screenYpixels-1);
            Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, pt.win.screenXpixels-1, pt.win.screenYpixels-1);
            [~,vs] = Eyelink('GetTrackerVersion');
            Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
            Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
            Eyelink('StartRecording');
            fprintf('Eyelink connected (ver %s)\nEye used: %d\n',vs,Eyelink('EyeAvailable'));
        catch
            disp('Proceed without eye tracking')
        end
    case 'calib'
        d = 0.25;
        gridpos = round([pt.win.xCenter       pt.win.yCenter;                   pt.win.xCenter       pt.win.yCenter*(2-d);
                   pt.win.xCenter*(2-d) pt.win.yCenter;                   pt.win.xCenter       pt.win.yCenter*d;
                   pt.win.xCenter*d     pt.win.yCenter;                   pt.win.xCenter*(2-d) pt.win.yCenter*(2-d);
                   pt.win.xCenter*d     pt.win.yCenter*(2-d);                   pt.win.xCenter*d     pt.win.yCenter*d;
                   pt.win.xCenter*(2-d) pt.win.yCenter*d;                   pt.win.xCenter       pt.win.yCenter]);
        for d = 1:size(gridpos,1)
            msg = sprintf('Calibration point %d: %dx %dy',d,gridpos(d,:)); 
            disp(msg); 
            try
                Eyelink('message', msg); 
            catch
                fprintf('%s\tError: Eyelink msg not logged\n',datestr(now,'yyyy-mm-dd HH:MM:SS')); 
            end
            Screen('DrawDots',pt.win.w,gridpos(d,:),pt.stim.crosswidth,[1 1 1],[],2);
            Screen('Flip', pt.win.w);
            WaitSecs(3);
            try
                Eyelink('message', [msg ' end']); 
            catch
                fprintf('%s\tError: Eyelink msg2 not logged\n',datestr(now,'yyyy-mm-dd HH:MM:SS')); 
            end
        end
    case 'stopr'
        Eyelink('StopRecording'); % Eyelink('Command', 'set_idle_mode'); 
        Eyelink('CloseFile');
        fprintf('Finished recording edf data file ''%s''\n', edfFile );
        
        % Transfer doesn't work -> transferred files are 0KB
%         try
%             fprintf('Receiving edf data file ''%s''\n', edfFile );
%             % status=Eyelink('ReceiveFile');
%             if status > 0, fprintf('ReceiveFile status %d\n', status); end
%             if 2==exist(edfFile, 'file'), fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd); end
%         catch
%             fprintf('Problem receiving data file ''%s''\n', edfFile );
%         end
end
fprintf('*******************************************\n')

%         i = Eyelink('Openfile', edfFile);
%         if i~=0, fprintf('Cannot create EDF file ''%s'' ', edfFile); return; end
%         Eyelink('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox for NSL experiment''');
%         Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, pt.win.screenXpixels-1, pt.win.screenYpixels-1);
%         Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, pt.win.screenXpixels-1, pt.win.screenYpixels-1);
%         [v,vs] = Eyelink('GetTrackerVersion');
%         Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
%         Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
%         ds.explog = update_elog(ds.explog,sprintf('Eyelink connected (ver %s)',vs)); disp(ds.explog{end});
%         Eyelink('StartRecording'); el.eye_used = Eyelink('EyeAvailable'); % Start recording a bit before the start of the scan