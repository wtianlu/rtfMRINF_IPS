# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
190701 With copied section sfrom TVAEcc_1B 
"""

from psychopy import visual, event, core, logging, gui
import os, datetime, math
from random import shuffle
from iViewXAPI import *
res = iViewXAPI.iV_Connect(c_char_p('127.0.0.1'),c_int(4444),c_char_p('127.0.0.1'),c_int(5555))
if res == 1:
    print 'iView connected'
    @WINFUNCTYPE(None,CSample)
    def sample_callback(sample):
        fsample.write('S {timestamp} {left_gazeX} {left_gazeY} {left_eyeposZ}\n'.format(
                timestamp = sample.timestamp,
                left_gazeX = sample.leftEye.gazeX,
                left_gazeY = sample.leftEye.gazeY,
                left_eyeposZ = sample.leftEye.eyePositionZ))
    @WINFUNCTYPE(None,CEvent)
    def event_callback(event):
        fevent.write('E {fixstart} {eye} {duration} {posX} {posY}\n'.format(
                fixstart = event.startTime,
                eye = event.eye,
                duration = event.duration,
                posX = event.positionX,
                posY = event.positionY))
    eye = True
else:
    print 'Connection failed'
    eye = False

# Read mask images
cwd = os.getcwd()
if not os.path.exists(cwd+"/data"): os.makedirs(cwd+"/data")
        
gd = gui.Dlg()
gd.addField("Subject ID:")
gd.addField('Task:', choices=["Practice", "Experiment"])
gd.addField('Session:', choices=[ "Pre", "T1", "T3", "T2"])
gd.show()
idP = gd.data[0] # raw_input('Enter participant number [xx]: ') # '01'
if int(idP)<10 and len(idP)<2:
    idP = '0'+idP

if gd.data[1]=='Practice':ispractice = True
else:ispractice = False
    
# initialize logging
if ispractice:
    logfname = cwd+'/data/NSL'+idP+gd.data[2]+'_ctvaP.log'
    runs=1
    trialcoderuns = range(24)
else:
    logfname = cwd+'/data/NSL'+idP+gd.data[2]+'_ctvaT.log'
    runs=9
    trialcoderuns = range(24)+range(6)*2

logging.LogFile(f=logfname,level=10,filemode='a')
logging.log(level=logging.DATA, msg='ll_Participant '+idP+' start task '+str(datetime.datetime.now()))
file=open(logfname[:-4]+".txt", "w") #create participant file (rn it's one file per phase)
file.write(str((runs==1)*24+(runs==9)*324) + "\n")
file.close()
print logfname

if eye:
    fsample = open(logfname[:-4]+'_eyeS.txt','w')
    fevent = open(logfname[:-4]+'_eyeE.txt','w') 
    res = iViewXAPI.iV_SetSampleCallback(sample_callback)
    res = iViewXAPI.iV_SetEventCallback(event_callback)
    #iViewXAPI.iV_PauseEyetracking()

# set paradigm specs
mask_imgs = []
for i in os.listdir(cwd+"/masks"):
    if i.startswith ("Mask") and i.endswith(".bmp"):
        mask_imgs.append(cwd+"/masks/"+i)
        
letters = ["A", "B", "D", "E", "F", "G", "H", "J", "K", "L","M", "N", "O", "P", "R", "S", "T", "V", "X","Z"]
maskframes = 30 # 500 ms
expdurs = [1,2,3,5,9,12] # 60 Hz 17 33 50 83 150 200
trialcodeposU = [[0,1],[0,2],[1,2],[3,4],[3,5],[4,5]] # if tr in range(6,12) or range(15,21)
trialcodeposBw = [[[0,3],[0,4],[0,5]],[[1,3],[1,4],[1,5]],[[2,3],[2,4],[2,5]]] # if tr in range(12,15) 
trialcodeposBp = [[[1,3],[1,4],[1,5]],[[2,3],[2,4],[2,5]],[[0,3],[0,4],[0,5]]]# if tr in range(21,24)

gnclock = core.Clock()
logging.setDefaultClock(gnclock)
event.clearEvents()


# Set up window and stimuli
viewDis = 500 # in mm
winRes = [1366, 768] 
#winSize = [(317.5**2/(16**2+9**2))**0.5*16,(317.5**2/(16**2+9**2))**0.5*9]
#pixSize = (winSize[0]**2+winSize[1]**2)**0.5 / (winRes[0]**2+winRes[1]**2)**0.5
pixSize = 317.5/(winRes[0]**2+winRes[1]**2)**0.5
visAng = [2.7, 2.7, 0.5, 7.5]
dimpx = [round(math.tan(math.radians(d))*viewDis/pixSize) for d in visAng]
stimsize,masksize,crosssize,radius = dimpx
posrad = [math.pi/3, 0, -math.pi/3, -math.pi*2/3, math.pi, math.pi*2/3]
pos = [[math.cos(r)*radius,math.sin(r)*radius] for r in posrad]

############################################################################################
# %% -------------------------- FINISHED PREP START EXPERIMENT -------------------------------
win = visual.Window(size=winRes,color=[-1,-1,-1],monitor='testMonitor', units ='pix', screen=0, winType= "pyglet", fullscr=True)
win.setBlendMode(win.blendMode, log=False)
win.mouseVisible = False

fixation = visual.ShapeStim(win, vertices=((0, -0.5*crosssize), (0, 0.5*crosssize), (0,0), (-0.5*crosssize,0), (0.5*crosssize, 0)),
        lineWidth=round(crosssize/6), closeShape=False, lineColor="red",autoLog=False)  
stim_display,mask_display=[],[]
for i in range(len(pos)):
    m = visual.ImageStim(win=win, units = 'pix', size = masksize,autoLog=False)
    s = visual.TextStim(win=win, units = 'pix', text = '', color= 'red', height = stimsize*1.316,font='Arial',bold=False,name="TextStim"+str(i+1))
    mask_display.append(m)
    mask_display[i].setPos((pos[i][0],pos[i][1]))  
    stim_display.append(s)
    stim_display[i].setPos((pos[i][0],pos[i][1]))  
stim_display.append(fixation)
mask_display.append(fixation)
resp_display = visual.TextStim(win=win, units = 'pix', text = '', color= 'white', height = stimsize*1.316,font='Arial',bold=False,autoLog=False)
feedback_display = visual.TextStim(win=win, units = 'pix', text = '', color= 'white', height = stimsize*.3,wrapWidth = winRes[0]*0.8,font='Arial',bold=False)

feedback_display.setText('Welcome to the visual attention task!\nYou will see a number of letters appear briefly on the screen.\n'+\
                         'Please report as many of the RED letters as possible. You will not be timed.\n'+\
                         'Aim to get between 80-90% of your responses correct.\n\nPlease press enter to start. Good luck!')
feedback_display.draw()
win.flip()
win.logOnFlip(level=logging.DATA, msg='wf_waiting for keypress to start the task. Start recording')
iViewXAPI.iV_StartRecording()
wfr = True
while  wfr:
    for key in event.getKeys():
        if key in ['return','space']:
            wfr = False   
        elif key in ['escape']:
            wfr = False  
            win.close()
            core.quit()

for run in range(runs):
#    if eye: iViewXAPI.iV_ContinueEyetracking()
    resp_total,resp_correct = 0,0
    shuffle(trialcoderuns)
    for s in range(3):
        resp_display.setText(str(3-s))
        for i in range(maskframes*2):
            resp_display.draw()
            win.flip()
    
    for tr in range(len(trialcoderuns)):
        ##############################################
        # --------------- TRIAL PREP ---------------
        trc = trialcoderuns[tr]
        savestr = str(trc+1)+"\t{}\t"
        logging.log(level=logging.DATA, msg='ll_trial '+str(tr+1)+' trialcode '+str(trc+1))
        shuffle(letters)
        stiml = letters[0:len(pos)]
        tarstr,disstr = "",""
        for i in range(len(stiml)):
            if (trc in range(6)) or (trc in range(6,12) and i in trialcodeposU[trc-6])\
            or (trc in range(15,21) and i in trialcodeposU[trc-15])\
            or (trc in range(12,15) and i in trialcodeposBw[run%3][trc-12])\
            or (trc in range(21,24) and i in trialcodeposBp[run%3][trc-21]):
                stim_display[i].setText(stiml[i])
                stim_display[i].setColor('red')
                tarstr = tarstr+stiml[i]
                disstr = disstr+"0"
            elif (trc in range(15,21) and not (i in trialcodeposU[trc-15]))\
            or (trc in range(21,24) and not i in trialcodeposBp[run%3][trc-21]):
                stim_display[i].setText(stiml[i])
                stim_display[i].setColor('blue')                
                tarstr = tarstr+"0" 
                disstr = disstr+stiml[i]
            else:
                stim_display[i].setText("")
                tarstr = tarstr+"0" 
                disstr = disstr+"0"
        savestr = savestr+tarstr + "\t" + disstr + "\t"
        
        stim_buffer = visual.BufferImageStim(win,stim=stim_display)
        shuffle(mask_imgs)
        [mask_display[i].setImage(mask_imgs[i]) for i in range(len(pos))]
        mask_buffer = visual.BufferImageStim(win,stim=mask_display)
        
        if trc < 6:
            stimframes = expdurs[trc]
        else:
            stimframes = expdurs[3]
                
        ##############################################  
        # --------------- TRIAL Display ---------------
        win.logOnFlip(level=logging.DATA, msg='wf_fixation run '+str(run+1)+' trial '+str(tr+1))
        for i in range(maskframes*2):
            fixation.draw()
            win.flip()

        if eye: fsample.write('\nTrial {ct} stimuli onset {ts}\n'.format(ct=tr+1,ts = gnclock.getTime()))
        win.logOnFlip(level=logging.DATA, msg='wf_stimuli run '+str(run+1)+' trial '+str(tr+1))
        expdurstart = gnclock.getTime()*1000
        
        for i in range(stimframes):
            stim_buffer.draw()
            win.flip()
        expdur = gnclock.getTime()*1000 - expdurstart

        win.logOnFlip(level=logging.DATA, msg='wf_mask run '+str(run+1)+' trial '+str(tr+1))
        for i in range(maskframes):
            mask_buffer.draw()
            win.flip()

        win.logOnFlip(level=logging.DATA, msg='wf_reporting run '+str(run+1)+' trial '+str(tr+1))
        win.flip()
        
        event.clearEvents()
        responsestr = ""
        wfr = True
        while  wfr:
            for key,kt in event.getKeys(timeStamped=gnclock):
                if key in ['return','space']:
                    if len(responsestr)>0: savestr = savestr + responsestr + "\n"
                    else: savestr = savestr + "-\n"
                    wfr = False
                elif key.upper() in letters and not key.upper() in responsestr:
                    responsestr = responsestr + key.upper()
                    resp_display.setText(responsestr)
                    resp_display.draw()
                    win.flip()
                elif key in ['backspace']:
                    responsestr = responsestr[:-1]
                    resp_display.setText(responsestr)
                    resp_display.draw()
                    win.flip()  
                elif key in ['escape']:
                    win.close()
                    core.quit()
                    wfr = False
#        win.getMovieFrame()
#        win.saveMovieFrames("response"+str(tr+1)+".png")
         
        # Get correct responses
        resp_total = resp_total + len(responsestr)
        resp_correct = resp_correct + sum([l in stiml for l in responsestr])
        file=open(logfname[:-4]+".txt", "a+")
        file.write(savestr.format(round(expdur,1)))
        file.close()
    ##############################################  
    # --------------- BLOCK end -----------------
#    if eye: iViewXAPI.iV_PauseEyetracking()
    win.logOnFlip(level=logging.DATA, msg='wf_feedback run '+str(run+1))
    perc_correct = round(100*resp_correct/resp_total)
    fbstr = "Percentage correct: "+str(perc_correct)+"%\n"
    if perc_correct > 90:
        fbstr = fbstr + "Try to report more letters and feel free to guess more"
    elif perc_correct<80:
        fbstr = fbstr + "Try to only report the letters you are relatively sure of"
    if run < runs-1:  
        fbstr = fbstr + "\n\nPress space to start the next run."
    else:
        fbstr = fbstr + "\n\nThe end."
    feedback_display.setText(fbstr)
    feedback_display.draw()
    win.flip()
    wfr = True
    while  wfr:
        for key,kt in event.getKeys(timeStamped=gnclock):
            if key in ['return','space']:
                wfr = False    
            elif key in ['escape']:
                win.close()
                core.quit()
                wfr = False

logging.log(level=logging.DATA, msg='ll_Participant '+idP+' end task '+str(datetime.datetime.now())) 
logging.flush()
win.close()
core.quit()
iViewXAPI.iV_StopRecording()

