#coding:utf-8

import os
import json
import time
# ffmpeg -i video.m4s -i audio.m4s -c:v copy -c:a aac -strict experimental output.mp4
superPath = os.getcwd()  
partDirs = []  
paths = os.listdir(superPath)  


for p in paths:
    if os.path.isdir(p):
        pp = os.listdir(os.path.join(superPath, p))
        for p2 in pp:
            partDirs.append(os.path.join(superPath, p,p2))

for eatchPath in partDirs[1:]:
    videoJsonDir = eatchPath+'/'+'entry.json'  

   
    videoTitle = ''
    with open(videoJsonDir, 'r', encoding='utf-8') as load_f:
        load_dict = json.load(load_f)
        videoTitle += load_dict['title']
    videoDir = 'video.m4s'
    audioDir = 'audio.m4s'
    videoTitle = "".join(x for x in videoTitle if x.isalnum())
    outDir = '"'+videoTitle+'.mp4"'

 
    # ffmpeg -i video.m4s -i audio.m4s -c:v copy -c:a aac -strict experimental output.mp4
    command = 'cd '+ eatchPath + '/64 && ' #&& 多名命令，根据自己具体情况是32还是64
    command += 'ffmpeg -i ' + videoDir + ' -i ' + audioDir + ' -c copy ' + outDir

    print('{}'+command)
    os.system(command)
