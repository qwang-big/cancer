import cv2
import math
import sys
import numpy as np
from skimage import color

f = sys.argv[1]
v = '78,11,3,17,11,6'
v = '69,17,-17,21,11,13'
v = [float(x) for x in v.split(',')]

img = cv2.imread(f)

img_shape = img.shape
tile_size = (300, 300)
offset = (300, 300)


def normalize_tile(tile, NormVec):
    Lab = color.rgb2lab(tile)
    TileMean = [0,0,0]
    TileStd = [1,1,1]
    newMean = NormVec[0:3] 
    newStd = NormVec[3:6]
    for i in range(3):
        TileMean[i] = np.mean(Lab[:,:,i])
        TileStd[i] = np.std(Lab[:,:,i])
        # print("mean/std chanel " + str(i) + ": " + str(TileMean[i]) + " / " + str(TileStd[i]))
        tmp = ((Lab[:,:,i] - TileMean[i]) * (newStd[i] / TileStd[i])) + newMean[i]
        if i == 0:
            tmp[tmp<0] = 0 
            tmp[tmp>100] = 100 
            Lab[:,:,i] = tmp
        else:
            tmp[tmp<-128] = 128 
            tmp[tmp>127] = 127 
            Lab[:,:,i] = tmp
    tile = (color.lab2rgb(Lab) * 255).astype(np.uint8)
    return tile


for i in range(int(math.ceil(img_shape[0]/(offset[1] * 1.0)))):
    for j in range(int(math.ceil(img_shape[1]/(offset[0] * 1.0)))):
        cropped_img = img[offset[1]*i:min(offset[1]*i+tile_size[1], img_shape[0]), offset[0]*j:min(offset[0]*j+tile_size[0], img_shape[1])]
        cropped_img = normalize_tile(cropped_img, v)
        cv2.imwrite("debug_" + str(j) + "_" + str(i) + ".jpg", cropped_img)
