from PIL import Image
import math
import os, sys
import numpy as np
from skimage import color

f = sys.argv[1]
#v = '57,22,-8,20,10,5'
#v = [float(x) for x in v.split(',')]

def get_labstat(img):
    v = np.zeros(6)
    Lab = color.rgb2lab(img)
    for i in range(3):
        v[i] = np.mean(Lab[:,:,i])
        v[i+3] = np.std(Lab[:,:,i])
    return v

img = Image.open(f)
v = get_labstat(img)
img = np.asarray(img)

img_shape = img.shape
tile_size = (150, 150)
final_size = (300, 300)
offset = tile_size


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
        cropped_img = Image.fromarray(cropped_img).resize(final_size, Image.ANTIALIAS)
        cropped_img = np.asarray(cropped_img)
        cropped_img = normalize_tile(cropped_img, v)
        Image.fromarray(cropped_img).save(os.path.splitext(f)[0] + '/' + str(j) + "_" + str(i) + ".jpg")
