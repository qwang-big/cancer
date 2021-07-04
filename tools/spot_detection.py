#!/usr/bin/env python3


########
#
#   spot_detection.py
#
#   Author: YU Hao  (yuhao@genomics.cn)
#
#   Date:   2019-5-20  (v0.1)
#
########


import sys

import cv2 as cv
import numpy as np
import scipy.stats as sst


class Detector:
    channel_list = []

    def __init__(self, f_img_taking):

        self.__img_mat = f_img_taking

        self.blob_greyscale_model = np.zeros(self.__img_mat.shape, dtype=np.float32)

    def blob_richly_detect(self):

        kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (15, 15))
        channel = cv.morphologyEx(self.__img_mat, cv.MORPH_TOPHAT, kernel, iterations=3)

        channle_list = (channel,)

        mor_kps1 = []
        mor_kps2 = []

        blob_params = cv.SimpleBlobDetector_Params()

        blob_params.thresholdStep = 2
        blob_params.minRepeatability = 2
        blob_params.minDistBetweenBlobs = 1

        blob_params.filterByArea = True
        blob_params.minArea = 1
        blob_params.maxArea = 10

        blob_params.filterByCircularity = False
        blob_params.filterByConvexity = False

        for img in channle_list:
            blob_params.minThreshold = sst.mode(np.floor(np.reshape(img, (img.size,)) / 2) * 2)[0][0]

            blob_params.filterByColor = False

            mor_detector1 = cv.SimpleBlobDetector.create(blob_params)
            mor_kps1.extend(mor_detector1.detect(255 - img))

            blob_params.filterByColor = True
            blob_params.blobColor = 255

            mor_detector2 = cv.SimpleBlobDetector.create(blob_params)
            mor_kps2.extend(mor_detector2.detect(img))

        mor_kps = mor_kps1 + mor_kps2

        mask_layer = np.zeros(self.__img_mat.shape, dtype=np.uint8)

        for key_point in mor_kps:
            r = int(key_point.pt[1])
            c = int(key_point.pt[0])

            mask_layer[r:(r + 2), c:(c + 2)] = 255

        mask_layer = cv.GaussianBlur(mask_layer, (3, 3), 0)

        detector = cv.SimpleBlobDetector.create(blob_params)
        kps = detector.detect(mask_layer)

        diff_list = []

        for key_point in kps:

            r = int(key_point.pt[1])
            c = int(key_point.pt[0])

            diff = np.sum(channel[r:(r + 2), c:(c + 2)]) / 4 - np.sum(channel[(r - 2):(r + 4), (c - 2):(c + 4)]) / 36

            if diff > 0:
                diff_list.append(int(np.around(diff)))

        outImg = cv.cvtColor(self.__img_mat.copy(), cv.COLOR_GRAY2RGB)

        diff_break = 10

        cut_off = int(sst.mode(np.around(np.divide(np.array(diff_list, dtype=np.float32), diff_break)))[0][0]) - \
                  diff_break / 2

        for key_point in kps:

            r = int(key_point.pt[1])
            c = int(key_point.pt[0])

            if np.sum(channel[r:(r + 2), c:(c + 2)]) / 4 - \
                    np.sum(channel[(r - 2):(r + 4), (c - 2):(c + 4)]) / 36 > cut_off:
                self.blob_greyscale_model = cv.circle(outImg, (c, r), 0, (0, 255, 0), 2)


if __name__ == '__main__':

    img_mat = cv.imread(sys.argv[1], cv.IMREAD_GRAYSCALE)

    detected_image_obj = Detector(img_mat)
    detected_image_obj.blob_richly_detect()

    detected_image = detected_image_obj.blob_greyscale_model

    cv.imwrite(sys.argv[1] + '.debug_det.tif', detected_image)
