import scipy.io
import numpy as np
import cv2
import glob
import os

matfiles = glob.glob('results/fgbg_vgg/*.mat')

for file in matfiles:

	mat = scipy.io.loadmat(file)
	base = os.path.basename(file)
	filename = os.path.splitext(base)[0]

	pred = mat['pred']
	pred = (pred + 1) / 2
	for i in range(15):
		tmp = pred[i,:,:,:]
		tmp = np.reshape(tmp, (256,256,3))
		cv2.imwrite('results/imgs/' + filename + '_' + str(i) + '.png', tmp*255)