import scipy.io
import numpy as np
import cv2

mat = scipy.io.loadmat('results/fgbg_vgg/IMG_0302.mat')
x = mat['X']
x = (x+1)/2
y = mat['Y']
y = (y+1)/2

tmpx = x[0,:,:,:]
tmpx = np.reshape(tmpx, (256,256,3))
cv2.imwrite('0302.png', tmpx*255)
tmpy = y[8,:,:,:]
tmpy = np.reshape(tmpy, (256,256,3))
cv2.imwrite('gd8.png', tmpy*255)

mat = scipy.io.loadmat('results/fgbg_vgg/IMG_0309.mat')
x = mat['X']
x = (x+1)/2
tmpx = x[0,:,:,:]
tmpx = np.reshape(tmpx, (256,256,3))
cv2.imwrite('0309.png', tmpx*255)
