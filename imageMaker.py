import cv2
import numpy as np
import imageanalyzer2
import matplotlib.pyplot as plt



width = 1000
average = 120
freq = 15
T = 5
dt = 0.0005

t = np.arange(0, T, dt)
wave = np.sin(2 * np.pi * freq * t) + np.random.randint(-5,5, len(t)) + average
# plt.plot(t, wave)
# plt.show()

wave_obj = imageanalyzer2.Time(wave, dt)
dst_img = wave_obj.image(width).T

cv2.imshow('dst_img',dst_img)
cv2.waitKey(0)


cv2.imwrite("ref.tif", dst_img)