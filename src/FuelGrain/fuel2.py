import cv2
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from skimage import io, color

img = cv2.imread('/Users/paniz/Desktop/file2.jpg',0)

thresh = cv2.adaptiveThreshold(img, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 255, 19)
thresh = cv2.bitwise_not(thresh)
plt.subplot(221)
plt.title('Local adapatative Threshold')
plt.imshow(thresh, cmap="gray", vmin=0, vmax=255)

kernel = np.array([[0, 0, 1, 0, 0],
                [1, 1, 1, 1, 1],
                [0, 1, 1, 1, 0],
              [1, 1, 1, 1, 1],
              [0, 0, 1, 0, 0]], dtype=np.uint8)

# Dilatation et erosion
img_dilation = cv2.dilate(thresh, kernel, iterations=1)
img_erode = cv2.erode(img_dilation,kernel, iterations=1)
img_erode = cv2.medianBlur(img_erode, 1)
plt.subplot(221)
plt.title('Dilatation + erosion')
plt.imshow(img_erode, cmap="gray", vmin=0, vmax=255)

ret, labels = cv2.connectedComponents(img_erode)
label_hue = np.uint8(179 * labels / np.max(labels))
blank_ch = 255 * np.ones_like(label_hue)

labeled_img = cv2.merge([label_hue, blank_ch, blank_ch])

labeled_img = cv2.cvtColor(labeled_img, cv2.COLOR_HSV2BGR)
labeled_img[label_hue == 0] = 0

plt.subplot(222)
plt.title('Objects counted:'+ str(ret-1))
plt.imshow(labeled_img)
print('objects number is:', ret-1)
plt.show()

import cv2
import numpy as np
analysis = cv2.connectedComponentsWithStats(img_erode)
(totalLabels, label_ids, values, centroid) = analysis
area_circle = values[0, cv2.CC_STAT_AREA]
inner_area = values[1, cv2.CC_STAT_AREA]
for i in range(100):
# Apply the Component analysis function
    
    area_circle = values[0, cv2.CC_STAT_AREA]
    inner_area = values[1, cv2.CC_STAT_AREA]
    print("circle "+ str(area_circle))
    print(inner_area)
    if(inner_area<area_circle):
        kernel = np.array([[0, 0, 1, 0, 0],
                [1, 1, 1, 1, 1],
                [0, 1, 1, 1, 0],
              [1, 1, 1, 1, 1],
              [0, 0, 1, 0, 0]], dtype=np.uint8)
        dialation = cv2.dilate(img_erode,kernel,iterations = i)
        plt.imshow(dialation,cmap='gray', interpolation='nearest') #collecting the final image
        plt.show()
         #plotting the image
        analysis = cv2.connectedComponentsWithStats(dialation)
        (totalLabels, label_ids, values, centroid) = analysis
        area_circle = values[0, cv2.CC_STAT_AREA]
        inner_area = values[1, cv2.CC_STAT_AREA]
        

