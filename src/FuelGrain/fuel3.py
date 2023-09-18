import cv2
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from skimage import io, color
import math


img = cv2.imread('/Users/paniz/Desktop/file.jpg',0)
thresh = cv2.adaptiveThreshold(img, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 253 ,0)
thresh = cv2.bitwise_not(thresh)

#plt.subplot(221)
#plt.title('Local adapatative Threshold')
plt.imshow(thresh, cmap="gray", vmin=0, vmax=255)

kernel = np.array([[0, 0, 1, 0, 0],
                [1, 1, 1, 1, 1],
                [0, 1, 1, 1, 0],
              [1, 1, 1, 1, 1],
              [0, 0, 1, 0, 0]], dtype=np.uint8)

# Dilatation et erosion
img_dilation = cv2.dilate(thresh, kernel, iterations=1)
img_erode = cv2.erode(img_dilation,kernel, iterations=1)
# clean all noise after dilatation and erosion
img_erode = cv2.medianBlur(img_erode, 1)
plt.subplot(221)
plt.title('Dilatation + erosion')
plt.imshow(img_erode, cmap="gray", vmin=0, vmax=255)
kernel2 =  np.ones((5,5),np.uint8)
opening = cv.morphologyEx(img_erode, cv.MORPH_OPEN, kernel2)
#closing = cv.morphologyEx(img_erode, cv.MORPH_CLOSE, kernel2)

ret, labels = cv2.connectedComponents(opening)




inner_obj = np.array(labels, dtype=np.uint8)
inner_obj[labels == 4] = 255
plt.imshow(inner_obj)
plt.show()



import cv2
import numpy as np
analysis = cv2.connectedComponentsWithStats(opening)
(totalLabels, label_ids, values, centroid) = analysis
bg = values[1, cv2.CC_STAT_AREA]
outter_cir = values[2, cv2.CC_STAT_HEIGHT]
inner_cir =  values[3, cv2.CC_STAT_HEIGHT]

obj =  values[4, cv2.CC_STAT_HEIGHT]

cir_area = inner_cir**2 * math.pi
print(cir_area)
i =1
print(obj)
print(inner_cir)
for i in range(50):

    if(obj < cir_area):
        kernel = np.array([[0, 0, 1, 0, 0],
                [1, 1, 1, 1, 1],
                [0, 1, 1, 1, 0],
              [1, 1, 1, 1, 1],
              [0, 0, 1, 0, 0]], dtype=np.uint8)
            
        erosion = cv2.dilate(inner_obj,kernel,iterations = i)
        
        plt.subplot(224)
        
        analysis = cv2.connectedComponentsWithStats(erosion)
        (totalLabels, label_ids, values, centroid) = analysis
        
        obj = values[1, cv2.CC_STAT_AREA]
        edged = cv2.Canny(erosion, 50, 200)
        contours,hierarchy= cv2.findContours(edged, 1, 2)
        cnt = contours[0]
        peri = cv2.arcLength(cnt, True)
        blank_image = np.zeros((img.shape[0], img.shape[1], 3))
        cv2.drawContours(blank_image, contours, -1, (0,255,0), 3)
        print(obj)
        print(peri)
        plt.imshow(erosion,cmap='gray') #collecting the final image
        plt.imshow(blank_image)
        plt.show()
        
         #plotting the image

        
        
        

