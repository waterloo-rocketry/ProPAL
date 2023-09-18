import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

class edge_detection():

    mtx = []

    def __init__(self,inputArray) :
        self.mtx = inputArray
    

    def edge_detector(self):
        fuel_grain = self.mtx.mean(axis =2)
        height,width = fuel_grain.shape
        size_of_square_matrix = 3 #size of submatrix
        try:
            
            vef_v = np.array([[1,0,-1],[1,0,-1],[1,0,-1]])
            vef_h = np.array([[1,1,1],[0,0,0],[-1,-1,-1]])
            #Calculating the total number of possible square submatrices can be formed from the given matrix and the output matrix.
            no_of_square_matrices = ((width-size_of_square_matrix)+1)*((height-size_of_square_matrix)+1)
            sosm=size_of_square_matrix
            row=(height-size_of_square_matrix)+1                     #It is the number of rows it will cover to find out the desired output.
            column=(width-size_of_square_matrix)+1                  #It is the number of columns it will cover to find out the desired output.
            max_val=0
            result_v = []
            result_h = []
            r = (height-3)+1                                         #row value of horizontal edge filter
            c = (width-3)+1                                         #column value of horizontal edge filter

            for i in range(row):
                for j in range(column):
                    sq = fuel_grain[i:i+sosm,j:j+sosm]
                    sum_v = 0
                    sum_h = 0
                    for k in range(3):
                        for l in range(3):
                            sum_v += (sq[k,l] * vef_v[k,l])
                            sum_h += (sq[k,l] * vef_h[k,l])
                    result_v.append(sum_v)
                    result_h.append(sum_h)
            result_matrix_v = np.asarray(result_v).reshape(r,c)     #reshaping the resultant matrix
            result_matrix_h = np.asarray(result_h).reshape(r,c)  
            image = np.sqrt((result_matrix_v**2) + (result_matrix_h**2)) #combining both the vector values.
            plt.imshow(image, cmap='gray', interpolation='nearest') #collecting the final image
            plt.show() #plotting the image
        except Exception as e:
            print("Invalid Input, Try again",str(e))



#img = mpimg.imread('/Users/paniz/Edge-Detection-of-an-ImageFile-using-Numpy/lena.png')  
#m_input = edge_detection(img)
#m_input.edge_detector()






    

    