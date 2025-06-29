import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import make_interp_spline, griddata, LinearNDInterpolator
from PIL import Image
import time
import random

def putIntoDF(file_addy,secnum):
    data = []
    section = 0
    with open(file_addy,'r') as file:
        for line in file:
            line = line.strip()

            if(not line):
                continue
            
            if(line[0] == '#' and line[1] != '#'):
                section += 1

            if(section == secnum):
                if(line[0] == '#' and "projection:" in line):
                    index = line.find("projection:")
                    cols = line[index+len("projection: "):].split()
                    continue
                data.append(list(map(float,line.split())))
            if(section>secnum):
                break
    df = pd.DataFrame(data,columns = cols)
    return(df)
    
s = time.time()

df = putIntoDF("normalOrder.txt",1)
'''
df = putIntoDF("normalOrder.txt",2)
df = putIntoDF("normalOrder.txt",3)
df = putIntoDF("normalOrder.txt",4)
df = putIntoDF("normalOrder.txt",5)
df = putIntoDF("normalOrder.txt",6)
df = putIntoDF("normalOrder.txt",7)
df = putIntoDF("normalOrder.txt",8)
df = putIntoDF("normalOrder.txt",9)
'''
e = time.time()
#cols = df.columns.values.tolist()
#xdata = df[cols[0]]
#ydata = df[cols[1]]
#zdata = df[cols[2]]
#print(zdata)



print(e-s)


'''
def makePlots(in_addy,out_addy,sec_list):
    os.makedirs(out_addy, exist_ok=True)
    
    for i in sec_list:
        df = putIntoDF(in_addy,i)
        x_col = df.columns[0]
        y_col = df.columns[1]
        df.plot(kind='scatter',x=x_col,y=y_col,title=x_col+'vs'+y_col)
        fileName = "{}.png".format(i)
        fileName = os.path.join(directory,fileName)
        plt.savefig(fileName)
        plt.close()

def side_by_side(addy1, addy2):
    img1 = Image.open(addy1).convert("RGB")
    img2 = Image.open(addy2).convert("RGB")
    maxH = max(img1.height, img2.height)
    img1 = img1.resize((int(img1.width * maxH / img1.height), maxH))
    img2 = img2.resize((int(img2.width * maxH / img2.height), maxH))

    # Combine the widths and max height for the new image
    combW = img1.width + img2.width
    combImg = Image.new("RGB", (combW, maxH))  # Correct size order: width first, height second

    combImg.paste(img1, (0, 0))  # Paste img1 at the left side
    combImg.paste(img2, (img1.width, 0))  # Paste img2 to the right of img1

    return combImg

#MADE EXPLICITLY FOR DATA TEMPLATE P1,P2,P3,....PN-1,CHISQ
#Go in order 3D,2D,1D section numbers


def ListCreate(file_addy, secnum):
    dfnD = putIntoDF(file_addy, secnum)
    cols = dfnD.columns.values.tolist()
    #Min = (dfnD[cols[-1]].min())

    if len(cols) == 2:
        x = dfnD[cols[0]].to_numpy()
        y = dfnD[cols[1]].to_numpy()
        func=make_interp_spline(x, y, k=1)
    else:
        L = [dfnD[col].to_numpy() for col in cols[:-1]]
        points = np.column_stack(L)
        values = dfnD[cols[-1]].to_numpy()
        func=LinearNDInterpolator(points, values)

    return func

def ChisQ(t12,t13,t23,dcp,dms,dma,funcList,minList):
    chiSq = (funcList[0](t23,dma,dcp) - minList[0])
    + (funcList[1](t12,dms) - minList[1])
    + (funcList[2](t13) - minList[2])
    return chiSq


s = time.time()

func3 = ListCreate("normalOrder.txt",1)
end3d = time.time()
print("3D", end3d - s)
func2 = ListCreate("normalOrder.txt",4)
func1 = ListCreate("normalOrder.txt",17)
print("2D,1D",time.time()-end3d)
print(func3(0.256,1.6,-130), func2(0.19,-3.27), func1(0.01))

e = time.time()

print(e-s)
'''
'''
df1 = putIntoDF("normalOrder.txt",17)
df = putIntoDF("normalOrder.txt",4)

cols = df.columns.values.tolist()
cols1 = df.columns.values.tolist()
x = tuple(df[cols[0]])
y = tuple(df[cols[1]])
x1 = tuple(df[cols1[0]])
y1 = tuple(df[cols1[1]])
#print(df)
#print(cols)
#print(x)
points = np.array(tuple(zip(x,y)))
values = np.array(df[cols[2]])

#z_interp = griddata(points, values, [(0.17, -6)], method='linear')
z_func = LinearNDInterpolator(points,values)
#z_func2 = LinearNDInterpolator(points,values)
oneD = make_interp_spline(x1,y1,k=1)
zlist = [z_func,oneD]
print(zlist[0]+zlist[1])

'''


'''

# Get the list of files from both directories
files1 = sorted(os.listdir("InvOrder"))
files2 = sorted(os.listdir("NormOrder"))

# Construct the full paths for both sets of images
images1 = [os.path.join("InvOrder", f) for f in files1]
images2 = [os.path.join("NormOrder", f) for f in files2]

combined_files = []

# Combine images side by side
for img1, img2 in zip(images1, images2):
    combI = side_by_side(img1, img2)
    combined_files.append(combI)

# Save the combined images as a PDF
first_image = combined_files[0]
rest_images = combined_files[1:]

# Save the combined images to a PDF, using the first image and appending the rest
first_image.save("output.pdf", save_all=True, append_images=rest_images)



'''




'''
def funcCreate(file_addy):
    df1D = putIntoDF(file_addy,17)
    df2D = putIntoDF(file_addy,4)
    df3D = putIntoDF(file_addy,1)

    cols1D = df1D.columns.values.tolist()
    cols2D = df2D.columns.values.tolist()
    cols3D = df3D.columns.values.tolist()

    x1D =tuple(df1D[cols1D[0]])
'''
