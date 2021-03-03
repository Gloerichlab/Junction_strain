# -*- coding: utf-8 -*-
"""
@author: Jooske Monster
first used in publication: Monster et al. 2021 (JCB), https://doi.org/10.1083/jcb.202001042
disclaimer: I am not an experienced programmer so it will definitely not be the smartest script ever
feel free to use and edit script, but please refer to us accordingly :) 


Added additional annotations for clarity 03/03/2021

----------------------------------------------------------------------------------------------------------
STATISTICS OF CELL AND JUNCTION PROPERTIES OF A SPECIFIC CELL AND ITS NEIGHBORS (2X, interphase and upon entry into mitosis)

input: 
     1. STACK of segmented (skeletonized) images of cells (segments) and junctions (outlines) 
     (used SeedWaterSegmenter:https://pypi.org/project/SeedWaterSegmenter/ )
     2. desired cell that you want to extract data from, over 2 frames 
     for this, check the labels (=intensity) in software such as ImageJ
     for instance, follow a cell upon entry into mitosis (cell #14) from frame 2 to 4


output: 
     (individual) images of nodes (tricellular junctions)
     excel sheet of properties of cells, with both the information on all cells of the desired frames and of the individual cell + neighbors
     excel sheet of properties of junctions, with both the information on all junctions of the desired frames and of the individual cell + neighbors
     
1. opens stack and asks for user input to determine which cells and frames to analyze   
2. using skimage.measure.regionprops extract properties of individual cells
3. detecting nodes (tricellular junctions) 
4. define junctions and couple this to nodes, to later be able to couple them back to your desired cell
5. measuring the distance between the nodes with a Breadth-first search (BFS)
6. determine the desired cell and its neighbors, and their corresponding junctions
7. repeat 2-6, but with frame 2
8. write everything to excel
     
"""

#from skimage import io
#from matplotlib import pyplot as plt
import numpy as np
from tifffile import imread, imwrite
import math
from skimage import measure
import pandas as pd
from collections import deque
import os


pixel_size = 0.2063492                       #check every time! # 1 pixel is 0.2 micron
image_c_ori = imread("pos_06_segments.tif")  #adjust your stack name here!
image_j_ori=imread("pos_06_outlines.tif")

#-----------------------------------------------------------------------------------------------------------------
#FUNCTIONS


#function to make a 2x2 bin
def surrounding_pixels(x): 
     pixels = []
     pixels.append((x[0],x[1]))
     pixels.append((x[0]+1,x[1]))
     pixels.append((x[0],x[1]+1))
     pixels.append((x[0]+1,x[1]+1))
     return (pixels)

#searches where three or more objects intercept
def find_nodes(img): 
     img_nodes = np.zeros(img.shape, dtype=np.uint16)                           #empty picture for nodes to fill (16-bit)
     NodeID = []                                                                
     NodeArray = []                                                             
     NodeCells =[]   
     NodeIndex = []                                                          
     nodecounter = 0                                                           
     for index, x in np.ndenumerate(img):                                       #goes through all pixels
          if (index[0] < img.shape[0]-1) and (index[1] < img.shape[1]-1):       #to get rid of borders where error will occur
               a = surrounding_pixels(index)                                    #takes surrounding 2x2 pixels, and retrieves their values (=cell ID)
               pix1=img[a[0]]
               pix2=img[a[1]]
               pix3=img[a[2]]
               pix4=img[a[3]]
               pix1234=[pix1, pix2, pix3, pix4]    
               if len(set(pix1234))>2:                                          #if there are more than 2 unique values, it means that three cells intersect
                    NodeIndex.append(nodecounter)
                    nodecounter +=1
#                    NodeID.append("Node%03d_TCJ%s"%(nodecounter, list(set(pix1234))))                       #makes list of NodeIDs ("Node001")
                    NodeID.append("Node%03d"%nodecounter)
                    NodeArray.append(a)                                         #makes list of locations of nodes(2x2 array xy location)
                    NodeCells.append(list(set(pix1234)))                        #makes list of cells in the node (removes doubles with set)     
                    img_nodes[a[0]]=1                                           #it then fills the four pixels with 1 in the new image
                    img_nodes[a[1]]=1 
                    img_nodes[a[2]]=1 
                    img_nodes[a[3]]=1
#     plt.imshow(img_nodes)
     return img_nodes, NodeID, NodeArray, NodeCells, NodeIndex                             #and it returns an image and the lists for the nodes

def find_junctions(Nodelist):
#junction classifications
     JunctionID = []          # "Junction001_cells1-2
     JunctionIndex = []       #"0"
     JunctionCells = []       # [1,2]     
#per node
     JNode1_Index = []        # "0" 
     JNode2_Index = []        # "1
     JNode1_ID = []           # "Node001" 
     JNode2_ID = []           # "Node002
     JNode1_Cells = []        # [1,2,4]
     JNode2_Cells = []        # [1,2,6]
     JNode1_Array = []        # [(16, 236), (17, 236), (16, 237), (17, 237)]
     JNode2_Array = []        # [(16, 236), (17, 236), (16, 237), (17, 237)]
#     JunctionNodesIndex = [] #[0,1]
#     JunctionNodesID = [] #[Node001, Node002]
#     JunctionNodesArray = [] 
     count=0
     for node1 in Nodelist:
          for y in range(Nodelist.index(node1)+1, len(Nodelist)):
               node2=Nodelist[y]
               if len(set(node1).difference(set(node2))) == 1:
                    index_1 = NodeCells.index(node1)                  #index node 1
                    index_2 = NodeCells.index(node2)                  #index node 2                   
                    JunctionIndex.append(count)
                    count +=1                    
 #                   name="Junction%03d_cells%d/%d"%(count, cells[0], cells[1] ) #name for junction (Junction001_cells1/22)
                    name ="Junction%03d"%(count)
                    JunctionID.append(name)
                    cells=list(set(NodeCells[index_1]) & set(NodeCells[index_2])) #makes list with the two cells where the junction is in between
                    JunctionCells.append(cells)                    
                    JNode1_Index.append(NodeIndex[index_1])
                    JNode2_Index.append(NodeIndex[index_2])                   
                    JNode1_ID.append(NodeID[index_1])
                    JNode2_ID.append(NodeID[index_2])     
                    JNode1_Cells.append(NodeCells[index_1])
                    JNode2_Cells.append(NodeCells[index_2])
                    JNode1_Array.append(NodeArray[index_1])
                    JNode2_Array.append(NodeArray[index_2])                    
     return (JunctionID,      #0
             JunctionIndex,   #1
             JunctionCells,   #2
             JNode1_Index,    #3   
             JNode2_Index,    #4 
             JNode1_ID,       #5  
             JNode2_ID,       #6
             JNode1_Cells,    #7
             JNode2_Cells,    #8
             JNode1_Array,    #9
             JNode2_Array,)   #10         
    
# Function to find the shortest path between a given source and destination.
# disclaimer: I didn't come op with this myself. I adapted it from: https://www.geeksforgeeks.org/shortest-path-in-a-binary-maze/
def BFS(mat, jna1, jna2):   
     
    class Point:                                                                # To store matrix cell cordinates                                             
         def __init__(self,coordinate:tuple): 
              self.x = coordinate[0] 
              self.y = coordinate[1]
      
    class queueNode:                                                            # A data structure for queue
         def __init__(self,pt: Point, dist: int): 
             self.pt = pt                                                       # The cordinates of the cell 
             self.dist = dist                                                   # Cell's distance from the source
         
    rowNum = [-1,-1,-1, 0, 0, 1, 1, 1]#[-1, 0, 0, 1] #[-1,-1,-1, 0, 0, 1, 1, 1] #row and column numbers 
    colNum = [-1, 0, 1,-1, 1,-1, 0, 1]#[0, -1, 1, 0] #[-1, 0, 1,-1, 1,-1, 0, 1] #of 8 neighbours of a given point  
        
    ROW = mat.shape[0]
    COL = mat.shape[1]
    
    
    src = Point(jna1[0])
    dest = Point(jna2[0])

    def isValid(row: int, col: int):                                            # Check whether given point(row,col)
         return (row >= 0) and (row < ROW) and (col >= 0) and (col < COL)       # is a valid point or not 

             
    while mat[src.x][src.y]!=1:                                                 # checks which of the points in srcnode has value 1
        for i in range (4):
             src=  Point(jna1[i])                               
    while mat[dest.x][dest.y]!=1:                                               # checks which of the points in destnode has value 1
        for i in range (4):
             dest=  Point(jna2[i])  
          
    visited = [[False for i in range(COL)] for j in range(ROW)]    
    visited[src.x][src.y] = True                                                # Mark the source cell as visited  
      
    q = deque()                                                                 # Create a queue for BFS 
    s = queueNode(src,0)                                                        # Distance of source cell is 0 
    q.append(s)                                                                 # Enqueue source cell 
         
    while q:                                                                    # Do a BFS starting from source cell
        curr = q.popleft()                                                      # Dequeue the front cell         
        pt = curr.pt 

        if pt.x == dest.x and pt.y == dest.y:                                   # stop when reached dest 
            return (curr.dist*pixel_size)                                                    # then returns distance
               
        for i in range(8):                                                      # Otherwise enqueue its adjacent cells 
            row = pt.x + rowNum[i] 
            col = pt.y + colNum[i] 
              
            if (isValid(row,col) and mat[row][col] == 1 and not visited[row][col]): # if adjacent cell is valid, has path   
                visited[row][col] = True                                        # and not visited yet, enqueue it. 
                Adjcell = queueNode(Point((row,col)),curr.dist+1) 
                q.append(Adjcell) 


          
    return -1                                                                   # Return -1 if destination cannot be reached  

#calculates direct distance between two points
def calculateDistance(src, dest):  
     x1=src[0]
     x2=dest[0]
     y1=src[1]
     y2=dest[1]
     dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)  
     return dist  

#checks which junctions are in, and attached to it
def junctions_incells(cell, ncell1, ncell2): 
     Junction_incell = [] 
     Junction_outcell = []
     for i in range (len(ncell1)):
          if cell in ncell1[i] and cell in ncell2[i]:
               Junction_incell.append(i)     
          elif cell in ncell1[i] or cell in ncell2[i]:                
               Junction_outcell.append(i)
     return (Junction_incell, Junction_outcell)



#-----------------------------------------------------------------------------------------------------------
#1. opening images and creating folders



print ("This script will detect the properties of a cell and its neighbors in 2 desired frames. For instance, before and upon entry into mitosis")
print ("For this, pick a cell (check its label in software like ImageJ) and decide which frames to analyse.")

cell1=int(input("(interphase) What cell are you interested in? "))
frame1=int(input("(interphase) What frame is the mitotic cell in? "))

cell2=int(input("(mitosis) What cell are you interested in? "))
frame2=int(input("(mitosis) What frame is the mitotic cell in? "))

if os.path.isdir("/Nodes")==False:    
     os.mkdir("Nodes") 
     
if os.path.isdir("/Measurements")==False:    
     os.mkdir("Measurements")

outputfile_junctions="Measurements/output_junctions_t%sc%s-t%sc%s.xlsx"%(frame1, cell1, frame2, cell2)
outputfile_cellshape="Measurements/output_cellshape_t%sc%s-t%sc%s.xlsx"%(frame1, cell1, frame2, cell2)


nodeimagename1="Nodes/Nodes_t%s.tiff"%(frame1)
nodeimagename2="Nodes/Nodes_t%s.tiff"%(frame2)


image_c=image_c_ori[frame1-1]          #correcting for that frame 1 = 0 in stack
image_j=image_j_ori[frame1-1]

#-----------------------------------------------------------------------------------------------------------
#2. extracting properties of desired cell shapes # frame 1


#which properties of cells to measure
propList = ["label",
            'area',
            "bbox_area",
            'convex_area', 
            'eccentricity', 
            'equivalent_diameter',
            "orientation",
            'minor_axis_length',
            'perimeter',
            'solidity',] 

#measures properties of propList for all cells
cells = measure.regionprops_table(image_c, properties = propList)
df_c1=pd.DataFrame(cells) 

#converts properties to correct pixel sizes etc.
df_c1=df_c1*[1,
         pixel_size**2,            #Convert pixel square to um square
         pixel_size**2,            #Convert pixel square to um square
         pixel_size**2,            #Convert pixel square to um square
         1,                        #ratio, so no conversion
         pixel_size,               #Convert pixel to um
         57.2958,                  #Convert to degrees from radians
         pixel_size,               #Convert pixel to um
         pixel_size,               #Convert pixel to um
         1]                        #ratio, so no conversion

#-----------------------------------------------------------------------------------------------------------
#3. detecting nodes (tricellular junctions) # frame 1


#creating image with nodes
nodeimage=find_nodes(image_c)[0]
imwrite(nodeimagename1, nodeimage)
NodeID = find_nodes(image_c)[1]
NodeIndex = find_nodes(image_c)[4]
NodeArray = find_nodes(image_c)[2]
NodeCells = find_nodes(image_c)[3]


#--------------------------------------------------------------------------------------------------------------------- 
#4. finds junctions, couples them to cell # and creates lists of these properties               
JunctionID = find_junctions(NodeCells)[0]
JunctionIndex = find_junctions(NodeCells)[1]
JunctionCells = find_junctions(NodeCells)[2]

JNode1_Index=find_junctions(NodeCells)[3]   
JNode2_Index=find_junctions(NodeCells)[4]
JNode1_ID=find_junctions(NodeCells)[5]
JNode2_ID=find_junctions(NodeCells)[6]
JNode1_Cells=find_junctions(NodeCells)[7]
JNode2_Cells=find_junctions(NodeCells)[8]
JNode1_Array=find_junctions(NodeCells)[9]
JNode2_Array=find_junctions(NodeCells)[10]

#-------------------------------------------------------------------------------------------------------------------
# 5. measures the junction length between two nodes with a Breadth-first search (BFS)

ShortestPath = []                                                        #junctions through BFS 
DistanceNodes = []                                                         #junctions through calculate distance
for jun in range(len(JNode1_Array)):                                     #loops through all junction positions to calculate BFS
     length = BFS(image_j, JNode1_Array[jun], JNode2_Array[jun])                                        #measures length with BFS
     ShortestPath.append(length)
     distance = calculateDistance(JNode1_Array[jun][0], JNode2_Array[jun][0])*pixel_size
     DistanceNodes.append(distance)


#makes a pandas dataframe of junctions information    
Junctions = list(zip(JunctionIndex, JunctionID, JunctionCells,
                     JNode1_ID, JNode2_ID, JNode1_Array, JNode2_Array, 
                     JNode1_Cells, JNode2_Cells, ShortestPath, DistanceNodes))

df_j1 = pd.DataFrame(data=Junctions, columns= ["JunctionIndex", "JunctionID", "Cells junction", 
                                            "NodeID1", "NodeID2", "NodeArray1", "NodeArray2",
                                            "NodeCells1", "NodeCells2", "ShortestPath", "DistanceNodes"])
 

#--------------------------------------------------------------------------------------------------------------------- 
# 6. Define the cells and junctions that are indicated, and its neighbors

junction_in = junctions_incells(cell1, JNode1_Cells, JNode2_Cells)[0]
junction_out = junctions_incells(cell1, JNode1_Cells, JNode2_Cells)[1] 

#defines cells that are mitotic cell and next to mitotic cells
cell_out=[]
for i in junction_in:
     cell_out.append(JunctionCells[i])                                     #takes the two cells from each junction
cell_out = list(set([item for sublist in cell_out for item in sublist]))   #and makes a unique values list of this



#selects the rows from all the junctions for the junctions that are in or attached to mitotic cell
df_j2 = df_j1[df_j1["JunctionIndex"].isin(junction_in)]
df_j3 = df_j1[df_j1["JunctionIndex"].isin(junction_out)]
df_c2 = df_c1[df_c1["label"].isin(cell_out)]


#----------------------------------------------------------------------------------------------------
# 7. same script, but for cell in frame #2


image_c=image_c_ori[frame2-1]
image_j=image_j_ori[frame2-1]



#measures properties of propList for all cells
cells = measure.regionprops_table(image_c, properties = propList)
df_c4=pd.DataFrame(cells) 

#converts properties to correct pixel sizes etc.
df_c4=df_c1*[1,
         pixel_size**2,            #Convert pixel square to um square
         pixel_size**2,            #Convert pixel square to um square
         pixel_size**2,            #Convert pixel square to um square
         1,                        #ratio, so no conversion
         pixel_size,               #Convert pixel to um
         57.2958,                  #Convert to degrees from radians
         pixel_size,               #Convert pixel to um
         pixel_size,               #Convert pixel to um
         1]                        #ratio, so no conversion


nodeimage=find_nodes(image_c)[0]
imwrite(nodeimagename2, nodeimage)
NodeID = find_nodes(image_c)[1]
NodeIndex = find_nodes(image_c)[4]
NodeArray = find_nodes(image_c)[2]
NodeCells = find_nodes(image_c)[3]
       
JunctionID = find_junctions(NodeCells)[0]
JunctionIndex = find_junctions(NodeCells)[1]
JunctionCells = find_junctions(NodeCells)[2]

JNode1_Index=find_junctions(NodeCells)[3]   
JNode2_Index=find_junctions(NodeCells)[4]
JNode1_ID=find_junctions(NodeCells)[5]
JNode2_ID=find_junctions(NodeCells)[6]
JNode1_Cells=find_junctions(NodeCells)[7]
JNode2_Cells=find_junctions(NodeCells)[8]
JNode1_Array=find_junctions(NodeCells)[9]
JNode2_Array=find_junctions(NodeCells)[10]

ShortestPath = []                                                          #junctions through BFS 
DistanceNodes = []
for jun in range(len(JNode1_Array)):                                     #loops through all junction positions to calculate BFS
     length = BFS(image_j, JNode1_Array[jun], JNode2_Array[jun])                                        #measures length with BFS
     ShortestPath.append(length)
     distance = calculateDistance(JNode1_Array[jun][0], JNode2_Array[jun][0])*pixel_size
     DistanceNodes.append(distance)
  

#makes a pandas dataframe of junctions information    
Junctions = list(zip(JunctionIndex, JunctionID, JunctionCells,
                     JNode1_ID, JNode2_ID, JNode1_Array, JNode2_Array, 
                     JNode1_Cells, JNode2_Cells, ShortestPath, DistanceNodes))

df_j4 = pd.DataFrame(data=Junctions, columns= ["JunctionIndex", "JunctionID", "Cells junction", 
                                            "NodeID1", "NodeID2", "NodeArray1", "NodeArray2",
                                            "NodeCells1", "NodeCells2", "ShortestPath", "DistanceNodes"])

 
junction_in = junctions_incells(cell2, JNode1_Cells, JNode2_Cells)[0]
junction_out = junctions_incells(cell2, JNode1_Cells, JNode2_Cells)[1] 

#defines cells that are mitotic cell and next to mitotic cells
cell_out=[]
for i in junction_in:
     cell_out.append(JunctionCells[i])                                     #takes the two cells from each junction
cell_out = list(set([item for sublist in cell_out for item in sublist]))   #and makes a unique values list of this



#selects the rows from all the junctions for the junctions that are in or attached to mitotic cell
df_j5 = df_j4[df_j4["JunctionIndex"].isin(junction_in)]
df_j6 = df_j4[df_j4["JunctionIndex"].isin(junction_out)]
df_c5 = df_c4[df_c4["label"].isin(cell_out)]

#----------------------------------------------------------------------------------------------------------------
# 8. writes everything to excel


with pd.ExcelWriter(outputfile_junctions) as writer:
     df_j1.to_excel(writer, sheet_name="i_junctions_all") 
     df_j2.to_excel(writer, sheet_name="i_junctions mitotic cell ") 
     df_j3.to_excel(writer, sheet_name="i_junctions neighboring cell") 
     df_j4.to_excel(writer, sheet_name="m_junctions_all") 
     df_j5.to_excel(writer, sheet_name="m_junctions mitotic cell ") 
     df_j6.to_excel(writer, sheet_name="m_junctions neighboring cell")      

with pd.ExcelWriter(outputfile_cellshape) as writer:     
     df_c1.to_excel(writer, sheet_name="i_cellshape_all")
     df_c2.to_excel(writer, sheet_name="i_cellshape_mitotic+neighbors")
     df_c4.to_excel(writer, sheet_name="m_cellshape_all")
     df_c5.to_excel(writer, sheet_name="m_cellshape_mitotic+neighbors")
