# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a python script for creating a data structure to study orders in which chromosome arm aneuploidies (CAAs) are acquired during tumour evolution.
"""
#from IPython import get_ipython
def reset():
    try:
        get_ipython().magic('reset -sf')  #analysis:ignore
    except NameError:
        pass
reset()
#import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [16,9]
import networkx as nx
import xlsxwriter

#import xlrd



global Arms, ChromLen, BuiltLinks, NumberOfPaths, BestChromArray, TempPaths

BestChromArray = [['Path Value', 'Starting Chromosome', 'Path']]
  


#Node class definitions
class HeadNode:
    def __init__(self, initdata):
        self.level = initdata[0]    # Level n means that the poiter 'chromlist' is the head of the chromosome list with n CAAs
        self.chromlist = initdata[1]
        self.next = None

    def getLevel(self):
        return self.level
        
    def getChromList(self):
        return self.chromlist
        
    def getNext(self):
        return self.next

    def setLevel(self, newlevel):
        self.level = newlevel
    
    def setChromList(self, newchromlist):
        self.chromlist = newchromlist

    def setNext(self,newnext):
        self.next = newnext


class ChromNode:
    def __init__(self,initdata):
        self.chrom = initdata
        self.mark = [None]*ChromLen # for each arm i, eventually mark[i] will be either None, True or False. If mark[i] = True, it means that the parent has no abnormal ith arm. 
        self.count = 1
        self.parent = False
        self.children = []
        self.transitions = []
        self.general_transitions = np.zeros(2*ChromLen, dtype=int)
        self.final_transitions_sum = np.zeros(2*ChromLen, dtype=int)
        self.final_transitions_sum4 = np.zeros(2*ChromLen, dtype=int)
        self.final_transitions_prod = np.zeros(2*ChromLen, dtype=int)
        #self.final_transitions_sum0 = np.zeros(2*ChromLen, dtype=int)
        self.next  = None

    def getChrom(self):
        return self.chrom
        
    def getCount(self):
        return self.count

    def getChildren(self):
        return self.children
    
    def getParent(self):
        return self.parent
    
    def getNext(self):
        return self.next

    def setChrom(self, newdata):
        self.chrom = newdata
        
    def setCount(self, newcount):
        self.count = newcount

    def incCount(self):
        self.count = self.count + 1        
        
    def setNext(self,newnext):
        self.next = newnext
        
    def addChild(self, newchild):
        self.children.append(newchild)
        
        
def printPathLengths(pt, n):
    global NumberOfPaths
    children = pt.getChildren()
    if children == []:
        NumberOfPaths[n] = NumberOfPaths[n] + 1
#        print("Path length:", n)
    else:
        for ptr in children:
            printPathLengths(ptr, n+1)
            
def secondTransition(pt, j, b):
    pathcount = 0
    datasetcount = 0
    if pt.chrom[j] == b:
        pathcount = pathcount + 1
        datasetcount = datasetcount + pt.count
    elif pt.chrom[j] == 0:    
        for child in pt.children:
            x  = secondTransition(child, j, b)
            pathcount = pathcount + x[0]
            datasetcount = datasetcount + x[1]

    return([pathcount, datasetcount])
            
            
def firstTransition(pt, i, a, j, b):
    pathcount = 0
    datasetcount = 0
    if pt.chrom[j] == 0:
        if pt.chrom[i] == a and pt.mark[i] == True:
            for child in pt.children:
                x  = secondTransition(child, j, b)
                pathcount = pathcount + x[0]
                datasetcount = datasetcount + x[1]
                
        elif pt.chrom[i] == 0: 
            for child in pt.children:
                x = firstTransition(child, i, a, j, b)
                pathcount = pathcount + x[0]
                datasetcount = datasetcount + x[1]
    return([pathcount, datasetcount])
    
def addPaths(path, j, b):
    global TempPaths
    for chpt in path[-1].children:
        if chpt.chrom[j] == 0:
            path.append(chpt)
            addPaths(path,j, b)
        elif chpt.chrom[j] == b and path[-1].chrom[j] == 0:
            path.append(chpt)
            TempPaths.append(path)
            

def pathValue_to_b(pt, j, b):
    global TempPaths
    TempPaths = []
    if pt.chrom[j] == 0:
        addPaths([pt], j, b)
    else:
        print("Error: in the function pathValue_to_b")
        exit()
        
    nodes = []
    for row in TempPaths:
        nodes = nodes + row
    nodes = set(nodes)
    
    pathvalue = 0
    for ppt in nodes:
        pathvalue = pathvalue + ppt.count
#    if TempPaths != []:
#        print('-------------------------------------\n',TempPaths)
    return pathvalue
        

def addBestPaths(pt, value, headChrom, pathVec):
    global BestChromArray
    if pt.children == []:
        if len(np.nonzero(headChrom)[0]) == 0:
            headChromConv = 'Normal'
        else:
            headChromConv = ''
            indSet = np.nonzero(headChrom)[0]
            for i in indSet:
                if headChrom[i] == +1:
                    headChromConv = headChromConv + ' +' + Arms[i]
                elif headChrom[i] == -1:
                    headChromConv = headChromConv + ' -' + Arms[i]

        BestChromArray.append([value, headChromConv, pathVec])
        return
    else:
        for chpt in pt.children:
            indSet = np.nonzero(pt.chrom - chpt.chrom)[0]
            if len(indSet) > 1:
                print('Error: Child should not have more than one extra abnormal norm compared to the parent')
                exit()
            i = indSet[0]
            
            if chpt.chrom[i] == 1:
                addBestPaths(chpt, value + chpt.count, headChrom, pathVec + ' -> +' + Arms[i] )
            elif chpt.chrom[i] == -1:
                addBestPaths(chpt, value + chpt.count, headChrom, pathVec + ' -> -' + Arms[i] )

            
            
            
            
#definition of List class      
class ChromList:
    def __init__(self):
        self.head = None
        
        
    #returns True if the list is empty            
    def isEmpty(self):
        return self.head == None
    
    #Adds item to the front of the list
    def nodeAdd(self,item):
        temp = ChromNode(item)
        temp.setNext(self.head)
        self.head = temp
    
        
    #Returns the size as a list (x,y), where x is the total number of times all the y CAAs appeared
    def size(self):
        current = self.head
        count_c = 0
        count_t = 0
        while current != None:
            count_t = count_t + 1
            count_c = count_c + current.getCount()
            current = current.getNext()
    
        return (count_c, count_t)
    
    
    # Returns vector with i^th elemnt indicating the number of times i^th CAA appeared in the data
    def chromSizes(self):
        current = self.head
        repData = []
        while current != None:
            repData.append(current.getCount())
            current = current.getNext()
    
        return repData
    
    
    #Searches item in the list
    def search(self,item):
        current = self.head
        
        while current != None:
            if np.array_equal(current.getChrom(), item):
                return current
            else:
                current = current.getNext()
    
        return None
    
    def printList(self):
        global Arms, ChromLen
        current = self.head
        count = 1
        chromCount = []
        nApp = []
        AbNormVec = []
        while current != None:
            chromCount.append(count)
#            print('------- Chromosome: ', count, ' ---------------')
            nApp.append(current.getCount())
            item = current.getChrom()
            AbArms = ' '
            for i in range(ChromLen):
                if item[i] == 1:
                    AbArms = AbArms+'+' + Arms[i]+' '
                elif item[i] == -1:
                    AbArms = AbArms+'-' + Arms[i]+' '
            AbNormVec.append(AbArms)    
            current = current.getNext()
            count = count + 1
        data = {'Chromosome': chromCount, '# appearances': nApp, 'Abnormal arms': AbNormVec}
        dff = pd.DataFrame(data = data)
        dff = dff[['Chromosome', 'Abnormal arms','# appearances']]
        print(dff.to_string(index=False))

    def printChildren(self):
        global ChromLen
        current = self.head
        count = 1

        while current != None:
            item = current.getChrom()
            temparray =[]
            for i in range(ChromLen):
                if item[i] == 1:
                    temparray.append('+' + Arms[i])
                elif item[i] == -1:
                    temparray.append('-' + Arms[i])
            print('\n Chromosome ', count,': ', ' '.join(temparray))
            print('\n Parent Status: ', current.parent)
            chcount = 1
            for pt in current.getChildren():
                item = pt.getChrom()
                
                temparray =[]
                for i in range(ChromLen):
                    if item[i] == 1:
                        temparray.append('+'+Arms[i])
                    elif item[i] == -1:
                        temparray.append('-' + Arms[i])
                print('\t Child ', chcount,': ', ' '.join(temparray), '(count :', pt.count,')')
                chcount = chcount + 1
                    
            current = current.getNext()
            count = count + 1
            
    def countChildren(self):

        current = self.head
        count = 1
        chromCount = []
        AbNormVec = []
        nChild = []
        tcount = 0
        while current != None:
            item = current.getChrom()
            AbArms =' '
            for i in range(ChromLen):
                if item[i] == 1:
                    AbArms = AbArms + '+' + Arms[i]+' '
                elif item[i] == -1:
                    AbArms = AbArms + '-' + Arms[i] + ' '
            chromCount.append(count)
            AbNormVec.append(AbArms)
            chcount = len(current.getChildren())
            tcount = tcount + chcount
            nChild.append(chcount)
            current = current.getNext()
            count = count + 1  
        data = {'Chromosome': chromCount, 'Abnormal arms': AbNormVec, '# children': nChild}
        dff = pd.DataFrame(data = data)
        dff = dff[['Chromosome', 'Abnormal arms','# children']]
        print(dff.to_string(index=False))        
        return tcount
        
#    #Removes the item from the list 
#    def remove(self,item):
#        current = self.head
#        previous = None
#        found = False
#        while not found:
#            if current.getData() == item:
#                found = True
#            else:
#                previous = current
#                current = current.getNext()
#    
#        if previous == None:
#            self.head = current.getNext()
#        else:
#            previous.setNext(current.getNext())
#    
    
class HeadList:
    def __init__(self):
        self.head = None
        
        
    #returns True if the list is empty            
    def isEmpty(self):
        return self.head == None
    
    # maximum number of CAAs
    def maxLevel(self):
        current = self.head
        level = 0
        while current != None:
            level = current.getLevel()
            current = current.getNext()
        return level
    
    #Adds item to the list making sure that the order is preserved
    def insertChrom(self, item):
        current = self.head
        previous = None
        stop = False
        abarms = np.count_nonzero(item)
        while current != None and not stop:
            if current.getLevel() >= abarms:
                stop = True
            else:
                previous = current
                current = current.getNext()
                
        if current != None and current.getLevel() == abarms:
                pt = current.getChromList().search(item)
                if pt != None:
                    pt.incCount()
                else:
                    current.getChromList().nodeAdd(item)
            
        else: 
            temp_l = ChromList()
            temp_l.nodeAdd(item)
            temp_n = HeadNode((abarms, temp_l))
            if previous == None:
                temp_n.setNext(self.head)
                self.head = temp_n
            else:
                temp_n.setNext(current)
                previous.setNext(temp_n)    
                
    # Summary of the given data set    
    def summary(self):
        current = self.head
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'lp = buildChromList()\'.')
            return
        nChrom = []
        nTotChrom =[]
        nAbNorm = []
        while current != None:
            nChrom.append(current.getChromList().size()[1])
            nTotChrom.append(current.getChromList().size()[0])
            nAbNorm.append(current.getLevel())
#            print('There are ', current.getChromList().size()[1],' groups of ', current.getChromList().size()[0],' chromosomes with ', current.getLevel(),' abnormal arms')
            current = current.getNext()
        data = {'# chromosomes': nChrom, '# appearances in the data':nTotChrom, '# abnormal arms': nAbNorm}
        dff = pd.DataFrame(data=data)
        dff = dff[['# abnormal arms', '# chromosomes', '# appearances in the data']]
        print(dff.to_string(index=False))
        
        
    # List of all the chromosomes with n CAAs
    def abnormalLevel(self, n):
        
        current = self.head
        stop = False
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'lp = buildChromList()\'.')
            return
            
        while current != None and not stop:
            abarms = current.getLevel()
            if  abarms >= n:
                stop = True
            else:
                current = current.getNext()
                
        if current == None or abarms > n:
            print('There is no chromosome with ', n, ' abnormal arms. Use \'<list pointer>.summary()\' to see the summary of the given chromosome dataset')
            return
        else:
            current.getChromList().printList()
            return current.getChromList()
        
    # prints the array of arms
    def arms(self):
        global Arms
#        print('--------------------------------------------------------------')
        print('Sequence of Arms: ')
#        print('--------------------------------------------------------------')
        return Arms
        
    # Builds a tree such that n^th level nodes have n CAAs
    def buildLinks(self):
        global BuiltLinks
        current = self.head
        if current == None:
            print('Error: List is empty. Build the chromosome list using \'<list pointer> = buildChromList()\'.')
            return
        
        nextlist = current.getNext()
        while nextlist != None:
            pt1 = current.getChromList().head
            while pt1 != None:
                chrom1 = pt1.getChrom()

                pt2 = nextlist.getChromList().head
                while pt2 != None:
                    chrom2 = pt2.getChrom()

                    if np.count_nonzero(chrom1 - chrom2) == 1:
                        pt2.parent = True
                        pt1.addChild(pt2)
                    pt2 = pt2.getNext()
                pt1 = pt1.getNext()
            current = nextlist
            nextlist = nextlist.getNext()
        BuiltLinks = True
        
    # Identifying the CAAs that appeared for the first time
    def markFirstTimeNodes(self):
        global BuiltLinks
        current = self.head
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'<list pointer> = buildChromList()\'.')
            return
        if BuiltLinks == False:
            print('Error: Links are not built. Use <list pointer>.buildLinks().')
            return
        while current != None:
            pt = current.getChromList().head
            while pt != None:
                chrom = pt.chrom
                nonzero_ind = np.nonzero(chrom)[0]
                mark = pt.mark
                for i in nonzero_ind:
                    if mark[i] == None:
                        mark[i] = True
                    for child in pt.children:
                        child.mark[i] = False
                
                pt = pt.next
            current = current.next

    # Children of CAA level n
    def childrenOfLevel(self, n):  
        global BuiltLinks
        if not BuiltLinks:
            print('Error: There are no links. Use \'<list pointer>.bulidLinks()\' to create links.')
            return          
        current = self.head
        stop = False
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'lp = buildChromList()\'.')
            return
        
        while current != None and not stop:
            abarms = current.getLevel()
            if  abarms >= n:
                stop = True
            else:
                current = current.getNext()
                
        if current == None or abarms > n:
            print('There is no chromosome with ', n, ' abnormal arms. Use \'<list pointer>.summary()\' to see the summary of the given chromosome dataset')
            return
        else:
            current.getChromList().printChildren()
            
    # Returns array with i^th element indicating the number of times CAA i of level n appeared in the given data
    def chromRepOfLevel(self, n):
        current = self.head
        stop = False
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'lp = buildChromList()\'.')
            return
        
        while current != None and not stop:
            abarms = current.getLevel()
            if  abarms >= n:
                stop = True
            else:
                current = current.getNext()
        if current == None or abarms > n:
            return 0
        else:
            return current.getChromList().chromSizes()
        
        
    # Children of CAA level n
    def countChildrenOfLevel(self, n):  
        global BuiltLinks
        if not BuiltLinks:
            print('Error: There are no links. Use \'<list pointer>.bulidLinks()\' to create links.')
            return           
        current = self.head
        stop = False
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'lp = buildChromList()\'.')
            return
            
        while current != None and not stop:
            abarms = current.getLevel()
            if  abarms >= n:
                stop = True
            else:
                current = current.getNext()
                
        if current == None or abarms > n:
            print('There is no chromosome with ', n, ' abnormal arms. Use \'<list pointer>.summary()\' to see the summary of the given chromosome dataset')
            return
        else:
            tcount = current.getChromList().countChildren()
            print('***********************************************************************')
            print('The total number of links between level ', n,' and level ', n+1,' is ', tcount)
            print('***********************************************************************')

    
    def graph(self):
        global BuiltLinks
        
        current = self.head
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'lp = buildChromList()\'.')
            return
        
        if not BuiltLinks:
            print('Error: There are no links. Use \'<list pointer>.bulidLinks()\' to create links.')
            return 
            
        G = nx.DiGraph() # CAA nodes

        temp = current
        longerChromList = 0
        while temp != None:
            longerChromList = max(longerChromList, temp.getChromList().size()[1])
            temp = temp.getNext()
        max_len = longerChromList*50
        curLevel = current.getLevel()
        curLen = current.getChromList().size()[1]
        curGap = max_len/(curLen + 1)
        pt1 = current.getChromList().head
        i = 1
        while pt1 != None:

            G.add_node(str(i)+':('+str(curLevel)+','+str(pt1.getCount())+')', pos=(curLevel,i*curGap))
            pt1 = pt1.getNext()
            i = i+1
            
            
        nextlist = current.getNext()
        while nextlist != None:

            
            nextLevel = nextlist.getLevel()
            nextLen = nextlist.getChromList().size()[1]
            nextGap = max_len/(nextLen + 1)
            
            pt2 = nextlist.getChromList().head 
            j = 1
            while pt2 != None:
                G.add_node(str(j)+':('+str(nextLevel)+','+str(pt2.getCount())+')', pos=(nextLevel, j*nextGap))
                pt2 = pt2.getNext()
                j = j + 1

                
            pt1 = current.getChromList().head
            i = 1
            while pt1 != None:   
                pt2 = nextlist.getChromList().head  
                j = 1
                while pt2 != None:
                    if pt2 in pt1.getChildren():
                        G.add_edge(str(i)+':('+str(curLevel)+','+str(pt1.getCount())+')', str(j)+':('+str(nextLevel)+','+str(pt2.getCount())+')')
                    pt2 = pt2.getNext()
                    j = j + 1
                    
                pt1 = pt1.getNext()
                i = i+1
            
            current = nextlist
            curLevel = nextLevel
            nextlist = nextlist.getNext()
        
        pos=nx.get_node_attributes(G,'pos')
        nx.draw(G, pos,font_size=8, with_labels=True, node_size=1000, node_color='#b3ecff')
        plt.show()
        return G
    

#### transitionGraph for method 4 
    def transitionGraphMethod4(self):
        global BuiltLinks, Arms, ChromLen
        
        current = self.head
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'lp = buildChromList()\'.')
            return
        
        if not BuiltLinks:
            print('Error: There are no links. Use \'<list pointer>.bulidLinks()\' to create links.')
            return 


            
        G = nx.DiGraph() # CAA nodes

        extArms1 = []
        extArms2 = []
        
        for x in Arms:
            extArms1.append(list(('+'+x, 0)))
            extArms1.append(list(('-'+x, 0)))
            extArms2.append('+'+x+' ')
            extArms2.append('-'+x+' ')
        
        maxCount = 0
        for i in range(ChromLen):            
            current = self.head
            
            while current.next != None:
                
                p_chrom = current.getChromList().head

                while p_chrom != None:
                    
                    for pt in p_chrom.children:
                        if pt.chrom[i] == +1:
                            extArms1[i][1] = extArms1[i][1] + pt.count
                        elif pt.chrom[i] == -1:
                            extArms1[i+1][1] = extArms1[i + 1][1] + pt.count
                    p_chrom = p_chrom.next
                current = current.next
                
                
            maxCount = max(maxCount, extArms1[i][1])
            
            
        gap = 50
        
        node_colors = []
        for i in range(len(extArms2)):

            G.add_node(extArms1[i][0], pos=(i*gap, 100 ))
            node_colors.append('#FF0000') ### Red color
                               
            G.add_node(extArms2[i], pos=(i*gap, 0))
            node_colors.append('#00CED1') ### Red color
     
            

        T = np.zeros((2*ChromLen, 2*ChromLen)) 


        
        for i in range(ChromLen):

            
            for j in range(ChromLen):            
                current = self.head
                
                while current.next != None:
                    
                    pt = current.getChromList().head

                    while pt != None:
                        row_sum4 = pt.count + sum(pt.final_transitions_sum4)
                        if row_sum4 != 0:
                            vec4 = np.divide(pt.final_transitions_sum4, row_sum4)
                            
                            if pt.chrom[i] == +1:    
                                     
                                T[2*i][2*j] = T[2*i][2*j] + vec4[2*j] 
                                T[2*i][2*j+1] = T[2*i][2*j+1] + vec4[2*j + 1]
                                
                            if pt.chrom[i] == -1:
    
                                T[2*i +1][2*j] = T[2*i +1][2*j] + vec4[2*j]
                                T[2*i +1][2*j+1] = T[2*i +1][2*j+1] + vec4[2*j + 1]                                
                                    

                        pt = pt.next

                    
                    current = current.next 

        for i in range(2*ChromLen):

            row_sum = sum(T[i])
            if row_sum > 0:
                T[i] = np.divide(T[i], row_sum)

            
            for j in range(2*ChromLen):

                if T[i][j] > 0:
                    G.add_edge(extArms1[i][0], extArms2[j], arrowstyle='-|>', color='#696969',  weight=2*T[i][j])
        nodeSizes = np.zeros(2*len(extArms2))
        for i in range(len(extArms2)):
            nodeSizes[2*i] = 20 + 5*extArms1[i][1]*gap/(2*maxCount)
            nodeSizes[2*i +1] = 20 + 5*extArms1[i][1]*gap/(2*maxCount)
        pos=nx.get_node_attributes(G,'pos')
        edges = G.edges()
        edge_colors=[G[i][j]['color'] for i, j in edges]
        weights=[G[i][j]['weight'] for i, j in edges]
#        arrowstyle=[G[i][j]['arrowstyle'] for i, j in edges]
        plt.title('BRCA3 - Method 4 \n Transitions are from red nodes to blue nodes')  
        nx.draw(G, pos, node_size=nodeSizes, font_size=6, with_labels=False, edges=edges, node_color=node_colors, edge_color=edge_colors, width=weights, arrows=False)                


        for p in pos:  # change text positions
            x = list(pos[p])
            if x[1] == 0:
                x[1] -= 2
            else:
                x[1] += 2 
            pos[p] = tuple(x)
        nx.draw_networkx_labels(G, pos, font_size=4)
        plt.show()
        return
    
###### transitionGraph for method 5
    def transitionGraphMethod5(self):       
        global BuiltLinks, Arms, ChromLen
        
        current = self.head
        if current == None:
            print('Error: List is empty. Bulid the chromosome list using \'lp = buildChromList()\'.')
            return
        
        if not BuiltLinks:
            print('Error: There are no links. Use \'<list pointer>.bulidLinks()\' to create links.')
            return 


            
        G = nx.DiGraph() # Chromosome nodes

        extArms1 = []
        extArms2 = []
        
        for x in Arms:
            extArms1.append(list(('+'+x, 0)))
            extArms1.append(list(('-'+x, 0)))
            extArms2.append('+'+x+' ')
            extArms2.append('-'+x+' ')
        
        maxCount = 0
        for i in range(ChromLen):            
            current = self.head
            
            while current.next != None:
                
                p_chrom = current.getChromList().head

                while p_chrom != None:
                    
                    for pt in p_chrom.children:
                        if pt.chrom[i] == +1:
                            extArms1[i][1] = extArms1[i][1] + pt.count
                        elif pt.chrom[i] == -1:
                            extArms1[i+1][1] = extArms1[i + 1][1] + pt.count
                    p_chrom = p_chrom.next
                current = current.next
                
                
            maxCount = max(maxCount, extArms1[i][1])
            
         
        gap = 50
        
        node_colors = []
        for i in range(len(extArms2)):

            G.add_node(extArms1[i][0], pos=(i*gap, 100 ))
            node_colors.append('#FF0000') ### Red color
                               
            G.add_node(extArms2[i], pos=(i*gap, 0))
            node_colors.append('#00CED1') ### Red color
        
            

        T = np.zeros((2*ChromLen, 2*ChromLen)) 


        
        for i in range(ChromLen):

            
            for j in range(ChromLen):            
                current = self.head
                
                while current.next != None:
                    
                    pt = current.getChromList().head

                    while pt != None:

                        if pt.chrom[i] == +1:    
                                 
                            T[2*i][2*j] = T[2*i][2*j] + pt.final_transitions_sum4[2*j]
                            T[2*i][2*j+1] = T[2*i][2*j+1] + pt.final_transitions_sum4[2*j + 1]
                            
                        if pt.chrom[i] == -1:

                            T[2*i +1][2*j] = T[2*i +1][2*j] + pt.final_transitions_sum4[2*j]
                            T[2*i +1][2*j+1] = T[2*i +1][2*j+1] + pt.final_transitions_sum4[2*j + 1]                                
                                

                        pt = pt.next

                    
                    current = current.next 


        for i in range(2*ChromLen):
            row_sum = sum(T[i])
            if row_sum > 0:
                T[i] = np.divide(T[i], row_sum)
            
            for j in range(2*ChromLen):

                if T[i][j] > 0:
                    G.add_edge(extArms1[i][0], extArms2[j], arrowstyle='-|>', color='#696969',  weight=2*T[i][j])
        nodeSizes = np.zeros(2*len(extArms2))
        for i in range(len(extArms2)):
            nodeSizes[2*i] = 20 + 5*extArms1[i][1]*gap/(2*maxCount)
            nodeSizes[2*i +1] = 20 + 5*extArms1[i][1]*gap/(2*maxCount)
            
        pos=nx.get_node_attributes(G,'pos')
        edges = G.edges()
        edge_colors=[G[i][j]['color'] for i, j in edges]
        weights=[G[i][j]['weight'] for i, j in edges]
#        arrowstyle=[G[i][j]['arrowstyle'] for i, j in edges]
        plt.title('BRCA3 - Method 5 \n Transitions are from red nodes to blue nodes')                
        nx.draw(G, pos, node_size=nodeSizes, font_size=6, with_labels=False, edges=edges, node_color=node_colors, edge_color=edge_colors, width=weights, arrows=False)                
        for p in pos:  # change text positions
            x = list(pos[p])
            if x[1] == 0:
                x[1] -= 2
            else:
                x[1] += 2 
            pos[p] = tuple(x)
        nx.draw_networkx_labels(G, pos, font_size=4)
        plt.show()
        return
    
    # Computes number of paths that starts at 0 and terminates at a CAA with n CAAs
    def pathLengths(self):
        global NumberOfPaths
        current = self.head
        root = current.getChromList().head
        NumberOfPaths = np.zeros(self.maxLevel())
        printPathLengths(root, 0)
        print(NumberOfPaths)
        return NumberOfPaths
    
    # Construct a matrix with (i,j) element being the number of paths that have i -> j transitions
    def transitionMatrix(self):
        global Arms, ChromLen
        extArms = []
        
        for x in Arms:
            extArms.append('+'+x)
            extArms.append('-'+x)
            
        TC = [[]]
        TC[0].append(' ')
        
        TP = [[]]
        TP[0].append(' ')
        for x in extArms:
            TC[0].append(x)
            TC.append([x])
            TP[0].append(x)
            TP.append([x])


        
        for i in range(ChromLen):
            TC[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            TC[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            TP[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            TP[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            for j in range(ChromLen):            
                current = self.head
                
                while current.next != None:
                    
                    p_chrom = current.getChromList().head

                    while p_chrom != None:
                        
                        for pt in p_chrom.children:
                            if p_chrom.chrom[i] == 0 and pt.chrom[i] != 0 and pt.chrom[j] == 0:
                                row_sign = pt.chrom[i]
                                
                                col_sign = [0, 0]
                                for ppt in pt.children:
                                    if ppt.chrom[j] == +1:
                                        col_sign[0] = +1
                                    elif ppt.chrom[j] == -1:
                                        col_sign[1] = -1

                                if row_sign == +1:
                                    if col_sign[0] == +1:
                                        TC[2*i+1][2*j+1] = TC[2*i+1][2*j+1] + 1
                                    if col_sign[0] == -1:
                                        TC[2*i+1][2*j+2] = TC[2*i+1][2*j+2] + 1
                                else:
                                    if col_sign[0] == +1:
                                        TC[2*i+2][2*j+1] = TC[2*i+2][2*j+1] + 1
                                    if col_sign[0] == -1:
                                        TC[2*i+2][2*j+2] = TC[2*i+2][2*j+2] + 1

                                col_count = [0, 0]
                                for ppt in pt.children:
                                    if ppt.chrom[j] == +1:
                                        col_count[0] = col_count[0] + ppt.count
                                    elif ppt.chrom[j] == -1:
                                        col_sign[1] = col_count[1] + ppt.count

                                if row_sign == +1:
                                    TP[2*i+1][2*j+1] = TP[2*i+1][2*j+1] + col_count[0]
                                    TP[2*i+1][2*j+2] = TP[2*i+1][2*j+2] + col_count[1]
                                else:
                                    TP[2*i+2][2*j+1] = TP[2*i+2][2*j+1] + col_count[0]
                                    TP[2*i+2][2*j+2] = TP[2*i+2][2*j+2] + col_count[1]

                        p_chrom = p_chrom.next

                    
                    current = current.next 


        workbook = xlsxwriter.Workbook('TransitionsAlongPaths.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, TC[row][col])
        workbook.close()
    
        workbook = xlsxwriter.Workbook('TransitionsInDataSet.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, TP[row][col])
        workbook.close()


    # Construct a matrix with (i,j) element being number of paths along which i appears before j
    def generalTransitionMatrix(self):
        global Arms, ChromLen
        extArms = []
        
        for x in Arms:
            extArms.append('+'+x)
            extArms.append('-'+x)
            
        TC = [[]]
        TC[0].append(' ')
        
        TP = [[]]
        TP[0].append(' ')
        for x in extArms:
            TC[0].append(x)
            TC.append([x])
            TP[0].append(x)
            TP.append([x])


        
        for i in range(ChromLen):
            TC[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            TC[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            TP[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            TP[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            for j in range(ChromLen):            
                current = self.head
                
                while current.next != None:
                    
                    pt = current.getChromList().head

                    while pt != None:
                        
                        if pt.parent == False:
                            x1 = firstTransition(pt, i, +1, j, +1)
                            x2 = firstTransition(pt, i, +1, j, -1)
                            x3 = firstTransition(pt, i, -1, j, +1)
                            x4 = firstTransition(pt, i, -1, j, -1)  
                            
                            TC[2*i +1][2*j+1] = TC[2*i +1][2*j+1] + x1[0]
                            TP[2*i +1][2*j+1] = TP[2*i +1][2*j+1] + x1[1]
                            
                            TC[2*i +1][2*j+2] = TC[2*i +1][2*j+2] + x2[0]
                            TP[2*i +1][2*j+2] = TP[2*i +1][2*j+2] + x2[1]
                            
                            TC[2*i +2][2*j+1] = TC[2*i +2][2*j+1] + x3[0]
                            TP[2*i +2][2*j+1] = TP[2*i +2][2*j+1] + x3[1]
                            
                            TC[2*i +2][2*j+2] = TC[2*i +2][2*j+2] + x4[0]
                            TP[2*i +2][2*j+2] = TP[2*i +2][2*j+2] + x4[1]

                        pt = pt.next

                    
                    current = current.next 

        workbook = xlsxwriter.Workbook('GeneralTransitionsAlongPaths.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, TC[row][col])
        workbook.close()   
        
        
        workbook = xlsxwriter.Workbook('GeneralTransitionsInDataSet.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, TP[row][col])
        workbook.close() 
        
        
    # Compute transition values
    def computeTransitions(self):
#        global Arms, ChromLen
        current = self.head
        while current != None:
            ppt = current.chromlist.head
            while ppt != None:
                if ppt.children != []:
                    parent_count = ppt.count
                    ch_total_count = 0
                    for chpt in ppt.children:
                        ch_total_count = ch_total_count + chpt.count + parent_count
                    for chpt in ppt.children:
                        ppt.transitions.append((chpt, parent_count + chpt.count, (parent_count + chpt.count)/ch_total_count ))
                ppt = ppt.next
            current = current.next
        return
    
                


        
    # Construct a matrix with (i,j) element being transition probability of i -> j 
    def transitionProbabilities(self):
        global Arms, ChromLen
        extArms = []
        
        for x in Arms:
            extArms.append('+'+x)
            extArms.append('-'+x)
            
        T1 = [[]]
        T1[0].append(' ')
        
        T2 = [[]]
        T2[0].append(' ')
        for x in extArms:
            T1[0].append(x)
            T1.append([x])
            T2[0].append(x)
            T2.append([x])


        
        for i in range(ChromLen):
            T1[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T1[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            T2[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T2[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            for j in range(ChromLen):            
                current = self.head

                
                while current.next != None:
                    
                    p_chrom = current.getChromList().head

                    while p_chrom != None:
                        
                        for pt in p_chrom.children:
                            if p_chrom.chrom[i] == 0 and pt.chrom[i] != 0 and pt.chrom[j] == 0:
                                row_sign = pt.chrom[i]
                                

                                for chtrans in pt.transitions:
                                    if chtrans[0].chrom[j] != 0:
                                        if row_sign == +1:
                                            if chtrans[0].chrom[j] == +1:
                                                T1[2*i+1][2*j+1] = T1[2*i+1][2*j+1] + chtrans[1]
                                                T2[2*i+1][2*j+1] = T2[2*i+1][2*j+1] + chtrans[2]                                                
                                            if chtrans[0].chrom[j] == -1:
                                                T1[2*i+1][2*j+2] = T1[2*i+1][2*j+2] + chtrans[1]
                                                T2[2*i+1][2*j+2] = T2[2*i+1][2*j+2] + chtrans[2]                                                
                                        else:
                                            if chtrans[0].chrom[j] == +1:
                                                T1[2*i+2][2*j+1] = T1[2*i+2][2*j+1] + chtrans[1]
                                                T2[2*i+2][2*j+1] = T2[2*i+2][2*j+1] + chtrans[2]
                                            if chtrans[0].chrom[j] == -1:
                                                T1[2*i+2][2*j+2] = T1[2*i+2][2*j+2] + chtrans[1]
                                                T2[2*i+2][2*j+2] = T2[2*i+2][2*j+2] + chtrans[2]
                                                
                        p_chrom = p_chrom.next

                    
                    current = current.next 
        for i in range(len(extArms)):
            row_sum = sum(T1[i+1][1:])
#            print('Row index = ', T1[i+1][0], 'with sum = ', row_sum)
            if row_sum > 0:
                T1[i+1][1:] = np.divide(T1[i+1][1:], row_sum) 
                
            row_sum = sum(T2[i+1][1:])
            if row_sum > 0:
                T2[i+1][1:] = np.divide(T2[i+1][1:], row_sum) 
            

        workbook = xlsxwriter.Workbook('Transition_prob_1.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T1[row][col])
        workbook.close()
    
        workbook = xlsxwriter.Workbook('Transition_prob_2.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T2[row][col])
        workbook.close()
    # fill the general_transitions vactor at each node
    
    def computeGneralTransitions(self):
        global ChromLen

        current = self.head
        while current != None:
            
            pt = current.chromlist.head
            while pt != None:
                for j in range(ChromLen):
                    if pt.chrom[j] == 0:
                        pt.general_transitions[2*j] = pathValue_to_b(pt, j, +1)
                        pt.general_transitions[2*j + 1] = pathValue_to_b(pt, j, -1)
                        
                pt = pt.next
                
            current = current.next

                        
        
    # Construct a matrix with (i,j) element being probability of i appearing before j
    def generalTransitionProbabilities(self):
        global Arms, ChromLen
        extArms = []
        
        for x in Arms:
            extArms.append('+'+x)
            extArms.append('-'+x)
            
        T1 = [[]]
        T1[0].append(' ')
        
        T2 = [[]]
        T2[0].append(' ')
        for x in extArms:
            T1[0].append(x)
            T1.append([x])
            T2[0].append(x)
            T2.append([x])


        
        for i in range(ChromLen):
            T1[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T1[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            T2[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T2[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            for j in range(ChromLen):            
                current = self.head
                
                while current.next != None:
                    
                    pt = current.getChromList().head

                    while pt != None:
                        row_sum = sum(pt.general_transitions)
#                        print('Row sum :-( = ', row_sum)
                        if row_sum != 0:
                            vec = np.divide(pt.general_transitions, row_sum)
                            if pt.mark[i] == True:
                                if pt.chrom[i] == +1:    
                                    
                                    T1[2*i +1][2*j+1] = T1[2*i +1][2*j+1] + pt.general_transitions[2*j]
                                    T2[2*i +1][2*j+1] = T2[2*i +1][2*j+1] + vec[2*j]
                                    
                                    T1[2*i +1][2*j+2] = T1[2*i +1][2*j+2] + pt.general_transitions[2*j+1]
                                    T2[2*i +1][2*j+2] = T2[2*i +1][2*j+2] + vec[2*j + 1]
                                    
                                if pt.chrom[i] == -1:
                                    
                                    T1[2*i +2][2*j+1] = T1[2*i +2][2*j+1] + pt.general_transitions[2*j]
                                    T2[2*i +2][2*j+1] = T2[2*i +2][2*j+1] + vec[2*j]
                                    
                                    T1[2*i +2][2*j+2] = T1[2*i +2][2*j+2] + pt.general_transitions[2*j + 1]
                                    T2[2*i +2][2*j+2] = T2[2*i +2][2*j+2] + vec[2*j + 1]

                        pt = pt.next

                    
                    current = current.next 

                 
        for i in range(len(extArms)):
            row_sum = sum(T1[i+1][1:])
#            print('------------------------------------\n Row = ', T1[i+1], 'with sum = ', row_sum)
            if row_sum > 0:
                T1[i+1][1:] = np.divide(T1[i+1][1:], row_sum) 
                
            row_sum = sum(T2[i+1][1:])
            if row_sum > 0:
                T2[i+1][1:] = np.divide(T2[i+1][1:], row_sum) 
                
        workbook = xlsxwriter.Workbook('GeneralTransition_prob_1.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T1[row][col])
        workbook.close()   
        
        
        workbook = xlsxwriter.Workbook('GeneralTransition_prob_2.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T2[row][col])
        workbook.close() 
    
    
    
    # This is the sub routine for the final transition probabilities
    def computeFinalTransitions(self):
        global ChromLen

        current = self.head
        while current != None:
            
            pt = current.chromlist.head
            while pt != None:
                for chpt in pt.children:
                    for j in range(ChromLen):
                        if chpt.chrom[j] == +1 and pt.chrom[j] == 0:
                            pt.final_transitions_sum[2*j] = pt.count + chpt.count
                            pt.final_transitions_prod[2*j] = pt.count*chpt.count
                            pt.final_transitions_sum4[2*j] = chpt.count
                        elif chpt.chrom[j] == -1 and pt.chrom[j] == 0:
                            pt.final_transitions_sum[2*j + 1] = pt.count + chpt.count
                            pt.final_transitions_prod[2*j + 1] = pt.count*chpt.count
                            pt.final_transitions_sum4[2*j + 1] = chpt.count
                            
                
                pt = pt.next
            current = current.next
    # Construct a matrix with (i,j) element being probability of i appearing before j
    def finalTransitionProbabilities(self):
        global Arms, ChromLen
        extArms = []
        
        for x in Arms:
            extArms.append('+'+x)
            extArms.append('-'+x)
            
        T1 = [[]]
        T1[0].append(' ')
        
        T2 = [[]]
        T2[0].append(' ')
        
        T3 = [[]]
        T3[0].append(' ')

        T4 = [[]]
        T4[0].append(' ')        
        
        T5 = [[]]
        T5[0].append(' ')  
        
        for x in extArms:
            T1[0].append(x)
            T1.append([x])
            T2[0].append(x)
            T2.append([x])
            T3[0].append(x)
            T3.append([x])
            T4[0].append(x)
            T4.append([x])
            T5[0].append(x)
            T5.append([x])

        
        for i in range(ChromLen):
            T1[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T1[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            T2[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T2[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            T3[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T3[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            T4[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T4[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            T5[2*i+1].extend(np.zeros(2*ChromLen, dtype=int))
            T5[2*i+2].extend(np.zeros(2*ChromLen, dtype=int))
            
            for j in range(ChromLen):            
                current = self.head
                
                while current.next != None:
                    
                    pt = current.getChromList().head

                    while pt != None:
                        row_sum = sum(pt.final_transitions_sum)
                        row_sum4 = pt.count + sum(pt.final_transitions_sum4)
                        row_sum5 = sum(pt.final_transitions_sum4)
#                        print('Row sum :-( = ', row_sum)
                        if row_sum5 != 0:
                            vec = np.divide(pt.final_transitions_sum, row_sum)
                            vec4 = np.divide(pt.final_transitions_sum4, row_sum4)

                            if pt.chrom[i] == +1:    
                                
                                T1[2*i +1][2*j+1] = T1[2*i +1][2*j+1] + pt.final_transitions_sum[2*j]
                                T2[2*i +1][2*j+1] = T2[2*i +1][2*j+1] + vec[2*j]
                                T3[2*i +1][2*j+1] = T3[2*i +1][2*j+1] + pt.final_transitions_prod[2*j]
                                T4[2*i +1][2*j+1] = T4[2*i +1][2*j+1] + vec4[2*j]     
                                T5[2*i +1][2*j+1] = T5[2*i +1][2*j+1] + pt.final_transitions_sum4[2*j]
                                
                                T1[2*i +1][2*j+2] = T1[2*i +1][2*j+2] + pt.final_transitions_sum[2*j+1]
                                T2[2*i +1][2*j+2] = T2[2*i +1][2*j+2] + vec[2*j + 1]
                                T3[2*i +1][2*j+2] = T3[2*i +1][2*j+2] + pt.final_transitions_prod[2*j+1]
                                T4[2*i +1][2*j+2] = T4[2*i +1][2*j+2] + vec4[2*j + 1]
                                T5[2*i +1][2*j+2] = T5[2*i +1][2*j+2] + pt.final_transitions_sum4[2*j + 1]
                                
                            if pt.chrom[i] == -1:
                                
                                T1[2*i +2][2*j+1] = T1[2*i +2][2*j+1] + pt.final_transitions_sum[2*j]
                                T2[2*i +2][2*j+1] = T2[2*i +2][2*j+1] + vec[2*j]
                                T3[2*i +2][2*j+1] = T3[2*i +2][2*j+1] + pt.final_transitions_prod[2*j]
                                T4[2*i +2][2*j+1] = T4[2*i +2][2*j+1] + vec4[2*j]
                                T5[2*i +2][2*j+1] = T5[2*i +2][2*j+1] + pt.final_transitions_sum4[2*j]
                                
                                T1[2*i +2][2*j+2] = T1[2*i +2][2*j+2] + pt.final_transitions_sum[2*j + 1]
                                T2[2*i +2][2*j+2] = T2[2*i +2][2*j+2] + vec[2*j + 1]
                                T3[2*i +2][2*j+2] = T3[2*i +2][2*j+2] + pt.final_transitions_prod[2*j + 1]
                                T4[2*i +2][2*j+2] = T4[2*i +2][2*j+2] + vec4[2*j + 1]
                                T5[2*i +2][2*j+2] = T5[2*i +2][2*j+2] + pt.final_transitions_sum4[2*j + 1]                                
                                

                        pt = pt.next

                    
                    current = current.next 


        for i in range(len(extArms)):
            row_sum = sum(T1[i+1][1:])
#            print('------------------------------------\n Row = ', T1[i+1], 'with sum = ', row_sum)
            if row_sum > 0:
                T1[i+1][1:] = np.divide(T1[i+1][1:], row_sum) 
                
            row_sum = sum(T2[i+1][1:])
            if row_sum > 0:
                T2[i+1][1:] = np.divide(T2[i+1][1:], row_sum) 
            
            row_sum = sum(T3[i+1][1:])
            if row_sum > 0:
                T3[i+1][1:] = np.divide(T3[i+1][1:], row_sum) 
                
            
            row_sum = sum(T4[i+1][1:])
            if row_sum > 0:
                T4[i+1][1:] = np.divide(T4[i+1][1:], row_sum)  
                
            row_sum = sum(T5[i+1][1:])
            if row_sum > 0:
                T5[i+1][1:] = np.divide(T5[i+1][1:], row_sum)
                
        workbook = xlsxwriter.Workbook('FinalTransition-prob-method-1.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T1[row][col])
        workbook.close()   
        
        
        workbook = xlsxwriter.Workbook('FinalTransition-prob-method-2.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T2[row][col])
        workbook.close() 

        workbook = xlsxwriter.Workbook('FinalTransition-prob-method-3.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T3[row][col])
        workbook.close()
        
        workbook = xlsxwriter.Workbook('FinalTransition-prob-method-4.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T4[row][col])
        workbook.close()        
        
        workbook = xlsxwriter.Workbook('FinalTransition-prob-method-5.xlsx')
        worksheet = workbook.add_worksheet()


        for row in range(2*ChromLen + 1):
            for col in range(2*ChromLen + 1):
                worksheet.write(row, col, T5[row][col])
        workbook.close()        
            
    def bestPaths(self):
        global BestChromArray
        current = self.head
        while current != None:
            pt = current.chromlist.head
            while pt != None:
                if pt.parent == False:
                    addBestPaths(pt, pt.count, pt.chrom, '')
                pt = pt.next
            current = current.next
        firstRow = BestChromArray[0]
        BestChromArray = sorted(BestChromArray[1:], reverse=True)
        BestChromArray.insert(0, firstRow)
        workbook = xlsxwriter.Workbook('BestPaths.xlsx')
        worksheet = workbook.add_worksheet()        
        for row in range(len(BestChromArray)):
            for col in range(3):
                worksheet.write(row, col, BestChromArray[row][col])
        workbook.close() 
#----------------------------------------------
#          Main program starts here           -
#----------------------------------------------
def buildChromList():        
    Input_file = 'CALSCNAS_BRCA3_csv.csv'                                          # Name of the input csv file
    global Arms, ChromLen, BuiltLinks
    
    #Ordered linked list for maintaining the levels. An item in node is a tuple (n, pt) where pt is the pointer to the list of 
    #n CAAs. The first node points to nothing.
    
    
    
    
    
    #workbook = xlrd.open_workbook('CALSCNAS_BRCA3.xlsx')
    #worksheet = workbook.sheet_by_index(0)
    data = pd.read_csv(Input_file, encoding = "ISO-8859-1", header=None)           # Read the input dataset
    Arms = np.array(data.iloc[0][2:])                                              # Vector of CAA names
    ChromLen = len(Arms)                                                           # Length of each CAA
    noDataPts = len(data) - 1                                                      # No of CAAs given
    BuiltLinks = False
    lL = HeadList()

    
    for i in range(noDataPts):
    #    print(data.iloc[i+1])
        chromosome = np.array(data.iloc[i+1][2:], dtype='int')                                  # i^th CAA
        lL.insertChrom(chromosome)
        

#            print(data.iloc[i+1][2:])
#            print(chromosome)
#            print(np.count_nonzero(chromosome))
            

#    print('--------------------------------------------------------------')
#    print('Sequence of Arms: ', Arms)
#    print('--------------------------------------------------------------')

    lL.buildLinks()
    lL.markFirstTimeNodes()
    lL.computeTransitions()
    lL.computeGneralTransitions()
    lL.computeFinalTransitions()

    return lL
            

    
          
            
            
            
            
            
            
            
