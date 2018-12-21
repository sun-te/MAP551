import numpy as np
import pandas as pd
import os,re
def write_input(list_input):
    assert(len(list_input)==9);
    input_data=''
    for i in range(9):
        input_data+=str(list_input[i])+'\t \n'
    with open("input.dat",'r+') as f:
        read_data=f.read()
        f.seek(0)
        f.truncate()
        f.write(input_data)
    data=''
    variable=['neq','tini','tend','nt','xmin','xmax','nx','method','tol']
    count=0
    with open("input.dat",'r') as f:
        for line in f.readlines():
            data+=str(variable[count])+' = '+line
            count+=1
    print(data)
    return     

def String2Float(x):
    str_list=x.split(' ');
    ans=[];
    for i in range(len(str_list)):
        try:
            if(len(str_list[i])>0 ):
                ans.append(float(str_list[i]))
        except:
            print(str_list[i], " can not be converted in to a float")
    return ans

def read_data(file_name):
    df_data=pd.read_table(file_name,sep=r",",engine='python',skip_blank_lines=1)
    num_temp=String2Float(df_data.columns[0])
    col=df_data.columns[0]
    data_plot=[]
    for i in range(df_data.shape[0]):
        list_tmp=String2Float(df_data[col][i])
        data_plot.append(list_tmp)
    data_plot=np.array(data_plot).T
    return data_plot, num_temp
