import numpy as np

def table2(filename):
    data_uv = np.loadtxt(filename, delimiter=",")
    if data_uv.shape == (len(data_uv),):
        x_ = data_uv[0]
        y_ = data_uv[1]
    else:
        x_ = data_uv[:, 0]
        y_ = data_uv[:, 1]    
    return x_, y_ 

def table3(filename):
    data_uv = np.loadtxt(filename, delimiter=",")
    if data_uv.shape == (len(data_uv),):
        x_ = data_uv[0]
        y_ = data_uv[1]
        z_ = data_uv[2]
    else:
        x_ = data_uv[:, 0]
        y_ = data_uv[:, 1]  
        z_ = data_uv[:, 2]        
    return x_, y_, z_ 

def table4(filename):
    data_uv = np.loadtxt(filename, delimiter=",")
    if data_uv.shape == (len(data_uv),):
        x_ = data_uv[0]
        y_ = data_uv[1]
        z_ = data_uv[2]
        t_ = data_uv[3]
    else:
        x_ = data_uv[:, 0]
        y_ = data_uv[:, 1]  
        z_ = data_uv[:, 2]  
        t_ = data_uv[:, 3] 
    return x_, y_, z_, t_       

def table5(filename):
    data_uv = np.loadtxt(filename, delimiter=",")
    if data_uv.shape == (len(data_uv),):
        x_ = data_uv[0]
        y_ = data_uv[1]
        z_ = data_uv[2]
        w_ = data_uv[3]
        t_ = data_uv[4]
    else:
        x_ = data_uv[:, 0]
        y_ = data_uv[:, 1]  
        z_ = data_uv[:, 2]  
        w_ = data_uv[:, 3] 
        t_ = data_uv[:, 4] 
    return x_, y_, z_, w_, t_ 
    
    
def table6(filename):
    data_uv = np.loadtxt(filename, delimiter=",")
    if data_uv.shape == (len(data_uv),):
        x_ = data_uv[0]
        y_ = data_uv[1]
        z_ = data_uv[2]
        u_ = data_uv[3]
        v_ = data_uv[4]
        w_ = data_uv[5]
    else:
        x_ = data_uv[:, 0]
        y_ = data_uv[:, 1]  
        z_ = data_uv[:, 2]  
        u_ = data_uv[:, 3] 
        v_ = data_uv[:, 4]  
        w_ = data_uv[:, 5] 
    return x_, y_, z_, u_, v_, w_ 

def table7(filename):
    data_uv = np.loadtxt(filename, delimiter=",")
    if data_uv.shape == (len(data_uv),):
        x_ = data_uv[0]
        y_ = data_uv[1]
        z_ = data_uv[2]
        u_ = data_uv[3]
        v_ = data_uv[4]
        w_ = data_uv[5]
        t_ = data_uv[6]
    else:
        x_ = data_uv[:, 0]
        y_ = data_uv[:, 1]  
        z_ = data_uv[:, 2]  
        u_ = data_uv[:, 3] 
        v_ = data_uv[:, 4]  
        w_ = data_uv[:, 5] 
        t_ = data_uv[:, 6] 
    return x_, y_, z_, u_, v_, w_, t_     

