    
def format_index(i):
        
    string = ""
    for j in range(3 - len(str(i))):
        string += "0"
    
    string += str(i)
    return string

def diff(li1, li2): 

    dif = []
    
    for i in range(len(li1)):
        dif.append(li1[i] - li2[i])
        
    
    return dif

