#This file will contain all the strings for the names and units to go in plots in one place. Everything will be read from here. 

def get_parameter_info(model="jorg"):
    '''Returns the info dictionary. For each key, the first element is the name and the second is the unit.
    
    Parameters
    ----------
    model - 'jorg' or 'cmem' 
    '''
    
    if model == 'jorg': 
        info = {
        "rmp":(r"${r_0}^{mp}$", "R"+r"$_E$"), 
        "rbs":(r"${r_0}^{bs}$", "R"+r"$_E$"), 
        "A1":(r"$A_1$", "eV cm"+r"$^{-3}$"+" s"+r"$^{-1}$"), 
        "A2":(r"$A_2$", "eV cm"+r"$^{-3}$"+" s"+r"$^{-1}$"), 
        "B":("B", "eV cm"+r"$^{-3}$"+" s"+r"$^{-1}$"), 
        "alpha":(r"$\alpha$", ""), 
        "beta":(r"$\beta$", ""), 
        "ay_mp":(r"${\alpha_y}^{mp}$", ""), 
        "az_mp":(r"${\alpha_z}^{mp}$", ""), 
        "ay_bs":(r"${\alpha_y}^{bs}$", ""), 
        "az_bs":(r"${\alpha_z}^{bs}$", ""),
        }
        return info 
    
    elif model == 'cmem': 
        info = { 
        "p0":(r"$p_0$", ""), 
        "rbs":(r"${r_0}^{bs}$", "R"+r"$_E$"), 
        "A1":(r"$A_1$", "eV cm"+r"$^{-3}$"+" s"+r"$^{-1}$"), 
        "A2":(r"$A_2$", "eV cm"+r"$^{-3}$"+" s"+r"$^{-1}$"), 
        "B":("B", ""), 
        "alpha":(r"$\alpha$", ""), 
        "beta":(r"$\beta$", ""), 
        "p1":(r"$p_1$", ""),
        "p2":(r"$p_2$", ""),
        "p3":(r"$p_3$", ""),
        "ay_bs":(r"${\alpha_y}^{bs}$", ""), 
        "az_bs":(r"${\alpha_z}^{bs}$", ""),
        }
        return info 
        
    else: 
        raise ValueError("Invalid model selected. 'jorg' or 'cmem'") 
    

