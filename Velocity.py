"""General operation with velocity of secondary element"""

def VelocityData(tau, F1):
    if v in globals():
        print('variable of velocity exists')
    elif sk in globals():
        print('variable of velocity exists')
    else:
        print('variables of slip and velocity exist')






def VelocityData2(tau, F1):
    Vsinhr = 2 * tau * F1   # Sinhronus velocity
    "Check for the existence of variables "
    try:
        v
    except NameError:
        try:
            sk
        except NameError:
            print('v and sk are not defined')
        else:
            print('v is not defined sk is defined')
            #v = (1 - sk) * Vsinhr  # Operating velocity
    else:
        #sk = 1 - v/Vsinhr
        print('v is defined')



def VelocityData4(Flag, tau, F1):
    if Flag == 'sk':
        global sk
        Vsinhr = 2 * tau * F1  # Sinhronus velocity
        v = (1 - sk) * Vsinhr  # Operating velocity
        return sk, v
    else:
        global v
        Vsinhr = 2 * tau * F1  # Sinhronus velocity
        sk = 1 - v/Vsinhr
        return sk, v

