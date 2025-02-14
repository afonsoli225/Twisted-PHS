import math
from sage.modules.free_module_integer import IntegerLattice
from sage.all import *
from lattice import *
from number_field import *
from twphs_algo import *
import numpy as np
from sage.modules.free_module_integer import IntegerLattice




def svp_exact_Real(Mv, work_prec=0):
    _R         = Mv.base_ring().fraction_field();
    _log_scale = _R.precision() if (work_prec == 0) else work_prec;
    _scale     = Integer(2)**(_log_scale);
    _MZ        = Mv.apply_map(lambda _mij:Integer(floor(_scale*_mij)))
    L=IntegerLattice(_MZ)
    v=L.shortest_vector()
    return v/((Integer(2)**_log_scale))


def Short_vector_ideal(P):
    K.<a> = NumberField(P)
    OK = K.ring_of_integers()
    signature_K = K.signature()
    print("We have " + str(signature_K[0]) + " real embeddings and " + str(signature_K[1]) + " pairs of complex embeddings")
    r_K=signature_K[0]
    s_K=signature_K[1]
    

    #the rank of Unit group is dim_Lambda_K
    dim_Lambda_K=r_K+s_K-1
    UK=UnitGroup(K); UK
    
    #define the class group of K, then search for the prime ideals such their classes generate the class group of K
    #it is for now a minimal set
    G = K.class_group(); 
    
    #the representative_prime function return a prime ideal in an ideal class.
    #gets prime ideals such that their classes generate the class group of K
    G.gens()
    S_0=[c.representative_prime() for c in G.gens()]
    #define the S_0 unit group
    SUK = UnitGroup(K,S=tuple(S_0)); SUK
    
    #size of S_0 
    k_0=len(S_0)
    
    #Returns the inifinite places of the number field K
    p_inf=get_inf_places(K, b_prec=fp.BIT_PREC_DEFAULT)
    print ("size of S_0 is ",len(S_0))

    A=SUK.fundamental_units()

    for i in range(0,r_K+s_K-1):
        un=UK.fundamental_units()[i]
        A.remove(un)
    
    Log_S_unit=twphs_get_matrix(UK.fundamental_units(), A, p_inf, S_0, method='TW', b_prec=fp.BIT_PREC_DEFAULT)


    print("This is the matrix basis of the log-S-unit lattice\n",Log_S_unit.n(5))
    print("and its volume is ",(Log_S_unit.det().abs()).n(20))
    print("This is the shortest vector of the log-S-unit lattice")
    svp_Log=svp_exact_Real(Log_S_unit)
    print(svp_Log.n(5))
    print("which is of the rounded norm ")
    print(svp_Log.norm().n(20))
    d=dim_Lambda_K+len(S_0)
    print("Dimension of S-unit lattice is ",d)
    print("The gaussian heuristic give the expected length of the shortest vector to be around")
    gh=(sqrt(d / (2 * pi * e)) * ((Log_S_unit.det().abs())^(1/d))).n(20)
    print(gh)
    print("Minkowski bound")
    print((sqrt(d)*((Log_S_unit.det().abs())^(1/d))).n(5))
    

    




