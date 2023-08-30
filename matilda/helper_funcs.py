import matilda
import numpy
import itertools
import math


def kth_combination(current_n,current_r,current_k,binom_lookup): # RETURNS THE KTH (r-1)-SIMPLEX BY LEX ORDER
    if current_r == 0:
        return []
    result = [0 for x in range(current_r)]
    current_entry = 0
    current_value = 0
    while current_r>0:
        c = binom_lookup[current_n-1, current_r-1]
        #print c
        while not (current_k<c or current_r == current_n):
            current_k-=c
            current_value += 1
            current_n -=1
            c = binom_lookup[current_n-1, current_r-1]
        result[current_entry] = current_value
        current_r -= 1
        current_entry += 1
        current_value += 1
        current_n -=1
    return result


def kth_combination_at_degree(n,k,binomial_accum_lookup,binom_lookup): # RETURNS KTH SIMPLEX BY LEX ORDER ON ALL POSSIBLE SIMPS
    upper_index = numpy.searchsorted(binomial_accum_lookup,k+1)-1
    #print binomial_accum_lookup,upper_index,k
    k = k - binomial_accum_lookup[upper_index]
    #print k
    return kth_combination(n,upper_index+1,k,binom_lookup)


def allfaces(simplex):
    l = len(simplex)
    return [simplex[0:i] + simplex[i + 1:-1] for i in range(l)]

def mod_inverse(a): # ASSUMES a INVERTIBLE
    if matilda.modulus == 2:
        return 1
    old_r, r = matilda.modulus, a
    old_t, t = 0, 1
    while r!=0:
        q = old_r // r
        old_t, t = t, old_t - q*t
        old_r, r = r, old_r - q*r
    if old_t < 0:
        old_t += matilda.modulus
    return old_t
    


def mod_add(a, b):
    return (a + b) % matilda.modulus


def combination_kth_at_degree(facet, n, binomial_accum_lookup, binom_lookup):
    # return inverse_simplices[facet]
    return binomial_accum_lookup[len(facet) - 1] + combination_kth(facet, n, binom_lookup)


def combination_kth(facet, n, binom_lookup):  # ZERO BASED, MAX ENTRY N-1, RETURNS THE INDEX OF FACET
    """
    Returns the index of the simplex `facet`, assuming we have simplices from vertices [0,n-1], in lexicographical
    order (0 is the index of the first simplex of length the same as `facet`).

    :param facet:
    :param n:
    :param binom_lookup:
    :return:
    """
    # facet = [x+1 for x in facet]
    result = 0
    k = len(facet)
    r = binom_lookup[n, k]
    for i in range(k):
        if n - facet[i] - 1 < k - i:
            continue
        r = r - binom_lookup[n - facet[i] - 1, k - i]
    return r - 1

def distance_of_addition_vietoris(simplex, metric):  # SIMPLICES ARE LISTS, METRIC IS NP MATRIX
    if len(simplex) == 1:
        return 0.
    return max([metric[x, y] for x, y in itertools.combinations(simplex, 2)])

def distance_of_addition_weighted_vietoris(simplex, metric, weight_sq):  # SIMPLICES ARE LISTS, METRIC IS NP MATRIX
    if len(simplex) == 1:
        return 0.
    alpha = -math.inf
    for x, y in itertools.combinations(simplex, 2):
        d = metric[x,y]
        b = abs((weight_sq[y] - weight_sq[x])/d)
        alpha =  max(alpha,math.sqrt(d**2+b**2+2*weight_sq[x]+2*weight_sq[y])/2)
    return alpha



def naive_intersection_sorted(array1, array2):  # intersection of two sorted numpy arrays with no repeated values
    i = 0
    j = 0
    result = []
    while i<array1.size and j<array2.size:
        while array1[i]<array2[j]:
            i+=1
        while array2[j]<array1[i]:
            j+=1
        if array1[i]==array2[j]:
            result.append(array1[i])
            i+=1
            j+=1
    return numpy.array(result)


def power_metric_distance(p1, p, metric):
    """

    :param p1: point
    :param p2: point-set p
    :return: power distance between points p1 and point set p: Eqn (2) of Ref [1]
    """

    return min(metric[p1][p])

def s_step_function(alpha,lambdas,epsilon):
    """

    :param alpha: filtration radius
    :param lambdas:
    :param epsilon:
    :return: step function
    """
    if alpha <= lambdas/epsilon:
        return 0
    elif alpha < lambdas/(epsilon*(1-epsilon)):
        return alpha-lambdas/epsilon
    else:
        return epsilon*alpha



def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def piecewise_linear_lambda(v,xs,ys):   
    def in_interval(x,a,b):
        return (a<=x) &(x<=b)
    def interp_lambda(x1,y1,x2,y2):
        if x1>x2:
            raise ValueError("x coordinates must be sorted.",x1,x2)
        if y1 == y2:
            return lambda x: y1
        else:
            return lambda x:numpy.interp(x,[x1,x2],[y1,y2])
    conditions = [in_interval(v,xs[i],xs[i+1]) for i in range(len(xs)-1)]
    definitions = [interp_lambda(xs[i],ys[i],xs[i+1],ys[i+1]) for i in range(len(xs)-1)]
    return lambda x: numpy.piecewise(x,conditions,definitions)
        
        
def sample_from_uniform_func_random(point_cloud,subsample,epsilon):
    dims = len(point_cloud[0])
    max_in_each_dim = numpy.ones(dims) * (-1) * math.inf
    min_in_each_dim = numpy.ones(dims) * math.inf
    for ind,num in enumerate(point_cloud):  # no of points
        for dim_iter,dim_val in enumerate(num):
            if dim_val<min_in_each_dim[dim_iter]:
                min_in_each_dim[dim_iter] = dim_val
            if dim_val>max_in_each_dim[dim_iter]:
                max_in_each_dim[dim_iter] = dim_val
    min_in_each_dim = numpy.array(min_in_each_dim) - epsilon
    max_in_each_dim = numpy.array(max_in_each_dim) + epsilon

    # print ('max_in_each_dim:', max_in_each_dim)
    # print('min_in_each_dim:', min_in_each_dim)
    subsample_points = []
    for _ in range(subsample):
        dist = epsilon + 1
        thispoint = []
        while (dist > epsilon):
            thispoint = []
            for dim_iter in range(dims):
                thispoint.append(numpy.random.uniform(min_in_each_dim[dim_iter], max_in_each_dim[dim_iter], 1)[0])

            for points in point_cloud:

                dist = numpy.sum((numpy.array(thispoint) - numpy.array(points)) ** 2)
                if dist<=epsilon:
                    break
                # print(dd, points, emin, qq)
            # exit(0)

        subsample_points.append(thispoint)
    return subsample_points


def choose_distant_objects(D):
    ob =  random.randint(0,N)
    # print(D,ob,len(D))
    # print(D[ob,:])
    oa = numpy.argmax(D[ob,:])
    ob = numpy.argmax(D[oa,:])
    return [oa,ob]

def distortion_parameter(M1,M2):
    #M1: square matrix of distances
    #M2: Nxk matrix euclidean points: each row is a point
    if len(M1)!=len(M2):
        raise ValueError("dimension mismatch")
    sum_num = 0
    sum_dem = 0
    for i in range(len(M1)):
        for j in range(i+1,len(M1)):

            d1 = 0
            for elem in range(len(M2[0])):
                d1 += (M2[i,elem]-M2[j,elem])**2
            d1 = math.sqrt(d1)
            d = M1[i, j]
            sum_num += (d1-d)**2
            sum_dem += d**2
    return math.sqrt(sum_num/sum_dem)

def fast_map(k,D):
    global col_N
    global X
    global pa
    if k<=0:
        return 1
    else:
        col_N += 1

    oab = choose_distant_objects(D)
    pa[0, col_N] = oab[0]
    pa[1, col_N] = oab[1]
    if D[oab[0],oab[1]]==0:
        for i in range(N):
            X[i,col_N] = 0
        return 1
    for i in range(N):
        x = (D[oab[0],i]**2 + D[oab[0],oab[1]]**2 - D[oab[1],i]**2)/(2*D[oab[0],oab[1]])
        X[i,col_N] = x
    D_new = numpy.zeros([N,N])
    for i in range(N):
        for j in range(N):
            # print('i',i,'j',j,'D[i,j]',D[i,j],'X',X[i,col_N],X[j,col_N],X[i,col_N]-X[j,col_N])
            val = D[i,j]**2-(X[i,col_N]-X[j,col_N])**2
            if val<1e-10:
                D_new[i,j] = 0
            else:
                D_new[i,j] = math.sqrt(val)

    fast_map(k-1,D_new)

    return 0

def sample_new_thm(point_cloud,subsample_size,epsilon):
    dims = len(point_cloud[0])

    norm_hist = []
    for i in range(len(point_cloud)):
        # for j in range(i+1,len(point_cloud)):
        norm_hist.append(numpy.linalg.norm(point_cloud[i,:]))
        #numpy.sqrt(numpy.sum((point_cloud[i,:])**2)))#

    print('norm hist calculated, dim: ',dims)

    subsample_pt =[]
    subsample_ind =[]
    for ind in range(subsample_size):
        found = False
        subsample_local_ind = []
        print(' going on ')
        while found==False:
            gen_point = numpy.random.normal(0, 1, dims)  # generate n-numbers for each dimension, need to scale back now
            x = numpy.random.uniform(min(norm_hist)-epsilon,max(norm_hist)+epsilon)
            gen_point = gen_point/numpy.linalg.norm(gen_point)
            gen_point = gen_point*x
            print ('generated p',gen_point)
            # min_p = math.inf
            # min_d = -1
            for inds,pointg in enumerate(point_cloud):
                dist = numpy.sqrt(numpy.linalg.norm(pointg-gen_point))
                # numpy.sqrt(numpy.sum((pointg - gen_point)**2))

                #
                if dist<=epsilon:
                    found = True
                    # if dist<min_p:
                        # min_p = dist
                        # min_d = inds
                    subsample_local_ind.append(inds)
                    # print(dist,epsilon)
                else:
                    print('not this one',dist,epsilon)

            if(found==True):
                subsample_local_ind=numpy.array(subsample_local_ind)
                subsample_ind.append(subsample_local_ind)
                subsample_pt.append(gen_point)
                print('found point ',ind,' with ind ',len(subsample_ind),subsample_local_ind)#,' dim ',len(subsample_local_ind))
                # print('closest points ',subsample_local_ind)
                # print('gen_point ',gen_point)sa
            print ('false',found,gen_point)

    return subsample_pt,subsample_ind


def param_parser(str_param, default_value, kwargs, required=False):
    """
    Simple kwargs parser for optional parameters with default values.
    """
    if str_param in kwargs.keys():
        param = kwargs[str_param]
        if param is None:
            param = default_value
    elif required:
        raise ValueError("Parameter "+str_param+" is required.")
    else:
        param = default_value
    return param