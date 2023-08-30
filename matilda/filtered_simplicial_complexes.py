import numpy
import itertools
# from matilda import helper_funcs
import math
import random
import matildacpp


class FilteredSimplicialComplex(object):
    """
    Class for modelling a filtered (finite) simplicial complex.

    Attributes
    ----------
    simplices
        List containing all simplices in the complex, regardless of filtration value.
        Simplices are themselves numpy arrays of integers.

    dimension: int
        Dimension of largest simplex. Equal to length of largest simplex minus one.
        number_of_vertices: Number of vertices present in the complex, regardless of filtration value.

    simplices_indices
        List of ints encoding a strict ordering of `simplices` according to filtration values of each simplex. It
        requires simplices to have larger indices than their faces. This means that to iterate over simplices
        in filtration order, we do so by iterating over simplices_indices, so simplices[simplices_indices[i]]
        is the ith simplex added to the filtration.

    appears_at
        List of all filtration values for all simplices in the complex. Also indexed by simplices_indices.
        
    """
    def __init__(self, 
                    dimension=0, 
                    simplices=None,
                    simplices_indices=None, 
                    appears_at=None):
        """
        Default initializer. Sets all main attributes of simplex to default values.

        Returns
        -------
        void
        """
        self.dimension = dimension
        if simplices_indices is None:
            self.simplices_indices=[]
        else:
            self.simplices_indices = simplices_indices
        if appears_at is None:
            self.appears_at = []
        else:            
            self.appears_at = appears_at
        if simplices is None:
            self.simplices = []
        else:
            self.simplices = simplices

    def print_simplices(self):
        """
        Prints vertex indices in the simplex along with the filtration value, sorted by filtration value.

        Returns
        -------
        void
        """
        for s_index in self.simplices_indices:
            print('simplex:',self.simplices[s_index],' filtration_value: ',self.appears_at[s_index])#,' index:',self.simplices_indices[s_index])

    def __str__(self):
        return "".join(['simplex: '+str(self.simplices[s_index])+' filtration_value: '+str(self.appears_at[s_index])+"\n" \
                    for s_index in self.simplices_indices])



    def has_valid_filtration(self):
        """
        Checks if any simplices appear before their faces. Note that this method is very time consuming.
        Checks for sanity of simplices 1. filtration value(appears_at) of faces 2. simplices_indices

        Returns
        -------
        self : object
            Returns self.
        """
        current_appears_at = -1
        for order_index,order_item in enumerate(self.simplices_indices):
            simplex_item = self.simplices[order_item]
            if self.appears_at[order_item]<current_appears_at:
                raise ValueError('Filtration value has decreasing value. ')
            else:
                current_appears_at = self.appears_at[order_item]

            if len(simplex_item) >= 2:  # check existence of corresponding (n-1)-simplices
                faces = numpy.array(list(itertools.combinations(simplex_item,len(simplex_item)-1)))
                for face_item in faces:

                    is_face = False
                    nn = numpy.asarray(face_item)

                    # print ("*")
                    for index,s_item in enumerate(self.simplices_indices): # going till indices of current simplex  ,range(order_item+1)
                        s_simplex_item = self.simplices[s_item]
                        # print (s_simplex_item)
                        if numpy.array_equal(s_simplex_item,nn):
                            if self.appears_at[order_item] < self.appears_at[s_item]:
                                # filtration value of coface is less
                                print (self.simplices[order_item],self.appears_at[order_item],order_index,self.simplices[s_item],self.appears_at[s_item],index)
                                raise ValueError('Simplex:',self.simplices[order_item],' has lower filtration value that its face:',self.simplices[s_item],'. Not permitted')
                            if  order_index< index:
                                print (self.simplices[order_item], self.appears_at[order_item], order_index,
                                       self.simplices[s_item], self.appears_at[s_item], index)
                                raise ValueError('Simplex:', self.simplices[order_item],
                                                 ' has lower simplices indices value that its face:', self.simplices[s_item],
                                                 '. Not permitted')
                            is_face = True
                            break

                    if not is_face:
                        raise ValueError("Simplex mismatch:", nn, " does not exist for ", simplex_item)


        print('Test successful, filtration is valid')
        return self
        
    def add_simplex(self, simplex, filtration_value, skip_checks = True):
        """
        Adds a simplex to the complex at the specified filtration value.

        Parameters
        ----------
        simplex: List
            List of vertices representing a simplex

        filtration_value: float
            Filtration value of addition of the simplex to the complex

        skip_checks: bool
            If `True`, it is assumed that the simplex is being added respecting the order induced by filtration values
            and dimension. If `False`, consistency checks are performed, raising an Exception
            if problems are found. Note that enabling checks increases run time.
            Default is `True`.

        Returns
        -------
        self : object
            Returns self.
        """
        if skip_checks:
            self.simplices.append(numpy.sort(numpy.array(simplex)))
            self.appears_at.append(filtration_value)
            self.simplices_indices.append(len(self.appears_at)-1)
        else:
            self.simplices.append(numpy.sort(numpy.array(simplex)))
            self.appears_at.append(filtration_value)
            self.simplices_indices.append(len(self.appears_at) - 1)
            has_valid_filtration()
            ### TODO: ADD MORE CHECKS
            assert(False)
        return self

    def check_duplicate_simplex(self,simplex):
        """
        Checks if simplex exist in complex

        Parameters
        ----------
        simplex: List
            List of vertices representing a simplex

        Returns
        -------
        bool
            True if simplex exists, False otherwise
        """
        if not numpy.array_equal(numpy.sort(simplex),simplex): #vertex indices are not in correct order
            return False
        for item in self.simplices_indices:
            if not len(simplex) == len(self.simplices[item]):
                continue
            elif numpy.array_equal(simplex, self.simplices[item]):
                return True

        return False


    def delete_simplex(self, simplex, skip_checks = True):
        """
        Deletes the specified simplex from the simplicial complex
        Deletes all of its co-faces to maintain sanity of the simplicial complex

        Parameters
        ----------
        simplex: List
            List of vertices representing a simplex

        skip_checks: bool
            If `True`, it is assumed that the simplex is being added respecting the order induced by filtration values
            and dimension. If `False`, consistency checks are performed, raising an Exception
            if problems are found. Note that enabling checks increases run time.
            Default is `True`.
        
        Returns
        -------
        void
            Error if simplex does not exist
            
        """

        found_simplex = False
        for index,order_item in enumerate(self.simplices_indices):
            simplex_item = self.simplices[order_item]
            if numpy.array_equal(simplex_item,simplex): # mark the simplex for deletion
                self.simplices[order_item] = numpy.array([])
                self.appears_at[order_item] = -1
                found_simplex = True
                continue

            if len(simplex_item) > len(simplex): # this simplex_item can be a face of simplex
                comb_face = (list(itertools.combinations(simplex_item,len(simplex)))) # all faces of simplex_item with dimension of simplex
                comb_face = [list(elem) for elem in comb_face]

                if simplex.tolist() in comb_face: #current simplex has a simplex to be removed as a face. Need to be removed as well
                    if not found_simplex:
                        raise ValueError("Invalid Filtration for: ",simplex,' Co-face: ',simplex_item, comb_face)
                    self.simplices[order_item] = numpy.array([])
                    self.appears_at[order_item] = -1

        # Actually delete the simplices
        self.appears_at = numpy.array([num for num in self.appears_at if num >= 0])
        self.simplices = numpy.array([num for num in self.simplices if len(num) > 0])
        self.simplices_indices = numpy.argsort(self.appears_at, kind="mergesort")  # USES mergesort BECAUSE IT'S STABLE


        if not found_simplex:
            raise ValueError("Simplex not found: ",simplex)
        if not skip_checks:
            has_valid_filtration()

        return self


    def make_consistent_indices(self):
        """
        Makes simplex indices consistent with their filtration values and dimensions and sets the correct dimension of
        complex. Assumes attributes `simplices` and `appears_at` represent a correct filtration.

        Returns
        -------
        self : object
            Returns self.
        """
        lengths = [len(x) for x in self.simplices]
        self.dimension = max(lengths)

        self.simplices_indices = numpy.argsort(self.appears_at,kind="stable")

        self.simplices = [self.simplices[i] for i in self.simplices_indices]
        self.appears_at = [self.appears_at[i] for i in self.simplices_indices]

        self.simplices_indices = [i for i in range(len(self.simplices))]

        return self

    def construct_vietoris_from_metric(self, matrix, dimension, upper_bound):
        """
        Constructs a filtered Vietoris-Rips simplicial complex from a given metric space, limited to points at most
        separated by a given distance. Based on incremental VR construction from [1]_.

        Parameters
        ----------
        matrix: Numpy 2D-array
            Metric matrix
        dimension: int
            Maximum dimension of simplices in complex
        upper_bound: float
            Maximum distance for edges

        Returns
        -------
        self : object
            Returns self.

        References
        ----------
        .. [1] Zomorodian, Afra. "Fast construction of the Vietoris-Rips complex." Computers \& Graphics 34 (2010): 263-271.
        """
        # TODO: Add automatic parsing of point cloud matrices.
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("Expected square matrix for input metric.")
        self.dimension = dimension
        temp = matildacpp.FilteredSimplicialComplex()
        temp.construct_vietoris_from_metric(matrix, dimension, upper_bound)
        self.simplices = temp.simplices
        self.simplices_indices = temp.simplices_indices
        self.appears_at = temp.appears_at
        return self.make_consistent_indices()


    def construct_full_weighted_vietoris_from_metric(self, matrix, weight_sq, dimension, verbose=False):  # DIMENSION IS LENGTH+1 OF LARGEST SIMP
        """
        Constructs a filtered weighted Vietoris-Rips simplicial complex from a given metric space.

        Parameters
        ----------
        matrix: Numpy 2D-array
            Metric matrix
        dimension: int
            Maximum simplex dimension of simplices in complex
        weight_sq: list
            Weights of vertices, squared

        Returns
        -------
        self : object
            Returns self.

        """
        self.number_of_vertices = matrix.shape[0]
        self.dimension = dimension
        #print "Initializing simplices... ",
        for i in range(1, dimension + 2):
            self.simplices.extend([numpy.array(x,dtype="int") for x in itertools.combinations(range(matrix.shape[0]), i)])
        #print "done. \n%d simplices initialized.\nDetermining filtration values..."%len(self.simplices)
        self.appears_at = numpy.zeros(len(self.simplices))
        i = 0
        edges_appears_at = {}
        for d in range(2):
            simplex_length = d + 1
            for s in itertools.combinations(range(matrix.shape[0]), simplex_length):
                self.appears_at[i] = helper_funcs.distance_of_addition_weighted_vietoris(s, matrix,weight_sq)
                edges_appears_at[s] = self.appears_at[i]
                i += 1
                if i % 100000 == 0 and verbose:
                    print (i, "of ",len(self.simplices))
        for d in range(2,self.dimension + 1):
            simplex_length = d + 1
            for s in itertools.combinations(range(matrix.shape[0]), simplex_length):
                self.appears_at[i] = max(edges_appears_at[t] for t in itertools.combinations(s,2))
                i += 1
                if i % 100000 == 0 and verbose:
                    print (i, "of ",len(self.simplices))
        #print "done. \nSorting simplices by filtration value...",
        self.simplices_indices = numpy.argsort(self.appears_at, kind="mergesort") # USES mergesort BECAUSE IT'S STABLE

    def simplex_index_length(self, simplex_index):
        result = numpy.searchsorted(self.binomial_accum_lookup, simplex_index + 1)
        return result

    def construct_witness_complex(self,v,dimension,matrix,type='random',ratio=20,radius=math.inf):
        """
        Constructs a Witness Complex from a given metric space and landmarks, limited to points at most
        separated by a given distance. Based on [1]_.

        Parameters
        ----------
        v: int
        Non-negative integer parameter defining the class of Witness complex [generally 0-2]

        dimension: int
            Dimension of largest simplex. Equal to length of largest simplex minus one.

        matrix: Numpy 2D-array
            Metric matrix


        type: String: 'random' or 'maxmin'
            Defines how landmarks are chosen

        ratio: float
            Ratio of original vertices to random vertices(rounded up). Initially it is 20 ( as suggested by paper)

        radius: float
            Maximum distance for edges

        Returns
        -------
        self : object
            Returns self.

        References
        ----------
        [1] De Silva and Carlsson Topological estimation using witness complexes

       """

        def random_witness_generator(matrix, ratio,witnessD):
            """
            Given original matrix, randomly finds N/20 witness vertices, where N is original number of vertices

            Parameters
            ----------
            matrix: Numpy 2D-array
                Metric matrix

            Returns
            -------
            List
                The list which gives indices of points chosen as landmarks

            """
            n = int(math.ceil(len(matrix) / ratio))  # number of witness

            N = len(matrix)
            # indices_of_points = numpy.array(range(N))
            witness_distance = numpy.zeros((n, N))
            maxmin = numpy.sort(random.sample(range(0, N), n))  # set of landmarks
            # print (maxmin)
            for index_maxmin, item_maxmin in enumerate(maxmin):
                for index_matrix in range(N):
                    witnessD[index_maxmin,index_matrix] = matrix[index_matrix][item_maxmin]  # matrix is symmetric

            # print('size here: ',len(self.witness_distance),len(self.witness_distance[0]))
            # witnessD = witness_distance
            return maxmin

        def maxmin_witness_generator(matrix, ratio, witnessD):
            """
            Given original matrix, finds N/20 witness vertices, where N is original number of vertices
            Iteratively finds points by maximising minimum distance to all other points. (Refer to paper)

            Parameters
            ----------
            matrix: Numpy 2D-array
                Metric matrix

            Returns
            -------
            List
                The list which gives indices of points chosen as landmarks

            """
            N = len(matrix)  # total number of datapoints, same variable name as paper
            indices_of_points = numpy.array(range(N))
            # print('indices_of_points',indices_of_points)
            maxmin = numpy.array(random.sample(range(0, N), 1))  # Generate first point
            n = len(maxmin)  # total number of landmarks: initially 1
            # print ('maxmin', maxmin)


            while n < N / ratio:  # Iteratively add points
                # print('indices_of_points2: ', indices_of_points, 'maxmin:',maxmin)
                z = []
                z_pos = []
                for z_index in indices_of_points:  # iterate through all remaining points
                    now_list = list(zip(itertools.repeat(z_index), maxmin))  # (z,l1),..(z,l_(i-1))
                    dzli = list(
                        matrix[now_list[x][0]][now_list[x][1]] for x in range(len(now_list)))  # D(z,li) section 2.3

                    z.append(min(dzli))
                    z_pos.append(z_index)

                val_to_append = z_pos[z.index(max(z))]
                # index_of_val = numpy.argwhere(indices_of_points == val_to_append)[0][0]
                maxmin = numpy.append(maxmin, val_to_append)
                n = len(maxmin)

            # print('indices_of_points2', indices_of_points,'total len:',len(matrix),'now chosen:',len(indices_of_points))
            witness_distance = numpy.zeros((n, N))
            for index_maxmin, item_maxmin in enumerate(maxmin):
                for index_matrix, item_matrix in enumerate(indices_of_points):
                    witness_distance[index_maxmin,index_matrix] = matrix[item_matrix][item_maxmin]  # matrix is symmetric


            witnessD[:] = witness_distance
            return maxmin
        ################### end maxmin ###########################


        n = int(math.ceil(len(matrix) / ratio))  # number of witness
        N = len(matrix)
        witness_distance = numpy.zeros((n, N))
        if type=='random':
            random_witness_generator(matrix, ratio,witness_distance)

        else:
            maxmin_witness_generator(matrix,ratio,witness_distance)

        mi = numpy.zeros(N)
        self.dimension = dimension
        if self.dimension<0:
            return self
        if v != 0:
            for i in range(N):
                mi[i] = numpy.partition(witness_distance[:,i],v)[v] #pick vth smallest entry of column D
        if self.dimension==0:
            return self

        for zero_simplex in range(n):
            self.add_simplex([zero_simplex],0)

        E = numpy.zeros((n,n))
        # print (self.witness_distance)
        for i in range(n):
            for j in range(i+1,n):
                # Find min over (k=1..N) {max(D(i,k),D(k,j))-mk}
                # Dlin = numpy.reshape(self.witness_distance,n*N)
                min_find = math.inf

                for k in range(N):
                    max_m = max(witness_distance[i,k],witness_distance[j,k])-mi[k] #Row major wise find max( D(i,k),D(k,j)) - mi(k)
                    # print (witness_distance[i,k],witness_distance[j,k],mi[k],i,j,k,min_find)
                    # print (witness_distance)
                    # print (mi)

                    if max_m < min_find:
                        min_find = max_m
                    if min_find<0:
                        raise ValueError('metric cannot be negative')

                E[i,j] = min_find
                E[j,i] = min_find
                if E[i,j]<=radius:
                    self.add_simplex([i,j],E[i,j])
                #     print (i,j,E[i,j])
                # print("***************")
                # exit(0)
                # print (i,j,min_find)


        # print (E)
        if self.dimension==1:
            self.simplices_indices = numpy.argsort(self.appears_at, kind="mergesort")
            return self


        # Dont need self witness distance after this

        # all edges are now created
        # all simplices will now be added dimension wise
        indices_of_landmarks = numpy.array(range(n))
        # print ('indices:',indices_of_landmarks)

        for current_dim in range(2,self.dimension+1): # start with triangles to highest simplices
            for comb in list(itertools.combinations(indices_of_landmarks,current_dim+1)): # iterate through all simplex in this dimension

                two_simplices = list(itertools.combinations(comb,2))
                lin_two_simplices =[E[x[0],x[1]] for x in two_simplices ]#
                #[two_simplices[x][0]*n+two_simplices[x][1] for x in range(len(two_simplices))]
                # Elin = numpy.reshape(E,n*n)

                if max(lin_two_simplices)<=radius:
                    self.add_simplex(comb,max(lin_two_simplices))

        # self.print_simplex()
        self.simplices_indices = numpy.argsort(self.appears_at, kind="mergesort")

        # print ('witness complex created')
        return E

    def greedy_permutation(self,matrix, lambdal):
        """
        Finds greedy permutation of a set of points and their correspondin lamda value

        Parameters
        ----------
        matrix: Numpy 2D array
            Input point metric
        lambdal: list
            The list which gives indices of points chosen as landmarks
            The list is filled with lamda_i giving the \lamda^p_i value for each point: monotonically decreasing with |lambda_0 = inf

        Returns
        -------
        Matrix
            New distance matrix such that 0,1,2,3,.... points correspond to (p0,p1,p2,....) according to greedy permutation

        """
        vertices_indices = [0] #append the first vertex
        vertices_left = []
        vertices_left.extend(range(1,len(matrix))) # rest of the points
        lambda_list = [math.inf] #lambda_0


        for index_range in range(1,len(matrix)):    # Find the next p_i, needs to iterate #points times

            hold_distance_metric = numpy.zeros(len(matrix) - index_range)
            for ind,counter in enumerate(vertices_left):
                hold_distance_metric[ind] = (helper_funcs.power_metric_distance(counter,vertices_indices,matrix))

            max_ind  = numpy.argmax(hold_distance_metric)

            lambda_list.append(hold_distance_metric[max_ind])
            vertices_indices.append(vertices_left[max_ind])

            vertices_left.remove(vertices_left[max_ind])

            # print (' outerloop: max_val: ', max_val, ' index of max in hold: ', hold_distance_metric.index(max_val))
            # print(' hold_dist_met: ', hold_distance_metric, ' vert_indices ', vertices_indices, ' vert_left ', vertices_left)




        matrix_new = numpy.zeros((len(matrix),len(matrix)))
        for i in range(len(matrix_new)):
            for j in range(len(matrix_new)):
                matrix_new[i][j] = matrix[vertices_indices[i]][vertices_indices[j]] # Physically changing the variable ordering

        matrix = matrix_new

        lambdal[:] = lambda_list

        return  matrix



    def construct_sparse_rips_complex(self,matrix,dimension,epsilon,radius=math.inf):
        """
        Constructs a Sparse Rips Complex from a given metric space, limited to points at most
        separated by a given distance. Based on [1]_.

        Parameters
        ----------
        matrix: Numpy 2D array
            input metric
        dimension: int
            maximum dimension to generate simplices
        epsilon(0-1): float
            sparsity level, bottleneck distance between rips and sparse rips epsilon=0(dense), epsilon=1(sparse)
        radius: float
            Maximum distance for edges

        Returns
        -------
        self : object
            Returns self.

        References
        ----------
        [1] De Silva and Carlsson Topological estimation using witness complexes

        """

        if not 0<=epsilon<=1:
            raise ValueError('Epsilon needs to be within 0-1')
        self.dimension = dimension
        n = len(matrix)
        lambda_list = []

        # alpha_matrix = numpy.ones((n,n))*math.inf
        list_of_alpha_values = []

        # print(len(lambda_list))


        matrix = self. greedy_permutation(matrix,lambda_list) #generates lambda_list and changes ordering of vertices in matrices


        for element in range(len(lambda_list)):
            self.add_simplex([element],0)
        # print (matrix)
        # print (lambda_list)

        # Find alpha values on list based on edges
        # sequence_of_alpha_values = []
        # store_list_of_simplices = []
        # list_of_alpha_values = []
        dict_edges = {}

        for i in range(len(matrix)):
            for j in range(i+1,len(matrix)):
                if lambda_list[j]>lambda_list[i]:
                    raise ValueError('Wrong Greedy Permutation')
                if matrix[i][j]<=2*lambda_list[j]/epsilon:
                    alpha = matrix[i][j]
                elif matrix[i][j] > (lambda_list[i]+lambda_list[j])/epsilon:
                    continue
                else:
                    alpha = (matrix[i][j]-lambda_list[j]/epsilon)*2
                    if epsilon<1 and alpha*epsilon * (1 - epsilon) / 2>lambda_list[j]:
                        continue
                if alpha < math.inf:
                    self.add_simplex([i,j],alpha)
                    # print('add',i,j)
                    # alpha_matrx[i,j] = alpha
                    dict_edges[tuple([i,j])] = alpha
                    # sequence_of_alpha_values.append(alpha)
                    # if alpha not in list_of_alpha_values:
                    #     list_of_alpha_values.append(alpha)
        del matrix
        for current_dim in range(2, self.dimension + 1):
    #     # start with triangles to highest simplices
    #     # for comb in list(itertools.combinations(indices_of_landmarks, current_dim + 1)):  # iterate through all simplex in this dimension
    #         store_list_of_simplices.extend(list(itertools.combinations(range(len(matrix)), current_dim + 1)))

            for simplices in itertools.combinations(range(n), current_dim + 1):
                two_simplex = list(itertools.combinations(simplices, 2))

                max_alpha = -1
                for edg in two_simplex:

                    if tuple(edg) not in dict_edges:
                        max_alpha=-1
                        break

                    if max_alpha<dict_edges[tuple(edg)]<radius:
                        # print ('alpha',alpha)
                        max_alpha = dict_edges[tuple(edg)]
                    # elif alpha_matrx[edg[0],edg[1]]==math.inf:
                    #     max_alpha=-1
                    #     break

                if max_alpha>=0:
                    # print ('add',simplices,max_alpha)
                    self.add_simplex(simplices,max_alpha)
                    # dict_edges[tuple(simplices)] = max_alpha

            print ('current_dim',current_dim)


                    #
        self.simplices_indices = numpy.argsort(self.appears_at, kind="mergesort")
        # print (alpha_matrx)

        count = [0,0,0,0]
        for sim in self.simplices_indices:
            if len(self.simplices[sim])==1:
                count[0] = count[0]+1
            elif len(self.simplices[sim]) == 2:
                count[1] = count[1] + 1
            elif len(self.simplices[sim]) == 3:
                count[2] = count[2] + 1
            else:
                count[3] = count[3]+1
        # print ('len of simplex:', max(self.simplices_indices),count)



