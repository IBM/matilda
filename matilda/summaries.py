import collections
import numpy
import itertools
import matplotlib
from matilda import homology
import matplotlib
import matplotlib.pyplot as plt
import math
import matilda.helper_funcs


class PersistenceBarcode(object):
    """
        Class for working with persistence barcodes.

        Attributes
        ----------
        bars
        Dictionary indexed by homology degree, its values are dictionaries indexed by simplex indices, and their values
        are lists of two elements corresponding to start and end of a bar.
        Concisely:  ``bars[degree][simplex_index] = [start,end]``
        
        filtration_start
        Starting value for filtration
        
        filtration_end
        Ending value for filtration
        
        
        max_death
        Ending coordinate of rightmost finite bar
    """

    def __init__(self):
        self.bars = collections.defaultdict(dict)
        self.filtration_start = None
        self.filtration_end = None
        self.max_death = None

    def plot(self,filename=None):
        """
        Parameters
        ----------
        filename: String
            Output file for image (optional: .png,.jpeg etc)


        Returns
        -------
        self : object
            Returns self.
        """
        
        plt.figure(1)
        max_death = [] #list of maximum death time in each dimension
        min_born = [] #list of minimum birth time in each dimension
        for index, dimen in enumerate(self.bars):       # if -1 death index, it is infinite, collect max death time
            max_death.append(0)
            min_born.append(math.inf)
            for y in self.bars[dimen]:
                # print('y:', self.bars[dimen][y], 'dimen:', dimen,' max: ',max,'self: ',self.bars[dimen][y][1])
                if self.bars[dimen][y][1]>max_death[len(max_death)-1] and self.bars[dimen][y][1]<math.inf:
                    max_death[len(max_death)-1] = self.bars[dimen][y][1]
                if self.bars[dimen][y][0]<min_born[len(max_death)-1]:
                    min_born[len(max_death) - 1] = self.bars[dimen][y][0]


        # print ('max_death:',max_death)
        for index,dimen in enumerate(self.bars):
            if not self.bars[dimen].values():
                continue
            if index==0:
                ax1 = plt.subplot(str(len(self.bars))+'1'+str(index+1))
                plt.title('H'+str(index))

            else:
                plt.subplot(str(len(self.bars))+'1'+str(index+1), sharex = ax1)
                plt.title('H' + str(index))
            plt.xlim([min(min_born), max(max_death)*(1.01)])
            for y_ind,y in enumerate(self.bars[dimen]):
                bars_here = self.bars[dimen][y]
                # print ('bars_here=: ',bars_here,' ')
                if bars_here[1]==math.inf:
                    plt.plot([bars_here[0],self.max_death*1.01],[y_ind,y_ind],color='red')
                else:
                    plt.plot([bars_here[0], bars_here[1]], [y_ind, y_ind],color='blue')

        plt.tight_layout()
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
        return self
    
    
class PersistenceLandscape(object):
    """
        Class for working with persistence landscapes.

        Attributes
        ----------
        lambdas
        Dictionary indexed by homology degree, its values are a persistence landscapes 
        at that degree, where each persistence landscape is stored as an increasing
        sequence of piecewise linear mappings, each stored as a list of vertices. 
        Concisely:  ``lambdas[degree][k] = [(x0,y0),...,(xn,yn)]``
        
        filtration_start
        Starting value for filtration
        
        filtration_end
        Ending value for filtration
    """

    def __init__(self):
        self.lambdas = collections.defaultdict(collections.defaultdict)
        self.filtration_start = None
        self.filtration_end = None
        
        
    def from_persistent_barcode(self,barcode):
        """
        Computes a persistence landscape from a persistence barcode
        Parameters
        ----------
        barcode: Persistent barcode
            PeristenceBarcode object
        """
        
        self.filtration_start = barcode.filtration_start
        self.filtration_end = barcode.filtration_end
        self.max_death = barcode.max_death
        
        ### IMPLEMENTATION OF ALGORITHM 1 OF BUBENIK,DLOTKO
        for degree,bars_dict in barcode.bars.items():
            A = sorted(sorted(list(bars_dict.values()),
                                 key=lambda x:x[1],reverse=True),
                            key=lambda x:x[0])
        ### SMALL MODIFICATION FOR INFINITE BARS
            A = [(b,d) if d<math.inf else (b,self.max_death*1.01) for (b,d) in A]
        ### END OF SMALL MODIFICATION
            k = 0
            while A:
                l = []
                (b,d) = A.pop(0)
                p = 0
                l.extend([(-math.inf,0),(b,0),((b+d)/2,(d-b)/2)])
                while l[-1] != (math.inf, 0):
                    if d > max([x[1] for x in A[p:]],default=-math.inf):
                        l.extend([(d,0), (math.inf,0)])
                    else:
                        to_pop = next((i for i,x in enumerate(A[p:]) if x[1]>d), 0)
                        to_pop += p
                        (bprime, dprime) = A.pop(to_pop)
                        p = to_pop
                        if bprime > d:
                            l.append((d,0))
                        if bprime >= d:
                            l.append((bprime,0))
                        else:
                            l.append(((bprime + d)/2,(d-bprime)/2))
                            A.insert(p,(bprime,d))
                            A[p:] = sorted(sorted(A[p:],key=lambda x:x[1],reverse=True),key=lambda x:x[0])
                            p+=1
                        l.append(((bprime+dprime)/2,(dprime-bprime)/2)) 
                        b,d = bprime,dprime
                    self.lambdas[degree][k] = sorted(l[:])
                k+=1
            
        return self
    
    def discrete_sampling(self,resolution_x, resolution_y,max_filtration_value = None, simplify = False):
        total_discrete_landscape = []
        for degree, bars_dict in self.lambdas.items():
            if simplify:
                threshold = 1e-3
                new_bars_dict = {k:bars_dict[k] for k in bars_dict if max(list(zip(*bars_dict[k]))[1])>threshold} 
                bars_dict = new_bars_dict
            discrete_landscape = numpy.zeros((len(bars_dict),resolution_x))
            for k, vertices in bars_dict.items():
                if not vertices:
                    continue
                xs,ys = zip(*vertices)
                if max_filtration_value is None:    ### OVERRIDE DEFAULT FILTRATION ENDING
                    upper_bound = self.max_death
                else:
                    upper_bound = max_filtration_value
                x_sampling = numpy.linspace(self.filtration_start,upper_bound,resolution_x)
                discrete_landscape[k] = matilda.helper_funcs.piecewise_linear_lambda(x_sampling,xs,ys)(x_sampling)
            total_discrete_landscape.append(discrete_landscape)
        total_result = []
        for discrete_landscape in total_discrete_landscape: ### SCALE UP TO RESOLUTION Y
            result = numpy.zeros((resolution_x,resolution_y))
            scale_y = resolution_y//discrete_landscape.shape[0]
            for i in range(scale_y * discrete_landscape.shape[0] ):
                result[:,i] = discrete_landscape[i//scale_y]
            total_result.append(result.transpose())
        return total_result




    def plot(self):
        """

        Parameters
        ----------
        dict_dim_list_landscape: Dictionary
            Index: Dimension, Value: Dictionary, which has:
            Index: Lambda-values(integers), Value: [x,y] values in the current landscape

        Returns
        -------
        self : object
            Returns self.

        """
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        import matplotlib.pyplot as plt
        import colorsys

        from matplotlib.tri import Triangulation
        
        fig = plt.figure()
        ll = len(self.lambdas.keys())
        rows = math.floor(math.sqrt(ll))
        cols = math.ceil(math.sqrt(ll))
        count = 1
        if rows*cols<ll:
            rows+=1

        for degree, bars_dict in self.lambdas.items():
            ax = fig.add_subplot(rows, cols, count, projection='3d')
            ax.set_title('H'+str(degree))
            count+=1
            xs,ys,zs = zip(*sorted([(k,max(0,y),min(z,1)) for k in bars_dict.keys() for (y,z) in bars_dict[k]  ]))
#            xs,ys =numpy.meshgrid(xs,ys)
            ax.scatter(xs, ys,zs, s=2,cmap='viridis');
        plt.show()
        return self
    
    def scalar_product(self, scalar):
        for degree, bars_dict in self.lambdas.items():
            for k in range(len(self.lambdas[degree])):
                self.lambdas[degree][k] = [(x,scalar*y) for (x,y) in self.lambdas[degree][k]]
        return self

def add_persistence_landscapes(landscape1,landscape2):
    """
    Parameters
    ----------

    landscape1
        PersistenceLandscape object

    landscape1
        PersistenceLandscape object
        
    Returns
    -------
    landscape1: PersistenceLandscape object
        Sum of input persistence landscapes
    """
    result = PersistenceLandscape()
    result.filtration_start = 0
    result.filtration_end = max(landscape1.filtration_end,landscape2.filtration_end)
    result.max_death = max(landscape1.max_death,landscape2.max_death)
    
    for degree, bars_dict in landscape1.lambdas.items():
        upper = max(len(landscape1.lambdas[degree]),len(landscape2.lambdas[degree]))
        for k in range(upper):
            try:
                xs1,ys1 = zip(*landscape1.lambdas[degree][k])
            except:
                xs1,ys1 = [],[]
            try:
                xs2,ys2 = zip(*landscape2.lambdas[degree][k])
            except:
                xs2,ys2 = [],[]
            xs = []
            xs.extend(xs1)
            xs.extend(xs2)
            xs = numpy.array(sorted(xs))
            _, unique = numpy.unique(xs.round(decimals = 3), return_index = True)
            unique = sorted(unique)
            xs = xs[unique]
            if not numpy.all(numpy.diff(xs) > 0):
                raise ValueError("Unique simplification not working",xs)
            if xs1:
                res1 = matilda.helper_funcs.piecewise_linear_lambda(xs,xs1,ys1)(xs)
            else:
                res1 = numpy.zeros((len(xs),))
            if xs2:
                res2 = matilda.helper_funcs.piecewise_linear_lambda(xs,xs2,ys2)(xs)
            else:
                res2 = numpy.zeros((len(xs),))
            ys = res1+res2
            result.lambdas[degree][k] = [(xs[i],ys[i]) for i in range(len(xs))]
    return result