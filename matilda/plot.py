import numpy
import itertools
import matplotlib
import matplotlib.pyplot as plt
import math



class Plotter(object):
    """
        Class for visualising topological outputs

        Attributes
        ----------
        
    """

    def plot_2skel(self, 
                    fsc, 
                    max_id, 
                    node_size=6, 
                    vertices_color = None, 
                    fontsize=None, 
                    offsets=None,
                    pos = None,
                    labels = None,
                    labels_ha = "center",
                    labels_va = "center",
                    t_alpha=0.5,
                    this_ax=None):
        """
        Plots the union of 2-skeletons of the subset of simplices of a filtered simplicial complex `fsc`
        whose id is less or equal `max_id` 
        """
        import networkx
        lines = []
        triangles = []
        vertices = [fsc.simplices[i][0] for i in fsc.simplices_indices if (len(fsc.simplices[i]) == 1) and (i<=max_id)]
        edges = [fsc.simplices[i] for i in fsc.simplices_indices if (len(fsc.simplices[i]) == 2) and (i<=max_id)]
        simplices_deg2 = [fsc.simplices[i] for i in fsc.simplices_indices if (len(fsc.simplices[i]) == 3) and (i<=max_id)]
        if vertices_color is None:
            vertices_colors = ["blue"]*len(vertices)
        else:
            vertices_colors = vertices_color*len(vertices)
        if offsets is None:
            offsets = {k:[0,0] for k in range(len(vertices))}
        G = networkx.Graph()
        G.add_nodes_from(vertices)
        G.add_edges_from(edges)
        if pos is None:
            pos = networkx.kamada_kawai_layout(G)
        if not isinstance(pos,dict):
            raise ValueError("Parameter `pos` is expected to be a dictionary indices->positions")
        for k, v in pos.items():
            if k in G.nodes:
                G.nodes[k]["pos"] = v
        for k in edges:
            lines.append([pos[k[0]], pos[k[1]]])
        xs = list(pos.values())
        triangles = [ [pos[x],pos[y],pos[z]] for  x,y,z in simplices_deg2]

        
        if this_ax == None:
            this_ax = plt.gca()

        lc = matplotlib.collections.LineCollection(lines, color="black",linewidths=.5)
        tc = matplotlib.collections.PolyCollection(triangles,facecolor=(.5,.5,.5,t_alpha))
        circles = matplotlib.collections.EllipseCollection(widths=[node_size]*len(xs), heights=[node_size]*len(xs), angles=[0]*len(xs),
                                                           facecolors=vertices_colors, units="points", zorder=2, offsets=xs, transOffset=this_ax.transData, alpha=.5)
        if not labels is None:
            for (k, v) in pos.items():
                plt.text(x=v[0]+offsets[k][0], y=v[1]+offsets[k][1],
                         s=labels[k], ha=labels_ha, va=labels_va, fontsize=fontsize)
        this_ax.add_collection(lc)
        this_ax.add_collection(tc)
        this_ax.add_collection(circles)
        this_ax.autoscale_view()
        return this_ax

    def plot_generators(self, 
                        phc, 
                        subplots_kw = None, 
                        pos = None, 
                        labels = None,
                        labels_ha = "center",
                        labels_va = "center",
                        fontsize=None,
                        label_offsets = None,
                        top_cycles=10,
                        networkx_layout_kw = None,
                        separate=False,
                        cycle_keys=None,
                        cycle_colors=None,
                        filename = None,
                        dpi=300):
        """
        Displays persistent homology generators associated to bars.

        Parameters
        ----------
        phc: Object
            Object of the class PersistentHomologyComputer
        subplots_kw: Dict
            Python dictionary with keyword parameters to be passed to `plt.subplots`
        pos: Dict
            Python dictionary with keys indices and values positions specified as pairs (x,y) of coordinates.
        separate: Bool
            Parameter to specify whether to plot each generator separately. Default is `False`.
        label_offsets: Dict
            Python dictionary with keys indices and values offsets to be added to node coordinates specified as .
        filename: String
            Output file for image (optional: .png,.jpeg etc)


        Returns
        -------
        self : object
            Returns self.
        """
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        import networkx
        import collections
        unfilled_markers = [m for m, func in Line2D.markers.items()
                            if func != 'nothing' and m not in Line2D.filled_markers and m!=","]
        if labels is None:
            labels = collections.defaultdict(str)
        if networkx_layout_kw is None:
            networkx_layout_kw = {"k":.3, "iterations":300, "weight":None}
        lines = []
        if not pos is None and not isinstance(pos,dict) :
            raise ValueError("Parameter `pos` is expected to be a dictionary indices->positions")
        if subplots_kw is None:
            subplots_kw = {}
        if not "dpi" in subplots_kw.keys():
            subplots_kw["dpi"] = dpi
        if separate == False:
            fig, ax = plt.subplots(**subplots_kw)
        cycles = []
        cycles_edges = []
        if cycle_keys is None:
            all_cycles = phc.persistent_cycles[1].values()
            colors = list(matplotlib.colors.TABLEAU_COLORS.values())
        else:
            all_cycles = [phc.persistent_cycles[1][k] for k in cycle_keys]
            top_cycles = len(all_cycles)   # if user specifies keys, plot all cycles
            if cycle_colors is None:
                colors = list(matplotlib.colors.TABLEAU_COLORS.values())
            else:
                colors = [cycle_colors[k] for k in cycle_keys]
        if label_offsets is None:
            label_offsets = {k[0]:numpy.array([0,0]) for k in phc.fsc.simplices if len(k)==1}
        for j,x in enumerate(all_cycles):
            cycle = []
            x_list = [phc.fsc.simplices[phc.fsc.simplices_indices[y]].copy() for y in x.keys()]
            cycles_edges.append([y for y in x.keys()])
            v = x_list[0]
            vv = v[1]
            cycle.append(vv)
            while True:
                try:
                    index, ix  = next((i,1) if (vv == y[0]) else (i,0) for i,y in enumerate(x_list) if vv in y)
                    v = x_list[index]
                    vv = v.pop(ix)
                    if vv != cycle[-1]:
                        cycle.append(vv)
                    else:
                        cycle.append(v.pop())
                except:
                    break
            cycles.append(cycle)
        
        cycles, cycles_edges = list(zip(*sorted(zip(cycles,cycles_edges), key = lambda x: len(x[0]), reverse=True)[:top_cycles]))
        if separate == False:
            G = networkx.Graph()
            prev_pos = None
            for j,x in enumerate(cycles_edges):
                G.add_weighted_edges_from([ [*phc.fsc.simplices[phc.fsc.simplices_indices[y]],
                                            phc.fsc.appears_at[phc.fsc.simplices_indices[y]]] for y in x])
                if not prev_pos is None:
                    fixed_nodes = prev_pos.keys()
                else:
                    fixed_nodes = None
                prev_pos = networkx.layout.kamada_kawai_layout(G)
                prev_pos = networkx.layout.spring_layout(G, pos=prev_pos, fixed=fixed_nodes, **networkx_layout_kw)
            if pos is None:
                if networkx_layout_kw is None:
                    networkx_layout_kw = {"k":.3, "iterations":300, "weight":None}
                pos = prev_pos
                # if offsets is None: TODO: ADD POSITION OPTION TO FUNCTION CALL FOR ALL CASES
                #     offsets = {k:[0,0] for k in pos.keys()}
        existing_labels = set([])
        for j, (cycle,cycle_edges) in enumerate(list(zip(cycles,cycles_edges))):
            if separate == True:
                pos = None #  TODO: ADD POSITION OPTION TO FUNCTION CALL FOR SEPARATE CYCLE PLOTS CASE
                # offsets = None # TODO: SAME AS ABOVE FOR OFFSETS
                fig, ax = plt.subplots(**subplots_kw)
                G = networkx.Graph()
                G.add_weighted_edges_from([ [*phc.fsc.simplices[phc.fsc.simplices_indices[y]],
                                            phc.fsc.appears_at[phc.fsc.simplices_indices[y]]] for y in cycle_edges])
                if pos is None:
                    if networkx_layout_kw is None:
                        networkx_layout_kw = {"k":.3, "iterations":300, "weight":None}
                    pos = networkx.layout.kamada_kawai_layout(G)
                    pos = networkx.layout.spring_layout(G, pos=pos, **networkx_layout_kw)
            to_plot = numpy.array([pos[i] for i in cycle])
            ax.tick_params(bottom=False, labelbottom=False,  # Generally coordinates are not meaningful
                        left=False, labelleft=False)    # but might add options in the future
            ax.plot(*to_plot.T, 
                    linestyle="-", 
                    marker=unfilled_markers[j%(len(unfilled_markers))], 
                    markersize=12,
                    markerfacecolor=colors[j%(len(colors))],
                    color=colors[j%(len(colors))])
            for i in cycle:
                if not i in existing_labels:
                    plt.text(*(numpy.array(pos[i])+label_offsets[i]).T, labels[i], ha=labels_ha, va=labels_va, fontsize=fontsize)
                    if not separate:
                        existing_labels.add(i)
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename,dpi=dpi)
        plt.close()

    def plot_persistent_diagram(self,phc,filename=None):
        """
        Parameters
        ----------
        phc: Object
            Object of the class PersistentHomologyComputer
        filename: String
            Output file for image (optional: .png,.jpeg etc)


        Returns
        -------
        self : object
            Returns self.
        """
        self.phc = phc
        plt.figure(1)
        max_death = []  # list of maximum death time in each dimension
        min_born = []  # list of minimum birth time in each dimension

        for index, dimen in enumerate(self.phc.bars):  # if -1 death index, it is infinite, collect max death time
            max_death.append(0)
            min_born.append(math.inf)
            for y in self.phc.bars[dimen]:
                # print('y:', self.hom.bars[dimen][y], 'dimen:', dimen,' max: ',max,'self: ',self.hom.bars[dimen][y][1])
                if self.phc.bars[dimen][y][1] > max_death[len(max_death) - 1]:
                    max_death[len(max_death) - 1] = self.phc.bars[dimen][y][1]
                if self.phc.bars[dimen][y][0] < min_born[len(max_death) - 1]:
                    min_born[len(max_death) - 1] = self.phc.bars[dimen][y][0]

        for index,dimen in enumerate(self.phc.bars):

            if index==0:
                ax1 = plt.subplot(len(self.phc.bars),1,index+1)

            else:
                plt.subplot(len(self.phc.bars),1,index+1, sharex = ax1)

            plt.title('H' + str(index))
            plt.plot([-1,max(max_death)+1],[-1,max(max_death)+1])
                # plt.title('H' + str(index))
            plt.xlim([min(min_born)-1, max(max_death)+1])
            plt.ylim([min(min_born) - 1, max(max_death) + 1])
            # plt.xlim([min(min_born) - 1, max(max_death) + 1])
            # plt.axis('scaled')
            plt.gca().set_aspect('equal')

            for y in self.phc.bars[dimen]:
                bars_here = self.phc.bars[dimen][y]
                # print ('bars_here=: ',bars_here,' ')
                if bars_here[1]==-1:
                    plt.plot(bars_here[0],bars_here[0]+max_death[dimen]+1,'ro')
                else:
                    plt.plot(bars_here[0], bars_here[1],'bo')

        plt.tight_layout()
        if filename==None:
            plt.show()
        else:
            plt.savefig(filename)
        return self

    def plot_barcode(self, phc, 
                    finite_bar_color = "blue", 
                    infinite_bar_color = "red", 
                    background_color = "white",
                    text_color = "black",
                    title_font_size = None,
                    font_size = None,
                    figsize = (10,10),
                    filename=None, 
                    dimension=None, 
                    **kwargs):
        """
        Parameters
        ----------
        phc: Object
            Object of the class PersistentHomologyComputer

        filename: String
            Output file for image (optional: .png,.jpeg etc)
        
        min_x: float
            Optional. If specified, limits plots to show only values above it in the x cooredinate.

        max_x: float
            Optional. If specified, limits plots to show only values below it in the x cooredinate.

        Returns
        -------
        self : object
            Returns self.
        """

        self.phc = phc
        fig = plt.figure(facecolor=background_color, figsize=figsize)
        min_xsize = 1 # If all bars are of inf length, use this size for plot
        max_death = [] #list of maximum death time in each dimension
        min_born = [] #list of minimum birth time in each dimension
        if dimension is None: # by default only show degree 1 barcode
            dimensions = [1]
        elif isinstance(dimension,int):
            dimensions= [dimension]
        else:
            dimensions = dimension
        for dimen in dimensions:      
            max_death.append(min_xsize)
            min_born.append(math.inf)
            for y in self.phc.bars[dimen]:
                if (self.phc.bars[dimen][y][1]>max_death[-1]) and not numpy.isinf(self.phc.bars[dimen][y][1]):
                    max_death[-1] = self.phc.bars[dimen][y][1]
                if self.phc.bars[dimen][y][0]<min_born[-1]:
                    min_born[-1] = self.phc.bars[dimen][y][0]
        max_finite_death = max([x for x in max_death if not numpy.isinf(x)],default=1)
        for plot_index, dimen in enumerate(dimensions):
            ax1 = None
            if dimen==0:
                ax1 = plt.subplot(len(dimensions),1,plot_index+1,facecolor=background_color)
                plt.title('H'+str(dimen),color=text_color,fontsize=title_font_size)
            else:
                ax1 = plt.subplot(len(dimensions),1,plot_index+1, sharex = ax1,facecolor=background_color)
                plt.title('H' + str(dimen),color=text_color,fontsize=title_font_size)
            if "min_x" in kwargs.keys():
                xlim_min = kwargs["min_x"]
            else:
                xlim_min = min(min_born)
            if "max_x" in kwargs.keys():
                xlim_max = kwargs["max_x"]
            else:
                xlim_max = max_finite_death
            for y_ind,v in enumerate(sorted(self.phc.bars[dimen].values(),key=lambda x:(x[0],x[1]))):
                bars_here = v
                if numpy.isinf(bars_here[1]):
                    plt.plot([bars_here[0],xlim_max],[y_ind,y_ind],color=infinite_bar_color)
                else:
                    plt.plot([bars_here[0], bars_here[1]], [y_ind, y_ind],color=finite_bar_color)
            ax1.spines["bottom"].set_color(text_color)
            ax1.spines["left"].set_color(text_color)
            ax1.spines["right"].set_color(text_color)
            ax1.spines["top"].set_color(text_color)
            ax1.set_yticks([])
            ax1.tick_params("x",colors=text_color)
            # for x in ax1.xaxis.get_ticklines():
            #     x.set_color(text_color)
            # for x in ax1.xaxis.get_ticklabels():
            #     x.set_color(text_color)
            for x in ax1.yaxis.get_ticklines():
                x.set_color(text_color)
            for x in ax1.yaxis.get_ticklabels():
                x.set_color(text_color)
            if font_size is not None:
                ax1.tick_params(axis="both", labelsize=font_size)
            plt.xlim([xlim_min, xlim_max])

        if filename:
            plt.savefig(filename, dpi=300, facecolor=fig.get_facecolor())
        return fig, ax1
    
    def my_fill_array(self,arr,pivotx1,pivotx2,pivoty1,pivoty2):
        pivotint_1 = math.ceil(pivotx1)
        pivotint_2 = math.floor(pivotx2)

        delta = (pivoty2-pivoty1)/(pivotx2-pivotx1)
        for enu in range(pivotint_1,pivotint_2+1):
            arr[enu] = pivoty1+(delta*(enu-pivotx1))
        return arr

    def landscape_heatmap_row(self, list_crit_pt_single_lambda, sq_dim, diam):
        """

        :param list_crit_pt_landscape_single_lambda: pair of (x,y) critical points
        :param max_val_filt:
        :param sq_dim:
        :return: stretch in a 1-D array of length sq_dim with interpolation
        """


        list_crit_pt_single_lambda = numpy.array(list_crit_pt_single_lambda)
        x_coords = list_crit_pt_single_lambda[:, 0]
        y_coords = list_crit_pt_single_lambda[:, 1]
        factor_to_mult = (sq_dim-1) /(diam)
        x_coords_stretch = x_coords*factor_to_mult
        new_row = numpy.ones((sq_dim))*-1

        for enu,ele in enumerate(range(len(x_coords_stretch[:-1]))):
            x1 = x_coords_stretch[ele]
            x2 = x_coords_stretch[ele+1]
            y1 = y_coords[enu]
            y2 = y_coords[enu+1]
            new_row = self.my_fill_array(new_row,x1,x2,y1,y2)

        new_row[new_row < 0] = 0
        new_row[new_row == numpy.nan] = 0
        return new_row

    def landscape_heatmap(self,dict_list_landscape,sq_dim0, sq_dim1, diam0, diam1):


        corr_mat = numpy.zeros((sq_dim0, sq_dim1))

        factor_to_mult = (sq_dim0-1)/diam0

        for keys,vals in dict_list_landscape.items():
            # numpy.set_printoptions(threshold=sys.maxsize)

            row_here = self.landscape_heatmap_row(vals,sq_dim1,diam1)
            corr_mat[int((keys*factor_to_mult)),:] = row_here
        prev = corr_mat[sq_dim0-1,:]
        for item in range(sq_dim0-2,-1,-1): #last row to first row
            row_here = corr_mat[item,:]
            if numpy.count_nonzero(row_here)==0:
                corr_mat[item,:] = prev
            else:
                prev = corr_mat[item,:]
        ss = corr_mat.shape
        count = 0
        for item1 in range(1,ss[0]):
            for item2 in range(ss[1]):
                if corr_mat[item1,item2]-corr_mat[item1-1,item2]>0.00001:
                    count+=1
                    print ('diff',corr_mat[item1,item2],corr_mat[item1-1,item2])
                    # corr_mat[item1, item2] = corr_mat[item1 - 1, item2]
        print(' count: ',count)
        return corr_mat


    def plot_landscape(self,dict_list_landscape):
        """

        Parameters
        ----------
        dict_list_landscape: Dictionary
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


        fig = plt.figure()
        ax = Axes3D(fig)
        N = len(dict_list_landscape)
        HSV_tuples = [(x * 1.0 / N, 0.5, 0.5) for x in range(N)]
        RGB_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))

        maxlist = [0, 0]
        minlist = [math.inf, math.inf]
        i = 0
        for key, value in dict_list_landscape.items():
            z = (numpy.ones(len(value)) * key).tolist()
            if max(value[:, 0]) > maxlist[0]:
                maxlist[0] = max(value[:, 0])
            if min(value[:, 0]) < minlist[0]:
                minlist[0] = min(value[:, 0])

            if max(value[:, 1]) > maxlist[1]:
                maxlist[1] = max(value[:, 1])
            if min(value[:, 1]) < minlist[1]:
                minlist[1] = min(value[:, 1])

            x = (value[:, 0]).tolist()
            y = (value[:, 1]).tolist()
            verts = [list(zip(x, z, y))]
            ax.add_collection3d(Poly3DCollection(verts, color=RGB_tuples[i]))
            i += 1

        ax.set_xlim3d(minlist[0], maxlist[0])
        ax.set_ylim3d(0, len(dict_list_landscape))
        ax.set_zlim3d(minlist[1], maxlist[1])
        plt.show()
        return self



    def plot_hierarchical_landscape(self,dict_dim_list_landscape):
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

        fig = plt.figure()
        ll = len(dict_dim_list_landscape)
        rows = math.floor(math.sqrt(ll))
        cols = math.ceil(math.sqrt(ll))
        count = 1
        if rows*cols<ll:
            rows+=1

        for dim in dict_dim_list_landscape:
            dict_list_landscape = dict_dim_list_landscape[dim]
            ax = fig.add_subplot(rows, cols, count, projection='3d')
            ax.set_title('H'+str(dim))
            count+=1


            N = len(dict_list_landscape)
            HSV_tuples = [(x * 1.0 / N, 0.5, 0.5) for x in range(N)]
            RGB_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
            maxlist = [0, 0]
            minlist = [math.inf, math.inf]
            i = 0
            for key, value in dict_list_landscape.items():
                z = (numpy.ones(len(value)) * key).tolist()
                if max(value[:, 0]) > maxlist[0]:
                    maxlist[0] = max(value[:, 0])
                if min(value[:, 0]) < minlist[0]:
                    minlist[0] = min(value[:, 0])

                if max(value[:, 1]) > maxlist[1]:
                    maxlist[1] = max(value[:, 1])
                if min(value[:, 1]) < minlist[1]:
                    minlist[1] = min(value[:, 1])
                x = (value[:, 0]).tolist()
                y = (value[:, 1]).tolist()
                verts = [list(zip(x, z, y))]
                ax.add_collection3d(Poly3DCollection(verts, color=RGB_tuples[i]))
                i += 1
            if minlist[0]==math.inf:
                ax.set_xlim3d(0, maxlist[0])
            else:
                ax.set_xlim3d(minlist[0], maxlist[0])
            ax.set_ylim3d(0, len(dict_list_landscape))
            if minlist[1] == math.inf:
                ax.set_zlim3d(0, maxlist[1])
            else:
                ax.set_zlim3d(minlist[1], maxlist[1])
        plt.show()
        return self


    def landscape_from_barcodes(self,var,k_list):
        """
        Parameters
        ----------
        var: List
            List of list suggesting intervals(barcodes)

        k_list: List of values of k needed

        Returns
        -------
        self : dictionary
            Returns dict.
        """

        dict = {}
        set_x_coords = set()
        diag_coord = []
        for generate_diag_coord in var:
            diag_coord.append([(generate_diag_coord[0]+generate_diag_coord[1])/2,(generate_diag_coord[1]-generate_diag_coord[0])/2])

        diag_coord = sorted(diag_coord, key=lambda x: x[0])

        for find_x_coord in diag_coord:
            set_x_coords.add(find_x_coord[0])
            set_x_coords.add(find_x_coord[0] + find_x_coord[1])
            set_x_coords.add(find_x_coord[0] - find_x_coord[1])

        for item in list(itertools.combinations(range(len(diag_coord)),2)):

            p0 = diag_coord[item[0]]
            p1 = diag_coord[item[1]]
            if p0[0]+p1[1] < p1[0] - p1[1]:
                continue
            set_x_coords.add((p0[0]+p0[1]+p1[0]-p1[1])/2)


        set_x_coords = list(set_x_coords)
        set_x_coords.sort()

        for k in k_list:
            for t in set_x_coords:
                m = []
                for items in var:
                    m.append(min(t-items[0],items[1]-t))
                m = numpy.array(m)
                m = m.clip(min=0)
                if len(m) < k:
                    continue
                if k not in dict:
                    dict[k] = [[t,numpy.sort(m)[-k]]]
                else:
                    dict[k].append([t,numpy.sort(m)[-k]])
            if k in dict:
                dict[k] = numpy.array(dict[k])
        print (dict)
        return dict

    def add_persistent_landscapes(self,landscape1,landscape2):
        """
        Parameters
        ----------

        landscape1: Dictionary 1
            Index: Lambda-values(integers), Value: [x,y] values in the current landscape

        landscape1: Dictionary 2
            Index: Lambda-values(integers), Value: [x,y] values in the current landscape

        Returns
        -------
        landscape1: Dictionary 3
            Sum of parameters
        """
        print('Work in progress',landscape1,landscape2)
        # for keys in landscape1.keys():



    def plot_1cycle(self, cycle, K, pointcloud):
        gen_lines = []
        w = []

        for i, (key, sign) in enumerate(cycle.items()):

            if sign != 0:

                pt0 = pointcloud[K.simplices[K.simplices_indices[key]][0]]
                pt1 = pointcloud[K.simplices[K.simplices_indices[key]][1]]

                gen_lines.append([pt0, pt1])

                w.append(sign)

        return gen_lines, numpy.array(w)
