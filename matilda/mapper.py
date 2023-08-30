import numpy
import itertools
import collections
from matilda import filtered_simplicial_complexes
import warnings
from matilda.helper_funcs import param_parser


node_increase_size_js_code = """
                                var j = 0;
                                for(j=0;j<sources.length;j++)
                                {
                                    var radius = sources[j].data.radius;
                                    var i = 0;
                                    for(i=0; i<radius.length; i++)
                                    {
                                        radius[i] = radius[i] * 1.25
                                    }
                                    sources[j].change.emit()  
                                }
                            """

node_decrease_size_js_code = """
                                var j = 0;
                                for(j=0;j<sources.length;j++)
                                {
                                    var radius = sources[j].data.radius;
                                    var i = 0;
                                    for(i=0; i<radius.length; i++)
                                    {
                                        radius[i] = radius[i] / 1.25
                                    }
                                    sources[j].change.emit()  
                                }
                            """
def generate_partitions(node_contents, 
                    original_element_values):
    """
    Given collections X_i of sets Y_ij and a function U Y_ij -> Z, computes for each i the
    count of elements in X_i that map to a given value in Z.
    
    
    node_contents:
        dict of lists such that `node_contents[node]` is a list of all elements contained in `node`.
    original_element_values:
        dict of values associated to each of the original set of elements. Not necessarily numeric.
    """
    import collections
    partitions = {}
    for k,v in node_contents.items():
        partitions[k] = dict(collections.Counter([original_element_values[x] for x in v]))
    return partitions
    
def pie_graph_plot(
                partitions,
                g=None,
                nodes = None,
                edges = None,
                graph_layout = None,
                radius = .008,
                node_labels = None,
                node_scaling = False,
                palette = None,
                background_fill_color = "white",
                title = None,
                match_aspect = True,
                edge_width = 1,
                edge_color = "black",
                plot_height = 600,
                plot_width = 600,
                sizing_mode = "fixed",
                show_node_labels = False,
                outline_line_width = 1,
                outline_line_color = "black"):
    """
        Produces a bokeh plot of a graph where each node is represented by a pie chart.

        Parameters
        ----------

        g:
            networkx graph

        partitions:
            dict of dicts such that partitions[node][class] == portion of class at node.

        graph_layout:
            dict of size 2 arrays such that graph_layout[node] == position of node in plane.
        
        nodes:
            list of nodes. If used, parameter `g` is ignored.
        
        edges:
            list of edges. Used in conjunction with `nodes`. If used, parameter `g` is ignored.

        radius:
            Radius of nodes. Upper bound on the radii of nodes if `node_scaling` is set to True.

        node_labels:
            list of text labels for nodes. If provided, must be in same order as iteration 
            over nodes in graph.

        node_scaling:
            If set to True, nodes are scaled from 
                `((smallest partition size)/(largest partition size))*radius` to `radius`, with
            a minimum size of nodes of radius/10 for enhanced viewability.

        palette:
            dict with keys nodes and values bokeh colors, or list with colors with same size as 
            number of classes present in `partitions`.
        
        title:
            Plot's title.
        
        edge_width:
            Graph's edge thickness. 
        
        edge_color:
            Graph's edge color. Used only if `sizing_mode` is "fixed".

        show_node_labels:
            If set to `True`, adds node labels to the plot.
        
        match_aspect:
            Preserves aspect ratio when changing size. Passed as-is to bokeh's figure initializer.

        plot_height:
            Plot's height in pixels. Used only if `sizing_mode` is "fixed". Passed as-is to bokeh's figure initializer.

        plot_width:
            Plot's width in pixels. Passed as-is to bokeh's figure initializer.

        sizing_mode:
            One of "fixed", "stretch_both", "stretch_width", "stretch_height", "scale_width", "scale_height". 
            Passed as-is to bokeh's figure initializer.

        outline_line_width:
            Thickness of the frame of the plot. Passed as-is to bokeh's figure initializer.

        outline_line_color:
            Color of the frame of the plot. Passed as-is to bokeh's figure initializer.
        
        background_fill_color:
            bokeh color to use for plot background. Passed as-is to bokeh's figure initializer.
        
        
    """
    import networkx,  numpy, pathlib, os, warnings
    import bokeh.plotting
    from bokeh.models import HoverTool, CustomAction
    from bokeh.models import ColumnDataSource, CustomJS, LabelSet
    if g is None and nodes is None and edges is None:
        raise ValueError("At least one of g, nodes, or edges must be provided.")
    if nodes is not None or edges is not None:
        g = networkx.Graph()
        if nodes is not None:
            g.add_nodes_from(nodes)
        if edges is not None:
            g.add_edges_from(edges)
    nodes_list = list(g.nodes())
    if node_labels is None:
        node_labels = list(map(str, g.nodes()))
    if not all("pos" in v.keys() for v in g.nodes.values()):
        if graph_layout is None:
            raise ValueError("Graph needs to have node position information stored in `pos` attribute "
                               " (you can specify position information with graph_layout).")
        else:
            for x in nodes_list:
                g.nodes[x]["pos"] = graph_layout[x]
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    p, node_x, node_y = init_bokeh_figure(g, 
                                            background_fill_color, 
                                            title, 
                                            match_aspect, 
                                            edge_width, 
                                            edge_color, 
                                            plot_height, 
                                            plot_width, 
                                            sizing_mode, 
                                            outline_line_width, 
                                            outline_line_color)
    nodes_enum = {n:i for i,n in enumerate(nodes_list)}
    node_sizes = numpy.array([sum([vv for kk,vv in partitions[n].items()]) for n in nodes_list ])

    if node_scaling:
        node_sizes_min = numpy.amin(node_sizes)
        node_sizes_max = numpy.amax(node_sizes)
        radii = numpy.interp(node_sizes, 
                            (node_sizes_min, node_sizes_max), 
                            (max(node_sizes_min*radius/node_sizes_max, radius/10), radius))
    else:
        radii = [radius]*len(node_sizes)
    all_node_sources = []
    factors = sorted(list(set.union(*[set(v.keys()) for  v in partitions.values()])))
    factors_enum = {f:i for i,f in enumerate(factors)}
    if palette is None:
        palette = bokeh.palettes.d3["Category10"][max(3, len(factors))]
    factor_sizes = numpy.zeros((len(factors), len(nodes_list)))
    for n in nodes_list:
        for kk, vv in partitions[n].items():
            factor_sizes[factors_enum[kk] , nodes_enum[n]] = vv
    end_angles = 2*numpy.pi * \
        numpy.cumsum(factor_sizes, axis=0)/node_sizes
    start_angles = numpy.roll(end_angles, 1, axis=0)
    start_angles[0, :] = 0
    # NOW PLOT ALL NODES AS MINI PIE CHARTS
    legend_items = []
    for l in factors:
        sourced = {}
        sourced["x"] = node_x
        sourced["y"] = node_y
        sourced["radius"] = radii
        # DEALS WITH CASE OF ONE CLASS PER NODE AND BOKEH'S BEHAVIOR
        for j in range(start_angles[factors_enum[l]].shape[0]):
            if start_angles[factors_enum[l], j] == 0 and end_angles[factors_enum[l], j] >= 2*numpy.pi - 1e-3:
                end_angles[factors_enum[l], j] = 2*numpy.pi - 1e-3
        sourced["start_angle"] = start_angles[factors_enum[l]]
        sourced["end_angle"] = end_angles[factors_enum[l]]
        sourced["factor_size"] = factor_sizes[factors_enum[l]]
        sourced["factor_label"] = [str(l) for _ in node_x]
        sourced["node_size"] = node_sizes
        sourced["node_label"] = node_labels
        if isinstance(palette,dict):
            sourced["color"] = [palette[l]]*len(node_x)
        else:
            sourced["color"] = [palette[factors_enum[l]]]*len(node_x)
        source = ColumnDataSource(sourced)
        all_node_sources.append(source)
        legend_items += [(l, p.wedge(x="x", y="y", start_angle="start_angle", end_angle="end_angle",
                                    source=source, radius="radius", fill_color="color", 
                                    alpha=.9, line_width=0,legend_label=str(l)))]
        p.arc(x="x", y="y", start_angle="start_angle", end_angle="end_angle",
                source=source, radius="radius", line_color="color", alpha=1)
        labels = LabelSet(x='x', y='y', text='node_label',
                        x_offset=5, y_offset=5, source=source, render_mode='canvas')
    tooltips = ("<div>Class <b>@factor_label</b> in node: <b>@factor_size</b> <br>"
                "Total in node: <b>@node_size</b> <br>"
                "Node Label: <b>@node_label</b>")
    hover_tool = HoverTool(tooltips=tooltips, renderers=[ x for _, x in legend_items])
    p.add_tools(hover_tool)
    if show_node_labels:
        p.add_layout(labels)
        
    params_node_size_increase_tool = dict(icon=pathlib.Path(os.path.join(dir_path,
                                                                "node_size_inc_tool_icon.png")),
                                            description = "Increase node size",
                                            callback=CustomJS(args = dict(sources=all_node_sources),
                                                                code=node_increase_size_js_code))
    try:
        node_size_increase_tool = CustomAction(**params_node_size_increase_tool)
    except ValueError as e:
        warnings.warn("\nWorkaround for bokeh<2.4 exception:\n" + str(e),RuntimeWarning)
        params_node_size_increase_tool["icon"] = os.path.join(dir_path,
                                                                "node_size_inc_tool_icon.png")
        node_size_increase_tool = CustomAction(**params_node_size_increase_tool)
    params_node_size_decrease_tool = dict(icon=pathlib.Path(os.path.join(dir_path,
                                                                "node_size_dec_tool_icon.png")),
                                            description = "Decrease node size",
                                            callback=CustomJS(args = dict(sources=all_node_sources),
                                                                code=node_decrease_size_js_code))
    try:
        node_size_decrease_tool = CustomAction(**params_node_size_decrease_tool)
    except ValueError as e:
        warnings.warn("\nWorkaround for bokeh<2.4 exception:\n" + str(e),RuntimeWarning)
        params_node_size_decrease_tool["icon"] = os.path.join(dir_path,
                                                                "node_size_dec_tool_icon.png")
        node_size_decrease_tool = CustomAction(**params_node_size_decrease_tool)
    p.add_tools(node_size_decrease_tool)
    p.add_tools(node_size_increase_tool)
    return p

def init_bokeh_figure(g, 
                        background_fill_color, 
                        title, 
                        match_aspect, 
                        edge_width, edge_color, 
                        plot_height, plot_width, 
                        sizing_mode, 
                        outline_line_width, 
                        outline_line_color):
    from bokeh.models import BoxZoomTool, WheelZoomTool
    import bokeh.plotting, networkx

    box_zoom_tool = BoxZoomTool(match_aspect=False)
    wheel_zoom_tool = WheelZoomTool(speed=.002, zoom_on_axis = False)
    p = bokeh.plotting.figure(title=title,
                            toolbar_location="right",
                            tools=[box_zoom_tool,
                                    wheel_zoom_tool,
                                    "pan", "reset", "save"],
                            active_scroll=wheel_zoom_tool,
                            match_aspect=match_aspect,
                            plot_height=plot_height,
                            plot_width=plot_width,
                            sizing_mode = sizing_mode,
                            outline_line_width = outline_line_width,
                            outline_line_color = outline_line_color,
                            background_fill_color = background_fill_color
                            )
    p.xaxis.visible = False
    p.yaxis.visible = False
    p.xgrid.visible = False
    p.ygrid.visible = False
    pos = networkx.get_node_attributes(g, 'pos')
    edge_x = []
    edge_y = []
    for edge in g.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_y.append(y0)
        edge_y.append(y1)
        # HACKY SOLUTION TO GET EDGES INSTEAD OF ONE VERY LONG LINE
        edge_x.append(float('nan'))
        edge_y.append(float('nan'))

    node_x = []
    node_y = []
    for node in g.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
    p.line(edge_x, edge_y, line_width=edge_width, color=edge_color)
    return p,node_x,node_y
    
class Mapper(object):

    """
    Class for applying the mapper procedure to data.

    Attributes
    ----------
    X: Numpy 2D Array
        Matrix containing data points of interest presented according to the `type_X` parameter.

    Y: Numpy 2D Array
        Matrix containing the images of data points of interest under the filter presented according to the `type_Y` parameter.

    clustering
        Type of clustering used in preimages of a covering set

    type_X: str, optional
        One of two options: "euclidean" and "metric".
        If "euclidean", X is expected to have R^n vectors as rows.
        If "metric", X is expected to be a square array representing a distance matrix.
        Default is "euclidean".


    type_Y: str, optional
        One of two options: "euclidean" and "metric".
        If "euclidean", Y is expected to have R^m vectors as rows.
        If "metric", Y is expected to be a square array representing a distance matrix.
        Default is "euclidean".

    mapper_complex:
        Simplicial Complex obtained as output of the mapper procedure

    cluster_indices:
        Dictionary with information about clusters. Keys are pairs (preimage_number, cluster_number),
        and values lists with points corresponding to that cluster.

    cluster_keys:
        Convenience default enumeration for clusters and vertices in complex.

    References
    ----------
    .. [1] PAPER 1
    .. [2] PAPER 2
    """

    def __init__(self, X, Y, clustering, type_X="euclidean", type_Y="euclidean"):
        self.X = X
        self.Y = Y
        self.clustering = clustering
        self.type_X = type_X
        self.type_Y = type_Y
        self.mapper_complex = filtered_simplicial_complexes.FilteredSimplicialComplex()

    # RIGHT NOW BALL MAPPER ASSUMES X VECTOR SPACE (OR METRIC IF CLU ADMITS PRECOMPUTED)
    def ball_mapper(self, epsilon, **kwargs):
        """
        Computes the ball mapper construction: the covering used is a collection comprised of 
        the clusters of data in original space enlarged by a given radius. 

        Parameters
        ----------
        epsilon
            Radius of cluster enlargement

        min_samples
            Number of minimal samples in each of the preimages

        """
        # AND Y EITHER VECTOR OR METRIC
        self.mapper_complex = filtered_simplicial_complexes.FilteredSimplicialComplex()

        # First compute a cover for the image under the filter
        print("Computing clusters...")
        Z = self.clustering.fit_predict(self.X).astype("int")
        idx_clusters = [numpy.where(Z == i)[0] for i in numpy.unique(Z)]
        if "min_samples" in kwargs.keys():
            min_samples = kwargs["min_samples"]
        else:
            min_samples = 1
        idx_clusters = [x for x in idx_clusters if len(x) >= min_samples]
        if self.type_Y == "euclidean":
            print("Precomputing distances...")
            import scipy.spatial.distance
            distances = scipy.spatial.distance.squareform(
                scipy.spatial.distance.pdist(self.Y))
        if self.type_Y == "metric":
            distances = self.Y

        # Returns indices of elements of Y closer than epsilon to Y[idxs]
        def epsilon_close(epsilon, idxs):
            idxs_complement = [i for i in range(
                self.Y.shape[0]) if i not in idxs]
            return numpy.r_[idxs, numpy.array([i for i in idxs_complement if numpy.amin(distances[i, :][idxs]) <= epsilon])].astype("int")
        print("Enlarging clusters...")
        idx_cover = [epsilon_close(epsilon, x) for x in idx_clusters]
        # Now compute the preimages (as indices), filtering out small preimages
        print("Computing preimages...")
        preimages = [x for x in idx_cover if len(x) > min_samples]
        preimage_cluster_labels = [numpy.array(
            [0 for _ in p]) for p in preimages]
        # Gather info for complex construction
        cluster_indices = collections.defaultdict(list)
        for i, (p, c) in enumerate(zip(preimages, preimage_cluster_labels)):
            for vector_index, label in zip(p, c):
                # Add index of preimage to make clusters unique across preimages
                cluster_indices[(i, label)].append(vector_index)
        # Construct complex from information
        self.mapper_complex.dimension = 1
        # Default enumeration for clusters and vertices in complex
        cluster_keys = list(cluster_indices.keys())
        print("Constructing complex...")
        # Vertices in complex follow the default enumeration
        for i, (k, v) in enumerate(cluster_keys):
            self.mapper_complex.add_simplex([i], 0.)
        for (k1, k2) in itertools.combinations(range(len(cluster_keys)), 2):
            if len(set.intersection(set(cluster_indices[cluster_keys[k1]]), set(cluster_indices[cluster_keys[k2]]))) > 0:
                self.mapper_complex.add_simplex([k1, k2], 0.)

        self.cover = idx_cover
        self.preimages = preimages
        self.preimage_cluster_labels = preimage_cluster_labels
        self.cluster_indices = cluster_indices
        self.cluster_keys = cluster_keys
        print("Done.")
        return self

    def rectangle_mapper(self, number_of_intervals, overlap, min_samples_preimage = None, cover = None):
        """
        Computes the usual mapper construction: the covering used is a collection comprised of 
        preimages of overlapping rectangles covering the image of the projection. 

        Parameters
        ----------
        number_of_intervals
            Number of division along axis, thus if d is the dimension of the image of projection, the number 
            of rectangles will be number_of_intervals**d.

        overlap
            Portion of axis interval overlap.
        
        min_samples_preimage
            Number of minimal samples in each of the preimages.

        cover
            Custom covering for Y. Useful for comparing with other mapper implementations.

        """
        # First compute a cover for the image under the filter
        self.mapper_complex = filtered_simplicial_complexes.FilteredSimplicialComplex()
        if cover is None:
            print("Computing cover...")
            cover = CoverRn.from_data(self.Y, number_of_intervals, overlap)
        # Now compute the preimages (as indices)
        print("Computing preimages...")
        preimages = [r.intersection(self.Y) for r in cover.sets]
        # Filter out empty preimages and elements of the cover that induced them and
        # also taking into account if the clustering procedure imposes restrictions on this
        if hasattr(self.clustering, "min_samples"):
            if min_samples_preimage is not None and min_samples_preimage < self.clustering.min_samples:
                raise ValueError("Number of elements in preimage must be at least the number of samples for clustering.")
            elif min_samples_preimage is None:
                min_samples_preimage = self.clustering.min_samples
        else:
            if min_samples_preimage is None:
                min_samples_preimage = 1
        cover, preimages = tuple(
            zip(*[(u, p) for (u, p) in zip(cover.sets, preimages) if len(p) >= min_samples_preimage]))
        # Apply clustering to each of the preimages
        print("Clustering...")
        if self.type_X == "euclidean":
            preimage_cluster_labels = [
                self.clustering.fit_predict(self.X[p, :]) for p in preimages]
        elif self.type_X == "metric":
            preimage_cluster_labels = [self.clustering.fit_predict(
                (self.X[p, :])[:, p]) for p in preimages]
        # Gather info for complex construction
        cluster_indices = collections.defaultdict(list)
        for i, (p, c) in enumerate(zip(preimages, preimage_cluster_labels)):
            for vector_index, label in zip(p, c):
                if label != -1:   # Discard outliers. For now we discard them. Might change in the future to be a toggle of some kind
                    # Add index of preimage to make clusters unique across preimages
                    cluster_indices[(i, label)].append(vector_index)
        # Construct complex from information
        self.mapper_complex.dimension = 1
        # Default enumeration for clusters and vertices in complex
        cluster_keys = list(cluster_indices.keys())
        # Vertices in complex follow the default enumeration
        for i, (k, v) in enumerate(cluster_keys):
            self.mapper_complex.add_simplex([i], 0.)
        for (k1, k2) in itertools.combinations(range(len(cluster_keys)), 2):
            if len(set.intersection(set(cluster_indices[cluster_keys[k1]]), set(cluster_indices[cluster_keys[k2]]))) > 0:
                self.mapper_complex.add_simplex([k1, k2], 0.)

        self.cover = cover
        self.preimages = preimages
        self.preimage_cluster_labels = preimage_cluster_labels
        self.cluster_indices = cluster_indices
        self.cluster_keys = cluster_keys

        return self

    def print_clusters(self):
        """
        Prints clusters and their contents (indices of elements in the original dataset) in human-readable form.
        """
        for k in self.cluster_keys:
            msg = "Cluster {} of preimage {} contains points:".format(
                k[1], k[0])
            print(msg, self.cluster_indices[k])

    def plot(self, **kwargs):
        """
        Produces a static plot of graph obtained through the mapper procedure.

        Parameters
        ----------
        colors: list, required
            List containing a color values for each point in the original dataset. Each value is assumed to be a real number.

        figsize: 2-tuple of floats, optional
            Pair of numbers specifying the size of the plot, as interpreted by matplotlib's subplots function.


        """
        import networkx
        import matplotlib
        import matplotlib.pyplot as plt
        colors = kwargs["colors"]
        if "figsize" in kwargs.keys():
            figsize = kwargs["figsize"]
        else:
            figsize = None
        vertices = [s[0] for s in self.mapper_complex.simplices if len(s) == 1]
        edges = [s for s in self.mapper_complex.simplices if len(s) == 2]
        G = self.visualization_graph()
        pos = networkx.get_node_attributes(G, 'pos')
        xs = [pos[i] for i in range(len(pos.items()))]
        fig, ax = plt.subplots(figsize=figsize)
        ax.tick_params(bottom=False, labelbottom=False,  # Generally coordinates are not meaningful
                       left=False, labelleft=False)    # but might add options in the future
        lc = matplotlib.collections.LineCollection([[pos[i], pos[j]] for [i, j] in edges],
                                                   color="white", linewidths=[.1]*len(edges), zorder=1, alpha=.5)
        cm = plt.get_cmap("viridis")
        ci = self.cluster_indices
        ck = self.cluster_keys
        vertices_colors = []
        for k in ck:
            vertices_colors += [numpy.mean([colors[x]
                                            for x in ci[k]]) for k in ck]
        vertices_colors = numpy.array(vertices_colors)
        vertices_colors -= numpy.amin(vertices_colors)
        vertices_colors /= numpy.amax(vertices_colors)
        vertices_colors = cm(vertices_colors)
        circles = matplotlib.collections.EllipseCollection(widths=[6]*len(xs), heights=[6]*len(xs), angles=[0]*len(xs),
                                                           facecolors=vertices_colors, units="points", zorder=2, offsets=xs, transOffset=ax.transData, alpha=.5)
        ax.add_collection(lc)
        ax.add_collection(circles)
        ax.autoscale_view()
        ax.set_facecolor("black")
        plt.axis("equal")
        plt.show()


    def visualization_graph(self, **kwargs):
        """
        Produces a 2D representation of the graph obtained through the mapper procedure. Uses a spring layout.

        Parameters
        ----------
        layout_kwargs:
            Dictionary containing parameters for the layout used (currently spring layout) to build the graph.
            If not used spring layout is called with parameters
            `k=None`, `seed=0`, and `iterations=50`
        Returns
        -------
        A (networkx) graph with position attributes for each node.
        """
        
        import networkx
        import rpack
        vertices = [self.cluster_keys[ s[0]] for s in self.mapper_complex.simplices if len(s) == 1]
        edges = [(self.cluster_keys[ s[0]], self.cluster_keys[ s[1]]) for s in self.mapper_complex.simplices if len(s) == 2]
        G = networkx.Graph()
        G.add_nodes_from(vertices)
        G.add_edges_from(edges)
        if "layout_kwargs" in kwargs.keys():
            layout_kwargs = kwargs["layout_kwargs"]
            if layout_kwargs is None:
                layout_kwargs = {"k": None, "seed": 0, "iterations": 50}
        else:
            layout_kwargs = {"k":  None, "seed": 0, "iterations": 50}
        pos = {}
        bboxes = []
        counts = []
        nodes_to_component = {}
        shapes = []
        for i, c in enumerate(networkx.connected_components(G)):
            init_pos = networkx.kamada_kawai_layout(G.subgraph(c))
            new_pos = networkx.spring_layout(G.subgraph(c), pos = init_pos, **layout_kwargs)  
            new_pos_values = numpy.array(list(new_pos.values()))
            bbox = [numpy.min(new_pos_values[:,0]),
                        numpy.max(new_pos_values[:,0]),
                        numpy.min(new_pos_values[:,1]),
                        numpy.max(new_pos_values[:,1])] 
            w, h  = (-bbox[0] + bbox[1],-bbox[2] +bbox[3])
            # MINIMUM SIZE OF 100 FOR VERY SMALL CONNECTED COMPONENTS
            original_points_in_component = max(100,sum( [len(self.cluster_indices[x]) for x in c] ))
            counts.append(original_points_in_component )
            bboxes.append(bbox)
            shapes.append((int(numpy.ceil(max(w,1)))*original_points_in_component,
                            int(numpy.ceil(max(h,1)))*original_points_in_component))
            pos = {**pos, **new_pos}
            for k in c:
                nodes_to_component[k] = i
        if not counts:
            raise ValueError("Graph is empty. Please check parameters.")
        max_ix = numpy.argmax(counts)
        bbox_corners = rpack.pack(shapes,max_width=1 + int(numpy.amax([x[0] for x in shapes])))
        bbox_centers = []
        for (x,y),(w,h) in zip(bbox_corners, shapes):
            bbox_centers.append((x+w/2, y+h/2))
        for k, v in pos.items():
            G.nodes[k]["pos"] = (v*counts[nodes_to_component[k]]*.8 + numpy.array(bbox_centers[nodes_to_component[k]]))/counts[max_ix]
        return G

    def plot_interactive(self, 
                        point_values, 
                        radius = .008, 
                        background_fill_color = "black", 
                        colorbar_kwargs = None, 
                        title = "Mapper output", 
                        plot_height = 600, 
                        plot_width = 600, 
                        sizing_mode = "stretch_width", 
                        match_aspect = True, 
                        colorbar_title = "", 
                        node_labels = None,
                        show_node_labels = False, 
                        edge_size = .5, 
                        palette = None, 
                        edge_color = '#808080', 
                        cluster_scaling = False, 
                        categorical = False, 
                        binning = False, 
                        binning_resolution = 10, 
                        bins = None,
                        layout_kwargs = None,
                        cluster_function = numpy.mean,
                        save_only = False,
                        output_filename = None):
        """
        Produces an interactive plot of graph obtained through the mapper procedure. Each
        node in the graph corresponds to a cluster.

        Parameters
        ----------
        point_values: list, required
            List containing values for each point in the original dataset.
            Each value is assumed to be a real number.

        categorical: boolean, optional
            If `True`, values in `point_values` are considered as categorical,
            and renders each node as a pie chart of the corresponding cluster.
            If `False`, values in `point_values` are considered as numerical.
            Default is `False`.

        binning: boolean, optional
            If `True`, uses each cluster's `point_values` to construct a histogram,
            rendering each node as a pie chart representation of the histogram.
            it as a pie chart per node.
            If `False`, interprets `point_values` consistently with other parameters.
            Default is `False`.

        labels = list, optional
            List containing text names for each element in point_values. Only used if
            `categorical` is set to `True`.

        cluster_function: bist of floats, optional
            Function to apply to point values in clusters. Must accept lists or 1D numpy arrays and output a real number.
            If not specified, defaults to average of point values per cluster.

        showscale: boolean, optional.
            Shows the colorbar if `True`. Default is True

        save_only: boolean, optional.
            If `True`, only generate the interactive plot and save it to `output_file`. If `False` generate the interactive plot and
            render it. Default is `False`.

        output_filename: str, required if `save_only` is True.
            Path to the output file to save the interactive plot.

        colorbar_kwargs: dict, optional.
            Dictionary with keywords passed as-is to bokeh's colorbar constructor

        title: str, optional.
            Plot title.

        radius: float, optional.
            Base radius for nodes in graph. Default is `.008`

        cluster_scaling: boolean, optional.
            If `True` scales nodes based on the number of elements contained in the corresponding cluster.
            Elements get scaled linearly to the interval (radius/2,radius*3).
            If `False` all nodes have the same radius.
            Default  is `False`.

        layout_kwargs:
            Dictionary containing parameters for the layout used (currently spring layout) to build the graph.
            If not used spring layout is called with parameters
            `k=.3`, `seed=0`, and `iterations=300`
        """

        import bokeh
        import bokeh.plotting
        from bokeh.models import ColorBar, ColumnDataSource, Title, CustomJS
        from bokeh.transform import linear_cmap
        from bokeh.models import HoverTool, BoxSelectTool,  CustomAction
        import os
        import pathlib
        ci = self.cluster_indices
        ck = self.cluster_keys
        if categorical:
            p = pie_graph_plot(generate_partitions(self.cluster_indices,
                                                    {i:point_values[i] for i in range(len(point_values))}),
                                                    g=self.visualization_graph(layout_kwargs = layout_kwargs),
                                                    radius=radius,
                                                    node_scaling=cluster_scaling,
                                                    palette=palette,
                                                    background_fill_color=background_fill_color,
                                                    title=title,
                                                    match_aspect=match_aspect,
                                                    sizing_mode=sizing_mode,
                                                    edge_width=edge_size,
                                                    edge_color=edge_color,
                                                    plot_height=plot_height,
                                                    plot_width=plot_width,
                                                    node_labels=node_labels,
                                                    show_node_labels=show_node_labels)
        elif binning:
            if palette is None:
                palette = bokeh.palettes.viridis
            if bins is None:
                bins = numpy.linspace(point_values.min(),point_values.max()*1.01,binning_resolution+1)
            digitized = numpy.digitize(point_values,bins)
            palette = palette(len(bins)-1)
            values = [ "{:.3f} - {:.3f}".format(bins[digitized[i]-1], bins[digitized[i]]) for i in range(len(point_values))]
            palette = {  "{:.3f} - {:.3f}".format(bins[i], bins[i+1]) :palette[i] for i in range(len(bins)-1) }
            p = pie_graph_plot(generate_partitions(self.cluster_indices,
                                                {i:values[i] for i in range(len(values))}),
                                                    g=self.visualization_graph(layout_kwargs = layout_kwargs),
                                                    radius=radius,
                                                    node_scaling=cluster_scaling,
                                                    palette=palette,
                                                    background_fill_color=background_fill_color,
                                                    title=title,
                                                    match_aspect=match_aspect,
                                                    sizing_mode=sizing_mode,
                                                    edge_width=edge_size,
                                                    edge_color=edge_color,
                                                    plot_height=plot_height,
                                                    plot_width=plot_width,
                                                    node_labels=node_labels,
                                                    show_node_labels=show_node_labels)
        else:  # DEFAULT CASE: SINGLE VALUE PER NODE, CONTINUOUS COLOR BAR
            if palette is None:
                palette = bokeh.palettes.viridis(256)
            all_node_sources = []
            p, node_x, node_y = init_bokeh_figure(
                                                    self.visualization_graph(layout_kwargs = layout_kwargs),
                                                    background_fill_color, 
                                                    title, 
                                                    match_aspect, 
                                                    edge_size, 
                                                    edge_color, 
                                                    plot_height, 
                                                    plot_width, 
                                                    sizing_mode, 
                                                    outline_line_width=1, 
                                                    outline_line_color="black")
            cluster_sizes = numpy.array([len([x for x in ci[k]]) for k in ck])
        
            if cluster_scaling:
                cluster_min = cluster_sizes.min()
                cluster_max = cluster_sizes.max()
                # SMALLEST CLUSTER WILL HAVE RADIUS AT LEAST 1/10TH OF LARGEST
                radii = numpy.interp(cluster_sizes, (cluster_min, cluster_max), 
                                    (max(cluster_min*radius/cluster_max, radius/10), radius))
            else:
                radii = numpy.full(cluster_sizes.shape, radius)
        
            sourced = {}
            sourced["x"] = node_x
            sourced["y"] = node_y
            sourced["radius"] = radii
            sourced["value"] = numpy.array(
                [cluster_function(numpy.array([point_values[x] for x in ci[k]])) for k in ck])
            sourced["cluster_size"] = cluster_sizes
            sourced["cluster_index"] = ck
            sourced["node_label"] = node_labels
            source = ColumnDataSource(sourced)
            all_node_sources.append(source)
            color_mapper = linear_cmap(field_name="value",
                                       low=numpy.amin(sourced["value"]),
                                       high=numpy.amax(sourced["value"]),
                                       palette=palette)
            labels = LabelSet(x='x', y='y', text='node_label',
                            x_offset=5, y_offset=5, source=source, render_mode='canvas')
            tooltips = ("Cluster value: <b>@value</b> <br>"
                         "Elements in cluster: <b>@cluster_size</b><br>"
                         "Cluster Index: <b>@cluster_index</b>"
                         "Node Label: <b>@node_label</b>")
            # NOW PLOT NODES
            r = p.circle(x="x", y="y", color=color_mapper, nonselection_alpha=.5,
                         source=source, radius="radius", alpha=.9, line_width=0)

            if show_node_labels:
                p.add_layout(labels)
            path = os.path.abspath(__file__)
            dir_path = os.path.dirname(path)
            source.selected.js_on_change("indices", CustomJS(args=dict(p=p), code="""
                    var inds = cb_obj.indices;
                    var kernel = IPython.notebook.kernel;
                    IPython.notebook.kernel.execute("Indices of selected elements = " + inds);
                    """))

            params_node_size_increase_tool = dict(icon=pathlib.Path(os.path.join(dir_path,
                                                                        "node_size_inc_tool_icon.png")),
                                                    description = "Increase node size",
                                                    callback=CustomJS(args = dict(sources=all_node_sources),
                                                                        code=node_increase_size_js_code))
            try:
                node_size_increase_tool = CustomAction(**params_node_size_increase_tool)
            except ValueError as e:
                warnings.warn("\nWorkaround for bokeh<2.4 exception:\n" + str(e),RuntimeWarning)
                params_node_size_increase_tool["icon"] = os.path.join(dir_path,
                                                                        "node_size_inc_tool_icon.png")
                node_size_increase_tool = CustomAction(**params_node_size_increase_tool)
            params_node_size_decrease_tool = dict(icon=pathlib.Path(os.path.join(dir_path,
                                                                        "node_size_dec_tool_icon.png")),
                                                    description = "Decrease node size",
                                                    callback=CustomJS(args = dict(sources=all_node_sources),
                                                                        code=node_decrease_size_js_code))
            try:
                node_size_decrease_tool = CustomAction(**params_node_size_decrease_tool)
            except ValueError as e:
                warnings.warn("\nWorkaround for bokeh<2.4 exception:\n" + str(e),RuntimeWarning)
                params_node_size_decrease_tool["icon"] = os.path.join(dir_path,
                                                                        "node_size_dec_tool_icon.png")
                node_size_decrease_tool = CustomAction(**params_node_size_decrease_tool)
            hover_tool = HoverTool(tooltips=tooltips, renderers=[r])
            box_select_tool = BoxSelectTool(renderers=[r])
            p.add_tools(node_size_decrease_tool)
            p.add_tools(node_size_increase_tool)
            p.add_tools(hover_tool)
            p.add_tools(box_select_tool)
            if colorbar_kwargs is None: 
                colorbar_kwargs = {}
            color_bar = ColorBar(
                color_mapper=color_mapper["transform"], **colorbar_kwargs)
            p.add_layout(color_bar, 'right')
            if colorbar_title:
                p.add_layout(Title(text=colorbar_title,
                                   align="center"), "right")
        if save_only:
            import bokeh.resources
            import os
            if output_filename is None:
                raise ValueError(
                    "Error: If saving plot a filename must be specified with the parameter \'output_filename\'. ")
            root, ext = os.path.splitext(output_filename)
            if ext == ".svg":
                p.output_backend = "svg"
                bokeh.io.export_svgs(p, filename=output_filename)
            elif ext == ".html":
                bokeh.io.save(p, filename=output_filename, resources=bokeh.resources.CDN, title = "MaTilDA output")
        else:
            bokeh.plotting.show(p)
        return p

class Rectangle(object):
    """
    Every Rectangle is determined by its center and sides' lengths divided by 2, called radius
    """

    def __init__(self, center, radius):
        self.center = numpy.asarray(center)
        self.radius = numpy.asarray(radius)
        self.dimension = self.center.shape[-1]

    def has_element(self, point):
        """
        Determines if the R^n point is in the rectangle
        """
        return numpy.all((self.center - self.radius) <= point) and numpy.all(point <= (self.center + self.radius))

    def intersection(self, X):
        """
        Determines and returns what elements of X are in rectangle
        """
        return numpy.where(numpy.all((self.center - self.radius <= X) & (X <= self.center+self.radius), axis=1))[0]

    def __str__(self):
        return "Center:" + str(self.center) + "\tRadius:" + str(self.radius) + "\n"


class CoverRn(object):
    """
    Class for covers of sets in R^n.

    Attributes
    ----------
    Sets:
        List of sets that make up the covering family
    References
    ----------
    .. [1] PAPER 1
    .. [2] PAPER 2
    """

    def __init__(self):
        """
        The cover is empty by default
        """
        self.sets = None

    @ classmethod
    def from_rectangle_list(cls, l):
        """
        Initialize the cover from list of rectangles
        """
        result = cls()
        result.sets = l
        return result

    @ classmethod
    def from_data(cls, Y, number_of_intervals, overlap):
        """
        Covers a point cloud Y by evenly distributed rectangles in Euclidean Space with a given overlap.

        Parameters
        ----------

        Y
            2D array containing the points that will be covered.

        number_of_intervals

        overlap
            Percentage of overlap between 2 consecutive intervals.
        """
        result = cls()
        Y = numpy.asarray(Y)
        upper_bounds = numpy.max(Y, axis=0)
        lower_bounds = numpy.min(Y, axis=0)
        diameters = upper_bounds - lower_bounds
        division_lengths = diameters / (number_of_intervals)
        rectangle_radius = division_lengths/ (2*(1-overlap))
        centers_per_axis = [[lb + division_lengths[axis]/2 + i*division_lengths[axis]
                             for i in range(number_of_intervals)] for axis, lb in enumerate(lower_bounds)]
        rectangles = [Rectangle(center, rectangle_radius)
                      for center in itertools.product(*centers_per_axis)]
        result.sets = rectangles
        return result

    def has_element(self, point):
        """
        Determines if the R^n point is in the cover
        """
        if not self.sets:
            return False
        result = numpy.any([x.has_element(point) for x in self.sets])
        return result

    def __str__(self):
        return "".join([str(x) for x in self.sets])


class Cover(object):
    """
    Class for an abstract representation of topological coverings
    of topological spaces (through an appropriate indexing method).


    Attributes
    ----------
    Sets:
        List of sets that make up the covering family
    References
    ----------
    .. [1] PAPER 1
    .. [2] PAPER 2
    """

    def __init__(self):
        """
        The cover is empty by default
        """
        self.sets = None

    @ classmethod
    def from_list(cls, l):
        """
        Initialize the cover from list of subsets of indices.
        """
        result = cls()
        result.sets = l
        return result
