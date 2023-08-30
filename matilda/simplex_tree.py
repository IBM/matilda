import collections
import numpy
import sympy
import itertools
# from itertools import combinations
# from matilda import helper_funcs
# from matilda import as boundary_coefs
import math
import random




class TrieNode(object):
    """
      Class for modelling simplex tree to store simplices

      Attributes
      ----------
      simplices
          List containing all simplices in the complex, regardless of filtration value.
          Simplices are themselves numpy arrays of integers.

      dimension: int
          Dimension of largest simplex. Equal to length of largest simplex minus one.
          number_of_vertices: Number of vertices present in the complex, regardless of filtration value.

      simplices_indices
          ID of the vertex corresponding to this node

      appears_at
          Filtration value for the simplex ending in this node in the complex


    """
    def __init__(self,simplices_indices,appears_at,dim):
        self.simplices_indices = simplices_indices
        self.appears_at = appears_at
        self.children = []
        # Dimension of this simplex: Init with -1
        self.dim = dim
        # How many times this character appeared in the addition process
        self.counter = 1
        self.parent = []



class SimplexTree(object):
    """
   Class for modelling simplex tree to store simplices

   Attributes
   ----------
   dimension_dict
   Circular dictionary: key-> dimension, value-> list containing link to nodes at each 'key' dimension level of the tree: not in order

   """


    def __init__(self, dim):
        self.root = TrieNode(-1, -1, -1) # Root node empty
        self.dim = dim
        self.dimension_dict = {}
        self.dimension_dict[-1] = [self.root]
        for i in range(dim+1):
            self.dimension_dict[i] = []

    def add(self, simplex, appears_at):
        node = self.root
        dim = len(simplex) - 1

        for ind, vertices in enumerate(simplex):
            # print ('Values:', node.dim,dim)
            children_s = [node.children[x].simplices_indices for x in range(len(node.children))]
            if node.dim + 1 == dim:  # will go in as children of this dimension
                if vertices in children_s or ind != dim:
                    raise ValueError("Duplicate simplex")
                new_node = TrieNode(vertices, appears_at, dim)
                node.children.append(new_node)
                node.children[len(node.children) - 1].parent = node
                # print('dim:',dim,'new_node:',new_node.simplices_indices,'simplex:',simplex)
                self.dimension_dict[dim].append(new_node)
                return self
            elif not vertices in children_s:
                # print ('vert:',vertices,children_s,node.simplices_indices)
                raise ValueError("Cannot insert coface before face")
            children_s = [node.children[x].simplices_indices for x in range(len(node.children))]
            node = node.children[children_s.index(vertices)]

        return self

    def locate_appears_at(self,simplex):
        """

        :param root:
        :param simplex:
        :return: Given simplex, finds its apears_at value
        """
        appears_at = -1
        node = self.root
        for vertices in simplex:
            x = 0
            print('vert:',vertices,node.simplices_indices)
            while x<len(node.children) and node.children[x].simplices_indices != vertices:
                print('vertx:', vertices, node.simplices_indices,node.children[x].simplices_indices)
                x = x+1
            if x==len(node.children):
                return -1
            appears_at = node.children[x].appears_at
            node = node.children[x]

        return appears_at


        # return self
    def print_tree(self,root_local):
        node = root_local
        for children in node.children:
            print ('parent:',node.simplices_indices,'child value:',children.simplices_indices,'child filt:',children.appears_at)
            self.print_tree(children)

        return self
            # for children in node.children:
            #     print_tree(children)

    def recursive_children(self,node):
        """

        :param node:
        :return: indices of all children of a node
        """
        all_child = []
        for child in node.children:
            simplex = [child.simplices_indices] + self.recursive_children(child)
            all_child.append(simplex)
        return all_child


    def locate_coface(self,simplex):
        """

        :param simplex: finds cofaces of the given simplex
        :return: list containing cofaces
        """

        list_j = [ self.dimension_dict[x] for x in range(len(simplex)-1,self.dim+1)]
        dict_coface = {} #Key: appears_at of simplex in that branch, value=actual_simplex
        list_coface = []
        for dim_iter,_ in enumerate(list_j):

            for item in list_j[dim_iter]: #L_j: Go up to check if simplex is satisfied, item is a self.root
                simplex_local = simplex
                if item.simplices_indices != simplex_local[len(simplex_local) - 1]:  # j does not match with last vertex of index
                    continue
                store_coface = [simplex_local[len(simplex_local)-1]]
                node_appears_at = item.appears_at
                node_traverse = item.parent
                simplex_local = simplex_local[:-1] # remove last simplex which matches
                while len(simplex_local)!=0:
                    if node_traverse.simplices_indices == -1: #We have reached the root
                        break
                    if node_traverse.simplices_indices == simplex_local[len(simplex_local)-1]: #indices match
                        simplex_local = simplex_local[:-1]
                        store_coface.append(node_traverse.simplices_indices)
                    node_traverse = node_traverse.parent
                if len(simplex_local)==0: # this co-face is contained
                    while node_traverse.simplices_indices !=-1:
                        store_coface.append(node_traverse.simplices_indices)
                        node_traverse = node_traverse.parent
                    store_coface.reverse()
                    dict_coface[node_appears_at] = [store_coface]
                    down = self.recursive_children(item)
                    copy_coface = dict_coface[node_appears_at][:]

                    for dcf in dict_coface[node_appears_at]:
                        for dd in down:
                            conc = dcf+dd
                            print (dcf, dd,dict_coface[node_appears_at],conc)
                            copy_coface.append(conc)
                    dict_coface[node_appears_at] = copy_coface

        for key in dict_coface:
            list_coface.extend(dict_coface[key])


        return  list_coface

    # def locate_face(self,simplex):
    #     """
    #
    #     :param simplex:
    #     :return: List containing appears_at of faces of given simplex
    #     """
    #     list_appears_at = []
    #     for dim_simplex in range(1,len(simplex)-1):
    #         for iter in itertools.combinations(simplex,dim_simplex): # (iter-1) dimensional faces
    #             node_array = self.root.children[root.children.index(iter[0])]
    #             while node_array






if __name__ == "__main__":
    root = SimplexTree(3)
    root.add([1],1)
    root.add([2],1)
    root.add([3], 1)
    root.add([2,3], 1)
    root.add([1, 2], 1)
    root.add([1, 2, 3], 10)
    root.add([1, 2, 3, 4], 1)
    root.add([2, 3, 4], 1)
    root.print_tree(root.root)
    # print(len(root.dimension_dict[3]))
    # print(root.dimension_dict[2][0].simplices_indices)
    # print(root.locate_appears_at([1,2,3]))
    print (root.locate_coface([2,3]))
