import numpy
import math

def three_and_one():
    return numpy.array([[0, 1, 1, 5],
                        [1, 0, 1, 5],
                        [1, 1, 0, 5],
                        [5, 5, 5, 0]])


def triangle_metric():
    return numpy.array([[0, 1, 1],
                        [1, 0, 1],
                        [1, 1, 0]])


def point_cloud_to_metric(cloud):
    result = numpy.array([[math.sqrt(sum([math.pow(w - z, 2) for w, z in zip(x, y)]))
                           for x in cloud] for y in cloud])
    return result


def circle_metric(res=4., radius=1.):
    cloud = [[radius * math.cos(float(i) * 2 * 3.14 / res),
              radius * math.sin(float(i) * 2 * 3.14 / res)] for i in range(int(res))]
    return point_cloud_to_metric(cloud)


def circle_cloud(res=4., radius=1.):
    cloud = [[radius * math.cos(float(i) * 2 * 3.14 / res),
              radius * math.sin(float(i) * 2 * 3.14 / res)] for i in range(int(res))]
    return numpy.array(cloud)


def sphere_cloud(res=50, radius=1., seed = None):
    if seed is not None:
        numpy.random.seed(seed)
    cloud = numpy.random.randn(res,3)  
    cloud /= numpy.linalg.norm(cloud, axis=1,keepdims=True)
    return cloud

def rectangle_boundary_cloud(lx,ux,ly,uy,res=10):
    pass

def rectangle_cloud(lx,ux,ly,uy,res=10,seed=None): 
    if seed is not None:
        numpy.random.seed(seed)
    cloud = numpy.array([(numpy.random.uniform(lx, ux),
                                numpy.random.uniform(ly, uy)) for _ in range(res)])
    return cloud
   
def sphere_metric(res=50, radius=1., seed = None):
    cloud = sphere_cloud(res=res, radius=radius, seed = seed)
    return point_cloud_to_metric(cloud)

def figure_eight_cloud(res=4., radius=1.):
    cloud = [[radius * math.cos(float(i) * 2 * 3.14 / res),
              radius * math.sin(float(i) * 2 * 3.14 / res)] for i in range(int(res))]
    cloud.extend([[2 * radius + radius * math.cos(float(i) * 2 * 3.14 / res),
                   radius * math.sin(float(i) * 2 * 3.14 / res)] for i in range(int(res))])
    return numpy.array(cloud)

def figure_eight_metric(res=4., radius=1.):
    cloud = figure_eight_cloud(res,radius)
    return point_cloud_to_metric(cloud)

def add_uniform_noise(cloud,lower,upper,seed=None): 
    if seed is not None:
        numpy.random.seed(seed)
    return cloud + numpy.random.uniform(lower,upper,cloud.shape)

def four_point_square():
    return numpy.array([[0,0],[1,0],[0,1],[1,1]])

def four_point_square_cloud_metric():
    return point_cloud_to_metric(four_point_square())