# -*- coding: utf-8 -*-
"""
Created on Thu May 08 10:19:53 2014
functiion for obtainning the first layer neighbor set

Updated on Tue Sep 02 09:29:51 2014
1. revise the judgement of first layer neighbor in function
poly_intersection(poly_a,poly_b,position_a,position_b) without
using the argument poly_b.

2. first eliminate the edges approximating the curves in poly_a, then
determine the first layer neighbor according to whether there is a egde
in poly_a perpendicular to the lind linkin the position_a and position_b.

@author: Chenlong.He
"""

# import libraries
import math, os
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon
from descartes import PolygonPatch

# defining constants
PARTICLE_NUM = 100 # the number of particles
SENSING_RANGE =  math.sqrt(1) # the sensing range of particle
CORE_RANGE = 0.01 # the radius of hard core of particle
NSTEPS = 200000 # the number of total steps of simulation
neighbor_num = [0] * PARTICLE_NUM # the number of neighbor of a particle
neighbor_set = [0] * PARTICLE_NUM # the list recording the indices of neighbor particle
upsilon = 1e-3 # tolerance for radius of particle's polar coordinate
RES = 100 # the number of vetices of polygons approximating curves
PI = math.pi
ZERO = 1e-4
ORIGIN = np.array([0,0])

# defining related functions
def is_equal(a, b):
    '''
    judging the equality of floating numbers
    :para a b: two floating number
    :rtype: boolean True or False
    '''
    return abs(a - b) < ZERO
    
def two_points_distance(p1, p2):
    '''
    calculating the distance between two points
    :para p1, p2: 2d position coordinates of two points
    :rtype: float distance
    '''
    return math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)

def two_parallel_segments_colinearity(p1,p2,p3,p4):
    if not point_on_segment(p1, p3, p4) and \
       not point_on_segment(p2, p3, p4) and \
       not point_on_segment(p3, p1, p2) and \
       not point_on_segment(p4, p1, p2):
            return False
    else:
        return True

def point_on_segment(p1, p2, p3):
    '''
    determing whether one polygon's vertex is
    on the edge of the other polygon
    :para p1, p2, p3: p1 is the point of one polygon
                      p2 and p3 are points of another
                      polygon, two of which form an edge
    :rtype: boolean True or False
    '''
    
    if is_equal(two_points_distance(p1,p2),0) or is_equal(two_points_distance(p1,p3),0): # whether p1 coincides with p2 or p3
        return True
    # p1 is in the region consisting of p2 and p3
    elif (min(p2[0],p3[0])<p1[0] and p1[0]<max(p2[0],p3[0])) and (min(p2[1],p3[1])<p1[1] and p1[1]<max(p2[1],p3[1])):
        s = 0.5*((p1[0]-p3[0])*(p2[1]-p3[1])-(p2[0]-p3[0])*(p1[1]-p3[1])) # area of the triangle consisting of three points
        if is_equal(s,0):
            return True
        else:
            return False
    else:
        return False

def circle_intersection(c_1, c_2, r):
    '''
    # algorithms for obtainning intersections of two circles
    algorithm refers to http://www.ambrsoft.com/TrigoCalc/Circles2/Circle2.htm
    :para c_1, c_2, r: centers of two circles
                       the same radius of two circle
    :rtype: lists for coordinates of two intersecting points
    '''
    a = c_1[0]
    b = c_1[1]
    c = c_2[0]
    d = c_2[1]
    d_center = math.sqrt((c-a)**2+(d-b)**2)
    if d_center > 2*r:
        intersect_1 = ORIGIN
        intersect_2 = ORIGIN
        print 'two circles have no intersection.'
    else:
        alpha = 0.25*math.sqrt((d_center+2*r)*d_center**2*(-d_center+2*r))
        x_1 = (a+c)/2.0 + 2*(b-d)/d_center**2*alpha
        x_2 = (a+c)/2.0 - 2*(b-d)/d_center**2*alpha
        y_1 = (b+d)/2.0 - 2*(a-c)/d_center**2*alpha
        y_2 = (b+d)/2.0 + 2*(a-c)/d_center**2*alpha
        intersect_1 = [x_1,y_1]
        intersect_2 = [x_2,y_2]
    return intersect_1, intersect_2

def segment_intersection(p1, p2, p3, p4):
    '''
    judging if two segments are intersected
    :para: four coordinates of intersection points. first two points are the first intersection line,
           and last two is the second line
    :rtype: if two segments are intersected, then return the coordinates of the intersection point,
             otherwise, return False
    '''
#    if p1[0] == p3[0] and p1[1] == p3[1] or p1[0] == p4[0] and p1[1] == p4[1]: # judging coincidence points
##        return p1
#        return True
#    elif p2[0] == p3[0] and p2[1] == p3[1] or p2[0] == p4[0] and p2[1] == p4[1]:
##        return p2
#        return True
#    else:
    A1 = (p1[1] - p2[1])
    B1 = (p2[0] - p1[0])
    C1 = -(p1[0]*p2[1] - p2[0]*p1[1])
    A2 = (p3[1] - p4[1])
    B2 = (p4[0] - p3[0])
    C2 = -(p3[0]*p4[1] - p4[0]*p3[1])
    
    D  = A1 * B2 - B1 * A2
    Dx = C1 * B2 - B1 * C2
    Dy = A1 * C2 - C1 * A2
    if D != 0:
        x = Dx / D # solving intersection point by Crammer rule
        y = Dy / D
        # judging intersection point in the boundary region of segments 
        if (x < max(min(p1[0],p2[0]), min(p3[0],p4[0]))) or (x > min(max(p1[0],p2[0]), max(p3[0],p4[0]))) and \
           (y < max(min(p1[1],p2[1]), min(p3[1],p4[1]))) or (y > min(max(p1[1],p2[1]), max(p3[1],p4[1]))):
            return False # segments intersecting but out of bound
        else:        
#                return [x,y]
            return True
    else:
        return False

def constrained_sensing_range(positions, r, res):
    '''
    plotting the constrained sensing range of voronoi-like region
    :para positions(3*2 matrix): the first line is the coordinates of particle, and the remainning lines are the
            coordinates of two intersection points of sensing range circles between particle and its neighbors
          r: SENSING_RANGE
          res: the resolution of polygon approximating curve
    :rtype: approx_positions: N*2 array for the vertices of constrained sensing range polygon
    '''
    # calculating theta1 belonging to [0,2pi]
    if  positions[1][0]-positions[0][0] > 0 and positions[1][1]-positions[0][1] > 0: # first quadrant
        theta1 = math.atan2(positions[1][1]-positions[0][1],positions[1][0]-positions[0][0])
    elif positions[1][0]-positions[0][0] < 0 and positions[1][1]-positions[0][1] > 0: # second quadrant
        theta1 = math.atan2(positions[1][1]-positions[0][1],positions[1][0]-positions[0][0])
    elif positions[1][0]-positions[0][0] < 0 and positions[1][1]-positions[0][1] < 0: # third quadrant
        theta1 = math.atan2(positions[1][1]-positions[0][1],positions[1][0]-positions[0][0])+2*PI
    else: # forth quadrant
        theta1 = math.atan2(positions[1][1]-positions[0][1],positions[1][0]-positions[0][0])+2*PI
    # calculating theta2 belonging to [0,2pi]
    if  positions[2][0]-positions[0][0] > 0 and positions[2][1]-positions[0][1] > 0: # first quadrant
        theta2 = math.atan2(positions[2][1]-positions[0][1],positions[2][0]-positions[0][0])
    elif positions[2][0]-positions[0][0] < 0 and positions[2][1]-positions[0][1] > 0: # second quadrant
        theta2 = math.atan2(positions[2][1]-positions[0][1],positions[2][0]-positions[0][0])
    elif positions[2][0]-positions[0][0] < 0 and positions[2][1]-positions[0][1] < 0: # third quadrant
        theta2 = math.atan2(positions[2][1]-positions[0][1],positions[2][0]-positions[0][0])+2*PI
    else: # forth quadrant
        theta2 = math.atan2(positions[2][1]-positions[0][1],positions[2][0]-positions[0][0])+2*PI
    # plotting the sector with the polygon by approximation (properly setting the starting and ending angle)
    # the center angle of sector is equal to 360 - gamma, which is calculated according to cosine law
    a = math.sqrt((positions[0][0]-positions[2][0])**2+(positions[0][1]-positions[2][1])**2)
    b = math.sqrt((positions[0][0]-positions[1][0])**2+(positions[0][1]-positions[1][1])**2) 
    c = math.sqrt((positions[1][0]-positions[2][0])**2+(positions[1][1]-positions[2][1])**2)
    gamma = math.acos((a**2+b**2-c**2)/(2*a*b))
    delta_theta = 2*PI - gamma # the rotating angle of the sector
    # starting and ending points are properly set
    theta_starting = theta1
    if not is_equal(abs((math.degrees(theta_starting) + math.degrees(delta_theta)) % 360) , math.degrees(theta2)): # converting radian to degree for judging rotating angle of curve
        theta_starting = theta2
    
    # approximating sector with polygon
    del_alpha = delta_theta / float(res-1) # incremental angle
    #start_point = [SENSING_RANGE*math.cos(math.radians(theta_starting)), SENSING_RANGE*math.sin(math.radians(theta_starting))]
    #end_point = [SENSING_RANGE*math.cos(math.radians(theta_ending)), SENSING_RANGE*math.sin(math.radians(theta_ending))]
    interpolation_points = []
    for i in range(res):
        interpolation_points.append([positions[0][0] + r*math.cos(theta_starting + i*del_alpha), positions[0][1] + r*math.sin(theta_starting + i*del_alpha)])
    approx_positions = np.asarray(interpolation_points)
    return approx_positions 
            
def poly_intersection(poly_a,poly_b,position_a,position_b):
    '''
    determine whether two polygons have a common edge:
    1. finding respective segments in poly_a and poly_b which are orthogonal to the segement linking
       the pos_a and pos_b as candidates of the common edge of two polygons
    2. if any, to exclude segments which approximate the curve, by calculating the distance between
       two points(p1, p2 for poly_a and p3, p4 for poly_b) on the candidate segment to pos_a and pos_b.
       If p1 to pos_a is equal to p1 to pos_b, and the same for p2 and for the candidate segment of 
       poly_b, then p1p2 and p3p4 are two the common edges. 
    3. Then, comparing the two candidate common edges by colinearity. If two parallel segments have no
       coincident point, then these two segments are not the common edge of two polygons.  
    :para: poly_a, poly_b: N*2 array represents the polygon
           positions_a, position_b: coordinates of particles of polygon_a and polygon_b
    :rtype: boolean True or False 
    '''
    
    com_seg_a1 = ORIGIN
    com_seg_a2 = ORIGIN
    com_seg_b1 = ORIGIN
    com_seg_b2 = ORIGIN
    # finding the segment in poly_a, which is orthogonal to the segment between pos_a and pos_b
    for i in range(len(poly_a)-1): # len(ploy)-1 is used to exclude the polygon's last point
        p1 = poly_a[i]
        p2 = poly_a[(i+1)%len(poly_a-1)]
        if abs(two_points_distance(p1, p2) > 1e-10): # excluding the coincdent ppint incurred by approximating errors 
            inner_product = (position_a[0]-position_b[0])*(p1[0]-p2[0])+(position_a[1]-position_b[1])*(p1[1]-p2[1])
            # judging whether two segments are orthogonal
            if abs(inner_product) < 1e-9 : # perpendicular bisector can be taken as the candidate common edge
                d1 = two_points_distance(position_a, p1) 
                d2 = two_points_distance(position_a, p2)
                d3 = two_points_distance(position_b, p1)
                d4 = two_points_distance(position_b, p2)
                if is_equal(d1,d3) and is_equal(d2,d4): # in order to exclude segments approximating curve
                    com_seg_a1 = p1
                    com_seg_a2 = p2
                    break

    # finding the segment in poly_b, which is orthogonal to the segment between pos_a and pos_b            
    for i in range(len(poly_b)-1): # len(ploy)-1 is used to exclude the polygon's last point
        p1 = poly_b[i]
        p2 = poly_b[(i+1)%len(poly_b-1)]
        if abs(two_points_distance(p1, p2)> 1e-10): # excluding the coincdent ppint incurred by approximating errors
            inner_product = (position_a[0]-position_b[0])*(p1[0]-p2[0])+(position_a[1]-position_b[1])*(p1[1]-p2[1])
            # judging whether two segments are orthogonal
            if abs(inner_product) < 1e-9: # perpendicular bisector can be taken as the candidate common edge
                d1 = two_points_distance(position_a, p1) 
                d2 = two_points_distance(position_a, p2)
                d3 = two_points_distance(position_b, p1)
                d4 = two_points_distance(position_b, p2)
                if is_equal(d1,d3) and is_equal(d2,d4): # in order to exclude segments approximating curve
                    com_seg_b1 = p1
                    com_seg_b2 = p2
                    break
     
    # according to the relationship between two common edge candidates, judging whether two polys are intersected 
    if not is_equal(two_points_distance(com_seg_a1, com_seg_a2),0) and not is_equal(two_points_distance(com_seg_b1, com_seg_b2),0):
        if two_parallel_segments_colinearity(com_seg_a1,com_seg_a2,com_seg_b1,com_seg_b2): # two segments have at least one coincidence  
             return True
        else:
             return False
    else: 
         return False

def first_layer_neighbor(positions):
    '''
    According to the current positions, the first layer neighbor set of particles
    is obtained.
    :para: positions: N*2 array representing N particles' 2d coordinates
    :rtype: first_layer_neighbor_set: list with size N, each row records the particle i's 
            first layer neighbor set.
    '''    
    # variable for recording intermediate data
    first_layer_neighbor_set = [0] * PARTICLE_NUM
    circle_inter_points = {}
    
    # graphic output
    fig=plt.figure()
    ax=fig.add_subplot(111)
    plt.axis('scaled') # equal axis
    i = 0
    for x,y in positions:
        plt.plot(x,y, 'ob',markersize=2) # plotting particles
        plt.text(x+0.005 ,y+0.005 , str(i)) # plotting partilces indices
        i += 1
    
    # obtainning info of neighbor particles in the sensing range and intersecting points of two sensing range circles
    for i in range(PARTICLE_NUM):
        neighbor_list = []
        circle_inter_list = []
        k = 0 # recording the number of neighbor particles
        for j in range(PARTICLE_NUM):
            if j != i:
                d = math.sqrt((positions[i][0]-positions[j][0])**2+(positions[i][1]-positions[j][1])**2) # distance between i and j
                if  d <= SENSING_RANGE: # particles i's neighbors 
                    k += 1
    #                pos_x = [positions[i][0], positions[j][0]]
    #                pos_y = [positions[i][1], positions[j][1]]
    #                plt.plot(pos_x, pos_y, '--b', alpha=0.2)# plotting the linkes between neighbor particles
                    cip_a, cip_b =  circle_intersection(positions[i], positions[j], SENSING_RANGE)
#                    cip_x = [cip_a[0], cip_b[0]]
#                    cip_y = [cip_a[1], cip_b[1]]
    #                plt.plot(cip_x, cip_y) # plot line between two circle intersecting points
                    neighbor_list.append(j)
                    circle_inter_list.append(cip_a)
                    circle_inter_list.append(cip_b)
        neighbor_num[i] = k # the number of particle i
        neighbor_set[i] = np.asarray(neighbor_list) #  the neighbor particles of particle i
        circle_inter_points[i] = circle_inter_list # the intersecting points of particle i's sensing circle with that of its neighbor particles
        
    # calculating and displaying constrained sensing ranges accroding to the info of sensing neighbors
    constrained_poly = {}
    for i in range(PARTICLE_NUM): # plotting the voronoi cell for each particle
        poly_points = {}
        fcolor = np.random.rand(3,1) # setting the color for filling the vn region of particle
        c_inter = circle_inter_points.get(i) # obtainning the intersection points from dictionary
        if neighbor_num[i] == 1: # particle has only one neighbor 
            m_points = np.array([positions[i],c_inter[0], c_inter[1]]) # save multiple points: particle's position and two intersection points as array for plotting polygons with the same color
            poly_points[0] = constrained_sensing_range(m_points, SENSING_RANGE, RES) # obtainning the constrained sehsin range by a polygon
            a = Polygon(poly_points[0].tolist())
            constrained_poly[i] = a # finally obtained polygons representing constrained sensing ranges for particles
            patch = PolygonPatch(a, fc=fcolor, ec=fcolor, alpha=0.6, zorder=1)
            ax.add_patch(patch)
        else:
            j = 1
            m_points = np.array([positions[i],c_inter[0], c_inter[1]])
            poly_points[j] = constrained_sensing_range(m_points, SENSING_RANGE, RES) # obtainning multiple polygon for multiple neighbor particles
            a = Polygon(poly_points[j].tolist())
            tmp_poly = a
            while j < neighbor_num[i]: # according to the number of neighbors of particle i, plotting the multiple parts of voronoi cell
                m_points = np.array([positions[i],c_inter[2*j+0], c_inter[2*j+1]])
                poly_points[j] = constrained_sensing_range(m_points, SENSING_RANGE, RES) # obtainning multiple polygons for multiple neighbor particles
                a = Polygon(poly_points[j].tolist()) # obtainning the intersection sets of particle i's all neighbors
                b = tmp_poly.intersection(a)
                tmp_poly = b
                j += 1
            constrained_poly[i] = b # finally obtained polygons representing constrained sensing ranges for particles
            patch = PolygonPatch(b, fc=fcolor, ec=fcolor, alpha=0.6, zorder=1)
            ax.add_patch(patch)
            
    # calculating the first layer neighbor particles
    for i in range(PARTICLE_NUM):
        first_layer_neighbor_list = []
        poly_array_a = np.array(constrained_poly[i].exterior)
        for j in range(neighbor_num[i]):
            poly_array_b = np.array(constrained_poly[neighbor_set[i][j]].exterior)
            if neighbor_num[i] == 1: # the only one particle in its sensing range is the voronoi-like neighbor
                first_layer_neighbor_list.append(neighbor_set[i][j])
                pos_x = [positions[i][0], positions[neighbor_set[i][j]][0]]
                pos_y = [positions[i][1], positions[neighbor_set[i][j]][1]]
                plt.plot(pos_x, pos_y, '--b', alpha=0.2) # plotting the linkes between voronoi-like neighbor particles
            elif poly_intersection(poly_array_a,poly_array_b,positions[i],positions[neighbor_set[i][j]]): # user-defined function to judge the intersection of two polygons
                first_layer_neighbor_list.append(neighbor_set[i][j])
                pos_x = [positions[i][0], positions[neighbor_set[i][j]][0]]
                pos_y = [positions[i][1], positions[neighbor_set[i][j]][1]]
                plt.plot(pos_x, pos_y, '--b', alpha=0.2) # plotting the linkes between voronoi-like neighbor particles
        first_layer_neighbor_set[i] = np.asarray(first_layer_neighbor_list)
        
    # setting the region for displaying graph
    x_max = max(positions[:,0])
    x_min = min(positions[:,0])
    y_max = max(positions[:,1])
    y_min = min(positions[:,1])
    plt.title(str(PARTICLE_NUM) + ' particles with their cosntrained sensing range')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(x_min-1.1*SENSING_RANGE,x_max+1.1*SENSING_RANGE)
    plt.ylim(y_min-1.1*SENSING_RANGE,y_max+1.1*SENSING_RANGE)
    plt.savefig(str(PARTICLE_NUM) +' particles at '+ str(len(t)) +' steps sensing range.png')
    plt.savefig(str(PARTICLE_NUM) +' particles at '+ str(len(t)) +' steps sensing range.eps')
    
#    # converting first layer neighbor set to adjacent matrix used to judge the symmetry of matrix
#    adjacent_matrix = np.zeros((len(first_layer_neighbor_set),len(first_layer_neighbor_set)),dtype=int)
#    for i in range(len(first_layer_neighbor_set)):
#        for j in range(len(first_layer_neighbor_set[i])):
#            adjacent_matrix[i][first_layer_neighbor_set[i][j]] = 1
#    sym_judging = (adjacent_matrix.transpose() == adjacent_matrix)
#    for i in range(len(first_layer_neighbor_set)):
#        for j in range(len(first_layer_neighbor_set)):
#            if sym_judging[i][j] == False:
#                print i,j
#    print 'adjacent matrix is symmetric!'
    
    # converting first layer neighbor set to adjacent matrix used to judge the symmetry and connectivity of the matrix
    degree_matrix =  np.zeros((PARTICLE_NUM,PARTICLE_NUM),dtype=int)
    adjacent_matrix = np.zeros((PARTICLE_NUM,PARTICLE_NUM),dtype=int)
    for i in range(len(first_layer_neighbor_set)):
        for j in range(len(first_layer_neighbor_set[i])):
            adjacent_matrix[i][first_layer_neighbor_set[i][j]] = 1
            degree_matrix[i][i] += 1
    sym_judging = (adjacent_matrix.transpose() == adjacent_matrix)
    flag_sym = True
    for i in range(len(first_layer_neighbor_set)):
        for j in range(len(first_layer_neighbor_set)):
            if sym_judging[i][j] == False:
                print 'asymmetric pair:'
                print i,j
                flag_sym = False
    if flag_sym == False:
        print 'adjacent matrix is asymmetric!'
    else:
        print 'adjacent matrix is symmetric!'
    
    # judging the connectivity of the particle system according to the eigenvalue of Laplacian matrix
    laplacian_matrix = degree_matrix -adjacent_matrix   
    eigenvalue = np.linalg.eigvals(laplacian_matrix)
    count_zero_eig = 0 # recording the number of 1 in eigenvalue
    for eig in eigenvalue:
        if abs(eig) < 1e-9:
            count_zero_eig += 1
    if count_zero_eig > 1:
        print 'particle system is not connected!'
    else:
        print 'particle system is connected!'
    
    plt.show()    
    
    return first_layer_neighbor_set

# reading order parameter for obtaining the size of order parameter
filename_ord_para = 'order_parameter_'+ str(NSTEPS) + '_steps_' + str(PARTICLE_NUM) +'_particles.txt' # reading data of order parameter from saved file
if os.path.isfile(filename_ord_para):
    f_op = open(filename_ord_para, 'r')
    order_para = []
    for d in f_op:
        order_para.append(float(d))
    f_op.close()
    print 'success to load file %s.' % filename_ord_para
else:
    print 'failure to load file %s.' % filename_ord_para
    
t = range(len(order_para))

# reading pariticles' positions to plot 
filename_config = 'config_'+ str(PARTICLE_NUM)  + '_particles.txt' # file recording postions
if os.path.isfile(filename_config):
    f_config = open(filename_config, 'r')
    ini_positions = []
    for a in f_config:
        x, y = a.split()
        ini_positions.append([float(x), float(y)])
    f_config.close()
    positions = np.asarray(ini_positions)
    print 'success to load file %s.' % filename_config
    # calling functions to plot graph
    fln = first_layer_neighbor(positions)
else:
    print 'failure to load file %s.' % filename_config