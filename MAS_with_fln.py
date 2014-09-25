# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 08:13:37 2014

@author: martianface

Multi-agent systems with varying velocities (migrating from Matlab codes)

improvements:
1. running is fast without recording the distance matrix
2. simulation is executed with more agents
3. initial configuration is saved in files for later usage
4. arrows is ploted with python-based arrow functions
5. the histogram and time series of order parameter are ploted to determine the type of phase transition
*6. trying to use Monte Carlo simulation near the critical point
"""

import math, random, os, time
import numpy as np
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from descartes import PolygonPatch
import shutil # file copy

# print the starting time
print 'Starting: ' + time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
#random.seed(1166) # with the same seed, we can obtain the same random numbers
#random.seed(1101120) # with the same seed, we can obtain the same random numbers
rs = 235013 # same random seed for all the parameter ua
random.seed(rs) # with the same seed, we can obtain the same random numbers

# defining parameters and constants
N = 100 # the number of particles
SENSING_RANGE = math.sqrt(1) # the sensing range of particle
CORE_RANGE = 0.01 # the radius of hard core of particle
NSTEPS = 200000 # the number of total steps of simulation
V_MAX = 0.03 # the max moving velocity  
RES = 100 # the number of vetices of polygons approximating curves
PI = math.pi # pi value
ZERO = 1e-4 # tolarance for checking two points coincidence
V_TOL = 1e-16 # tolerance for the absolute velocity
ORIGIN = np.array([0,0]) # original point in 2D 
TEST_INTERVAL = 1 # time interval for checking the symmetry and connectivity 
SAVE_INTERVAL = 20000 # time interval for saving the intermediate results

# defining the tmp variables
neighbor_num = [0] * N # the number of neighbor of a particle
positions = np.array([[0.0,0.0]] * N) # positions of particles
u_m = 1.0 # user defined magnitude coeffient for regions of potential forces
rep_margin = 2 * CORE_RANGE + u_m * 2 * V_MAX
att_margin = SENSING_RANGE - u_m * 2 * V_MAX
u_a = 7 # magnitude of potential forces
u_b = 1e0 # magnitude of alignment forces
upsilon = 1e-3 # tolerance for radius of particle's polar coordinate]

# user-defined functions
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
    plotting the triangles plus sectors for constrained sensing range of voronoi-like region
    :para positions(n*2 matrix): the first line is the coordinates of particle, and the remainning lines are the
            coordinates of intersection points of sensing range circles between particle and its neighbors
          r: SENSING_RANGE
          res: the resolution of polygon approximating curve
    :rtype: approx_positions: N*2 array for polygon of constrained sensing range
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
    gamma = math.acos((a**2+b**2-c**2)/(2*a*b)) # cosine theorem
    delta_theta = 2*PI - gamma # the rotating angle of the sector
    # starting and ending points are properly set
    theta_starting = theta1
    if not is_equal(abs((math.degrees(theta_starting) + math.degrees(delta_theta)) % 360) , math.degrees(theta2)): # converting radian to degree for judging rotating angle of curve
        theta_starting = theta2
    
    # approximating sector with polygon
    del_alpha = delta_theta / float(res-1) # anggular increment
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
    Theoretically, particle lying in position_a only need the information of poly_a, position_a and 
    position_b to judge the first layer neigbhoring relationship with particle lying in position_b. 
    But in order to avoid the numerical inaccuracy, here we use the information of poly_b
    1. finding respective segments in poly_a and poly_b which are orthogonal to the segement linking
       the pos_a and pos_b as candidates of the common edge of two polygons
    2. if any, to exclude segments which approximate the curve, by calculating the distance between
       two points(p1, p2 for poly_a and p3, p4 for poly_b) on the candidate segment to pos_a and pos_b.
       If p1 to pos_a is equal to p1 to pos_b, and the same for p2 and for the candidate segment of 
       poly_b, then p1p2 and p3p4 are two the common edges. 
    3. Then, comparing the two candidate common edges by colinearity. If two parallel segments have no
       coincident point, then these two segments are not the common edge of two polygons.  
    :para: poly_a, poly_b: N*2 array represents the polygon
           positions_a, position_b: generators of polygon_a and polygon_b, respectively
    :rtype: boolean True or False 
    '''
    
    com_seg_a1 = ORIGIN
    com_seg_a2 = ORIGIN
    com_seg_b1 = ORIGIN
    com_seg_b2 = ORIGIN
    # finding the segment in poly_a, which is orthogonal to the segment between pos_a and pos_b
    for i in range(len(poly_a)-1): # len(ploy)-1 is used to exclude the polygon's last point
        p1 = poly_a[i]
        p2 = poly_a[(i+1)%len(poly_a)]
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
        p2 = poly_b[(i+1)%len(poly_b)]
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
    neighbor_set = [0] * N
    first_layer_neighbor_set = [0] * N
    circle_inter_points = {}
    
    # graphic output
#    fig=plt.figure()
#    ax=fig.add_subplot(111)
#    plt.axis('scaled') # equal axis
#    i = 0
#    for x,y in positions:
#        plt.plot(x,y, 'ob',markersize=3) # plotting particles
#        plt.text(x+0.005 ,y+0.005 , str(i) ) # plotting partilces indices
#        i += 1    
    
    # obtainning info of neighbor particles in the sensing range
    for i in range(N):
        neighbor_list = []
        circle_inter_list = []
        k = 0 # recording the number of neighbor particles
        for j in range(N):
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
    for i in range(N): # plotting the voronoi cell for each particle
        poly_points = {}
#        fcolor = np.random.rand(3,1) # setting the color for filling the vn region of particle
        c_inter = circle_inter_points.get(i) # obtainning the intersection points from dictionary
        if neighbor_num[i] == 1: # particle has only one neighbor 
            m_points = np.array([positions[i],c_inter[0], c_inter[1]]) # save multiple points: particle's position and intersection points as array for plotting polygons with the same color
            poly_points[0] = constrained_sensing_range(m_points, SENSING_RANGE, RES) # obtainning the constrained sehsin range by a polygon
            a = Polygon(poly_points[0].tolist())
            constrained_poly[i] = a # finally obtained polygons representing constrained sensing ranges for particles
#            patch = PolygonPatch(a, fc=fcolor, ec=fcolor, alpha=0.6, zorder=1)
#            ax.add_patch(patch)
        else:
            j = 1
            m_points = np.array([positions[i], c_inter[0], c_inter[1]])
            poly_points[j] = constrained_sensing_range(m_points, SENSING_RANGE, RES) # obtainning multiple polygon for multiple neighbor particles
            a = Polygon(poly_points[j].tolist())
            tmp_poly = a
            while j < neighbor_num[i]: # according to the number of neighbors of particle i, plotting the multiple parts of voronoi cell
                m_points = np.array([positions[i],c_inter[2*j+0], c_inter[2*j+1]])
                poly_points[j] = constrained_sensing_range(m_points, SENSING_RANGE, RES) # obtainning multiple polygon for multiple neighbor particles
                a = Polygon(poly_points[j].tolist()) # obtainning the intersection sets of particle i's all neighbors
                b = tmp_poly.intersection(a)
                tmp_poly = b
                j += 1
            constrained_poly[i] = b # finally obtained polygons representing constrained sensing ranges for particles
#            patch = PolygonPatch(b, fc=fcolor, ec=fcolor, alpha=0.6, zorder=1)
#            ax.add_patch(patch)
            
    # calculating the first layer neighbor particles
    for i in range(N): 
        first_layer_neighbor_list = []
        poly_array_a = np.array(constrained_poly[i].exterior)
        for j in range(neighbor_num[i]):
            poly_array_b = np.array(constrained_poly[neighbor_set[i][j]].exterior)
            if neighbor_num[i] == 1: # the only one particle in its sensing range is the voronoi-like neighbor
                first_layer_neighbor_list.append(neighbor_set[i][j])
#                pos_x = [positions[i][0], positions[neighbor_set[i][j]][0]]
#                pos_y = [positions[i][1], positions[neighbor_set[i][j]][1]]
#                plt.plot(pos_x, pos_y, '--b', alpha=0.2) # plotting the linkes between voronoi-like neighbor particles
            elif poly_intersection(poly_array_a,poly_array_b,positions[i],positions[neighbor_set[i][j]]): # user-defined function to judge the intersection of two polygons
                first_layer_neighbor_list.append(neighbor_set[i][j])
#                pos_x = [positions[i][0], positions[neighbor_set[i][j]][0]]
#                pos_y = [positions[i][1], positions[neighbor_set[i][j]][1]]
#                plt.plot(pos_x, pos_y, '--b', alpha=0.2) # plotting the linkes between voronoi-like neighbor particles
        first_layer_neighbor_set[i] = np.asarray(first_layer_neighbor_list)
        
#    # setting the region for displaying graph
#    x_max = max(positions[:,0])
#    x_min = min(positions[:,0])
#    y_max = max(positions[:,1])
#    y_min = min(positions[:,1])
#    plt.xlim(x_min-1.1*SENSING_RANGE,x_max+1.1*SENSING_RANGE)
#    plt.ylim(y_min-1.1*SENSING_RANGE,y_max+1.1*SENSING_RANGE)
#    plt.savefig(str(N) +'_particles_sensing_range.png')
#    plt.savefig(str(N) +'_particles_sensing_range.svg')
    
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
    
    return first_layer_neighbor_set

# ----------------------------initial process----------------------------------
# 1. reading from existing files containing initial positions, or else 
#    generating initial positions to form a connected graph
# 2. initial file include config_N.txt, ini_theta.txt, ini_v.txt
# 3. corresponding to the file reading process, the file writing process is
#    given after each run to save the positions and velocities into files
# -----------------------------------------------------------------------------
saving_path = './N'+ str(N) +'/ua'+ str(u_a)+'/' # the directory for saving current simualtion results
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
    print 'directory ' + saving_path + ' is created.'

filename_config = 'config_'+ str(N)  + '_particles.txt' # file recording postions
if os.path.isfile(saving_path + filename_config):
    f_config = open(saving_path + filename_config, 'r')
    ini_positions = []
    for a in f_config:
        x, y = a.split()
        ini_positions.append([float(x), float(y)])
    f_config.close()
    positions = np.asarray(ini_positions)
    print 'starting from file', filename_config
else:
    ini_positions = [[0,0]] # the first particle lies in the origin
    i = 0
    while(i<N - 1):
        index = random.randint(0,i) # % randomly choosing a neighboring particle from exisitng ones
        rand_r = random.uniform(CORE_RANGE+upsilon, SENSING_RANGE-upsilon) # particle position is given by the polar coordinate, [CORE_RANGE+upsilon, SENSING_RANGE-upsilon]
        rand_theta = random.uniform(-math.pi, math.pi) # [-pi, pi]
        tmp_x = ini_positions[index][0] + rand_r * math.cos(rand_theta) # the newly generated particles keeps connectivity with the particle(index)
        tmp_y = ini_positions[index][1] + rand_r * math.sin(rand_theta)
        
        # judging the collision between old and new particles
        for j in range(i+1):
            distance = math.sqrt((tmp_x-ini_positions[j][0])**2+(tmp_y-ini_positions[j][1])**2)
            if distance < 2 * CORE_RANGE:
                isCollision = True # flag for whether newly generated partcle collide with exsiting particle
                break
            else:
                isCollision = False
        if isCollision == False:
            ini_positions.append([tmp_x, tmp_y])
            i += 1
    positions = np.asarray(ini_positions)
    print 'configuration starting from scratch'

filename_theta = 'theta_'+ str(N) + '_particles.txt' # file recording theta
if os.path.isfile(saving_path + filename_theta):
    f_theta = open(saving_path + filename_theta, 'r')
    theta = []
    for b in f_theta:
        theta.append(float(b))
    f_theta.close()
    print 'starting from file', filename_theta
else:
    theta = []
    for m in range(N):
        theta.append(random.uniform(-math.pi, math.pi)) # initial headings in [-pi, pi]
    print 'theta starting from scratch'
    
filename_v = 'v_'+ str(N) + '_particles.txt' # file recording velocity
if os.path.isfile(saving_path + filename_v):
    f_v = open(saving_path + filename_v, 'r')
    v = []
    for c in f_v:
        v.append(float(c))
    f_v.close()
    print 'starting from file', filename_v
else:
    v = [V_MAX] * N # Initial absolute velocities are all V_MAX
    print 'v starting from scratch'
    
# reading the current running steps in the file running_steps.txt, whose content is the total steps the code running
filename_running_steps = 'running_steps.txt' # file recording the total running steps. For the first run, 0 is written into the file 
if  not os.path.isfile(saving_path + filename_running_steps):
    print filename_running_steps + ' is created for the first time.'
    f_rs = open(saving_path + filename_running_steps, 'w') # writing 0 into file
    f_rs .write(str(0))
    f_rs.close()

f_rs = open(saving_path + filename_running_steps, 'r') # reading total running steps from file
running_steps = int(f_rs.read())
f_rs.close()

# defining file names later used for file operation
filename_ord_para = 'order_parameter_'+ str(NSTEPS) + '_steps_' + str(N) +'_particles.txt'  # file recording order parameter 
filename_v_list = 'v_list_' + str(N) +'_particles.txt'  # file recording absolute velocities
filename_theta_list = 'theta_list_' + str(N) +'_particles.txt'  # file recording moving headings
filename_fln_list = 'fln_list_' + str(N) +'_particles.txt'  # file recording the number of first layer neighbors

# ----------------------- main loop for NSTEPS---------------------------------
# evolution steps
# 1. planning the absolute moving velocities according to the 
#    current configuration (relative distance)
# 2. calculating the average heading according to current neighbors' 
#    information, including current moving directions and absolute velocities 
# 3. updating the positions of particles with planed headings and absolute
#    velocities
    
order_para = [] # defining varible for calculating order parameter
fln_num = [0] * N # the number of first layer neighbor of a particle
v_list = np.zeros((NSTEPS,N)) # recording the absolute velocities of all particles for calculating the order parameter
theta_list = np.zeros((NSTEPS,N)) # recording headings of all particles for calculating the order parameter
fln_num_list = np.zeros((NSTEPS, N))
#pos_off_norm_list = np.zeros((NSTEPS,N)) # recording the absolute velocitie of all particles
for steps in range(NSTEPS):
    # obtainning the first layer neighbor set according to the current positiions of all particles
    first_layer_neighbor_set = first_layer_neighbor(positions)
    tmp_v = [0.0] * N
    tmp_theta = [0.0] * N
    positions_offset = np.array([[0.0,0.0]] * N)
#    norm_positions_offset = [0.0] * N
    alig_force = np.array([[0.0,0.0]] * N) # forces generated by the aligment rule
    rep_force = np.array([[0.0,0.0]] * N) # forces generated by neighbor particles in therepulsive region 
    att_force = np.array([[0.0,0.0]] * N) # forces generated by neighbor particles in the attractive region
    resultant_force = np.array([[0.0,0.0]] * N) # resultant forces 'initializing forces at each steps'
    norm_moving_vectors = 0.0
#    tmp_v_list = []
#    po_norm_list = []
    for i in range(N): # looping for all particles
        d_min = float('inf') # initial max and min distance between two particles
        d_max = 0.0 # for planning the absolute velocity
        k = 0 # recording the number of neighbor particles
        for j in first_layer_neighbor_set[i]: # looping for all the first layer neighbors of particle i
            d = math.sqrt((positions[i][0]-positions[j][0])**2+(positions[i][1]-positions[j][1])**2) # sistance between i and j 
            if d < 2*CORE_RANGE: # for debugging
                print 'two particles collide at ' + str(steps) + ' step!' 
                print d, i, j            
            # calculating the min and max distance between particle i and its neighbors                    
            if d < d_min:
                d_min = d
            if d > d_max:
                d_max = d
            k += 1
            alig_force[i] += [u_b * v[j] * math.cos(theta[j]), u_b * v[j] * math.sin(theta[j])] # alignment force without particle i itself
            if d < rep_margin:
                rep_force[i] += [u_a*(((d-rep_margin)**2)*(positions[i][0]-positions[j][0])/d), \
                                 u_a*(((d-rep_margin)**2)*(positions[i][1]-positions[j][1])/d)]
            if d > att_margin:
                att_force[i] += [u_a*(-((d-att_margin)**2)*(positions[i][0]-positions[j][0])/d), \
                                 u_a*(-((d-att_margin)**2)*(positions[i][1]-positions[j][1])/d)]
        fln_num[i] = int(k)
        alig_force[i] += [u_b * v[i] * math.cos(theta[i]), u_b * v[i] * math.sin(theta[i])] # alignment force including particle i itself
        resultant_force[i] = alig_force[i] + rep_force[i] + att_force[i]
        if abs((d_min-2*CORE_RANGE)/2.0) < 1e-14 or abs((d_max-SENSING_RANGE)/2.0) < 1e-14: # setting too small absolute velocity to zero
            tmp_v[i] = 0.0
        else:
            tmp_v[i] = min([V_MAX, (d_min-2*CORE_RANGE)/2.0, (SENSING_RANGE-d_max)/2.0]) # truncting only the 16 number due to the accuracy of floating number
        tmp_theta[i] = math.atan2(resultant_force[i][1], resultant_force[i][0]) # theta = atan2(y,x)
        # Here, the norm of positions_offset is not equal to the tmp_v because of the inaccuracy of float,
        # therefore, lead to disconnectivity and collisiion between two particles
#        positions_offset[i] = [tmp_v[i]*math.cos(tmp_theta[i]), tmp_v[i]*math.sin(tmp_theta[i])] # updating particles' positions
#        norm_positions_offset = math.sqrt(positions_offset[i][0]**2 + positions_offset[i][1]**2) # the norm of position_offset
#       correciton of small distance
        positions_offset[i] = [tmp_v[i]*math.cos(tmp_theta[i]), tmp_v[i]*math.sin(tmp_theta[i])] # updating particles' positions           
#       correction of vector decompositions
#        if (norm_positions_offset - tmp_v[i]) < 1e-15 :
##            tmp_vx = str(tmp_v[i]*math.cos(tmp_theta[i]))
##            tmp_vy = str(tmp_v[i]*math.sin(tmp_theta[i]))
#            tmp_vx = '%.36f' %(tmp_v[i]*math.cos(tmp_theta[i]))
#            tmp_vy = '%.36f' %(tmp_v[i]*math.sin(tmp_theta[i]))
#            positions_offset[i] = [float(tmp_vx[0:20]), float(tmp_vy[0:20])] # norm_positions_offset[i] = math.sqrt(positions_offset[i][0]**2 + positions_offset[i][1]**2) # the norm of position_offset
#        po_norm_list.append(norm_positions_offset[i])
#        # due to the accuracy of float, the norm of positions offset might be larger than the planned absolute velocity tmp_v
#        # therefore, the piece of code is used to solve this problem. 
#        # ---------------------------- 2014.07.07 -----------------------------         
#        norm_position_offset = math.sqrt(positions_offset[i][0]**2 + positions_offset[i][1]**2) # the norm of position_offset
#        while norm_position_offset > tmp_v[i] and (d_max+2*norm_position_offset > SENSING_RANGE or d_min+2*norm_position_offset < 2*CORE_RANGE):
#            diff_norm = norm_position_offset - tmp_v[i]            
#            positions_offset[i] = [(norm_position_offset-diff_norm)*math.cos(tmp_theta[i]), (norm_position_offset-diff_norm)*math.sin(tmp_theta[i])]
#            norm_position_offset = math.sqrt(positions_offset[i][0]**2 + positions_offset[i][1]**2)
#            print 'the procedure of replanning position offset starts for partcle ' + str(i) + 'at ' + str(steps) + ' step' 
#        # ---------------------------------------------------------------------
    for i in range(N): # updating the state variables according to the tmp variables
        v[i] = tmp_v[i]
        theta[i] = tmp_theta[i]
        positions[i] += positions_offset[i]
    norm_sum_moving_vectors = math.sqrt(sum(v[i]*math.cos(theta[i]) for i in range(N))**2 + sum(v[i]*math.sin(theta[i]) for i in range(N))**2)
    order_para.append(1 / float(V_MAX * N) * norm_sum_moving_vectors)# calculating order parameter
    v_list[steps] = v # recording absolute velocities of current step
    theta_list[steps] = theta # recording headings of current step
    fln_num_list[steps] = fln_num
#    pos_off_norm_list[steps] = po_norm_list

    total_steps = steps + running_steps # current running steps plus already running steps
    # checking the symmetry and connectivity of MAS at every TIME_INTERVAL step
    if steps % TEST_INTERVAL == 0:    
        # judging the symmetry of the first layer neighbor set 
        degree_matrix =  np.zeros((N,N),dtype=int)
        adjacent_matrix = np.zeros((N,N),dtype=int)
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
            print 'adjacent matrix is asymmetric at ' + str(steps) + ' step'
        
        # judging the symmetry and connectivity of the first layer neighbor set 
        laplacian_matrix = degree_matrix -adjacent_matrix   
        eigenvalue = np.linalg.eigvals(laplacian_matrix)
        count_zero_eig = 0 # recording the number of 1 in eigenvalue
        for eig in eigenvalue:
            if abs(eig) < 1e-9:
                count_zero_eig += 1
        if count_zero_eig > 1:
            print 'adjacent matrix is not connected at ' + str(steps) + ' step'
            break
        
    # saving initial state of variables and random seed intp files
    if total_steps == 0:         
        filename_config_si = 'config of '+ str(N) + ' particles at '+ str(total_steps) + ' steps.txt' # saving positions data into intermediate file
        f_config_si = open(saving_path + filename_config_si, 'w') # writing positions
        for a in positions:
           f_config_si.write(str(a[0]) + ' ' + str(a[1]) + '\n')
        f_config_si.close()
#        print 'saving initial configuration'
        
        filename_theta_si = 'theta of '+ str(N) + ' particles at '+ str(total_steps) + ' steps.txt' # saving theta into intermediate file
        f_theta_si = open(saving_path + filename_theta_si, 'w') # writing theta
        for b in theta:
           f_theta_si.write(str(b) + '\n')
        f_theta_si.close()
#        print 'saving initial theta'
        
        filename_v_si = 'v of '+ str(N) + ' particles at '+ str(total_steps) + ' steps.txt' # saving theta into intermediate file
        f_v_si = open(saving_path + filename_v_si, 'w') # writing v
        for c in v:
           f_v_si.write(str(c) + '\n')
        f_v_si.close()
#        print 'saving initial v'
    
        filename_ord_para_si = 'order param of '+ str(N) + ' particles at '+ str(total_steps) + ' steps.txt' # saving order parameter into intermediate file
        f_op_si = open(saving_path + filename_ord_para_si, 'a') # writing order parameter with additional mode
        for op in order_para:
            f_op_si.write(str(op) + '\n')
        f_op_si.close()
#        print 'saving initial order parameter'
        print 'saving initial states'
        
        filename_random_seed = 'rand_seed.txt'
        f_rand_seed = open(saving_path + filename_random_seed, 'w') # writing random seed
        f_rand_seed.write(str(rs))
        f_rand_seed.close()
        print 'saving random seed'
    
    # saving the intermediate results for variables
    # 1. configuration of system
    # 2. theta
    # 3. v
    # 4. order paramter
    # plotting the intermediate figures and saving them into files
    # 1. configuration of MAS
    # 2. time evolution of order parameter
    # 3. histogram of order parameter
    # 4. first layer neighbors graph
    if (total_steps+1) % SAVE_INTERVAL == 0: # saving files at each SAVE_INTERVAL
        filename_config_si = 'config of '+ str(N) + ' particles at '+ str(total_steps) + ' steps.txt' # saving positions data into intermediate file
        f_config_si = open(saving_path + filename_config_si, 'w') # writing positions
        for a in positions:
           f_config_si.write(str(a[0]) + ' ' + str(a[1]) + '\n')
        f_config_si.close()
#        print 'saving intermediate config at ' + str(total_steps) + ' steps'
        
        filename_theta_si = 'theta of '+ str(N) + ' particles at '+ str(total_steps) + ' steps.txt' # saving theta into intermediate file
        f_theta_si = open(saving_path + filename_theta_si, 'w') # writing theta
        for b in theta:
           f_theta_si.write(str(b) + '\n')
        f_theta_si.close()
#        print 'saving intermediate theta at ' + str(total_steps) + ' steps'
        
        filename_v_si = 'v of '+ str(N) + ' particles at '+ str(total_steps) + ' steps.txt' # saving theta into intermediate file
        f_v_si = open(saving_path + filename_v_si, 'w') # writing v
        for c in v:
           f_v_si.write(str(c) + '\n')
        f_v_si.close()
#        print 'saving intermediate v at ' + str(total_steps) + ' steps'
        
        # order parameter should be attached to the saved file at the NSTEPS steps
        # step1: copy file 'order param of 10 particles at (NSTEPS-1) steps' as order param of 10 particles at (total_steps) steps
        # step2: attach the order_para data to the file
        filename_ord_para_si = 'order param of '+ str(N) + ' particles at '+ str(total_steps) + ' steps.txt' # saving order parameter into intermediate file
        if os.path.exists(saving_path + filename_ord_para):
            shutil.copyfile(saving_path + filename_ord_para, saving_path + filename_ord_para_si)
            f_op_si = open(saving_path + filename_ord_para_si, 'a') # writing order parameter with additional mode
            for op in order_para:
                f_op_si.write(str(op) + '\n')
            f_op_si.close()
            print 'saving intermediate states at ' + str(total_steps) + ' steps'
        else:
            f_op_si = open(saving_path + filename_ord_para_si, 'a') # writing order parameter with additional mode
            for op in order_para:
                f_op_si.write(str(op) + '\n')
            f_op_si.close()
#        print 'saving intermediate order parameter at ' + str(total_steps) + ' steps'
            print 'saving intermediate states at ' + str(total_steps) + ' steps'
        
        t = range(total_steps+1)
        # figure_1 time evolution of order parameter
        if os.path.isfile(saving_path + filename_ord_para_si):
            f_op_tol = open(saving_path + filename_ord_para_si, 'r')
            order_para_tol = []
            for g in f_op_tol:
                order_para_tol.append(float(g))
            f_op_tol.close()
        else:
            print 'failure to load file %s.' % filename_ord_para
        plt.figure(num=1)
        plt.plot(t, order_para_tol, 'b-') # plotting particles in the planea = Arrow(0, 0, 5*np.sin(z), 5*np.cos(z))    
        plt.title('Time evolution of order parameter of ' + str(N) + ' at '+ str(total_steps)+ ' steps')
        plt.xlabel('t')
        plt.ylabel('order parameter')
        plt.xlim(0,len(order_para_tol)*1.1) # plotting region is scaled by parameter 1.2 for clearly displaying
        plt.ylim(0,1*1.1)
        plt.grid('on')
        plt.savefig(saving_path + str(N) + ' particles at '+ str(total_steps) + ' steps order parameter.png')
        plt.savefig(saving_path + str(N) + ' particles at '+ str(total_steps) + ' steps order parameter.eps')
        plt.close()        
        
        # figure_2 particles configuration
        plt.figure(num=2, figsize=(9,9), dpi=100)
        ax = plt.subplot(111)
        for i in range(N):
            x = positions[i][0] # parameters for plotting particles and arrows
            y = positions[i][1]
            dx = v[i] * math.cos(theta[i])
            dy = v[i] * math.sin(theta[i])
            plt.plot(x,y, 'bo', markersize=2) # plotting particles in the planea = Arrow(0, 0, 5*np.sin(z), 5*np.cos(z))
        #    arrow(x, y, 2*dx, 2*dy, fc='red', alpha=0.75, length_includes_head=True, width=0.005, head_width=0.02, \
        #            head_length=0.01)# arrows represent the velocity vectors(properly scaled)
            arrows = FancyArrowPatch(posA=(x, y), posB=(x+35*dx, y+35*dy),
                                    color = 'r',
                                    arrowstyle='-|>',
                                    mutation_scale=100**.5,
                                    connectionstyle="arc3")
        
            ax.add_patch(arrows)
        plt.title('Configurtaion of '+ str(N)+ ' particles at '+ str(total_steps)+ ' steps')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('scaled')
        plt.grid('on')
        plt.savefig(saving_path + str(N) + ' particles at '+ str(total_steps) +' steps configuration.png')
        plt.close()

        # figure_3 histogram of order parameter
        plt.figure(num=3)
        plt.hist(order_para_tol, bins=50, normed='True', label='order parameter $\phi$', color='blue')
        plt.title('Histogram of order parameter in ' + str(total_steps) + ' steps with u_a = ' + str(u_a))
        plt.xlabel('$\phi$')
        plt.ylabel('$\Psi(\phi)$')
        plt.xlim(-0.1, 1.1) # plotting region is scaled by parameter 1.2 for clearly displaying
        plt.grid('on')
        plt.savefig(saving_path + 'Histogram of order parameter with u_a = '+ str(u_a) + ' at ' + str(total_steps)+ ' steps.png')
        plt.savefig(saving_path + 'Histogram of order parameter with u_a = '+ str(u_a) + ' at ' + str(total_steps)+ ' steps.eps')
        plt.close()
        
        # figure_4 first layer neighbors of MAS
        plt.figure(num=4)
        neighbor_set = [0] * N
        first_layer_neighbor_set = [0] * N
        circle_inter_points = {}
        
        # graphic output
        ax = plt.subplot(111)
        plt.axis('scaled') # equal axis
        i = 0
        for x,y in positions:
            plt.plot(x,y, 'ob',markersize=2) # plotting particles
            plt.text(x+0.005 ,y+0.005 , str(i)) # plotting partilces indices
            i += 1
        
        # obtainning info of neighbor particles in the sensing range
        for i in range(N):
            neighbor_list = []
            circle_inter_list = []
            k = 0 # recording the number of neighbor particles
            for j in range(N):
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
        for i in range(N): # plotting the voronoi cell for each particle
            poly_points = {}
            fcolor = np.random.rand(3,1) # setting the color for filling the vn region of particle
            c_inter = circle_inter_points.get(i) # obtainning the intersection points from dictionary
            if neighbor_num[i] == 1: # particle has only one neighbor 
                m_points = np.array([positions[i],c_inter[0], c_inter[1]]) # save multiple points: particle's position and intersection points as array for plotting polygons with the same color
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
                    poly_points[j] = constrained_sensing_range(m_points, SENSING_RANGE, RES) # obtainning multiple polygon for multiple neighbor particles
                    a = Polygon(poly_points[j].tolist()) # obtainning the intersection sets of particle i's all neighbors
                    b = tmp_poly.intersection(a)
                    tmp_poly = b
                    j += 1
                constrained_poly[i] = b # finally obtained polygons representing constrained sensing ranges for particles
                patch = PolygonPatch(b, fc=fcolor, ec=fcolor, alpha=0.6, zorder=1)
                ax.add_patch(patch)
                
        # calculating the first layer neighbor particles
        for i in range(N):
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
        plt.title(str(N) + ' particles with their cosntrained sensing range at ' + str(total_steps) + ' steps')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xlim(x_min-1.1*SENSING_RANGE,x_max+1.1*SENSING_RANGE)
        plt.ylim(y_min-1.1*SENSING_RANGE,y_max+1.1*SENSING_RANGE)
        plt.savefig(saving_path + str(N) +' particles at '+ str(total_steps) +' steps sensing range.png')
        plt.savefig(saving_path + str(N) +' particles at '+ str(total_steps) +' steps sensing range.eps')
        plt.close()
        print 'saving intermediate figures at ' + str(total_steps) + ' steps'
        
# ------------------------------file saving------------------------------------ 
# saving results of each run into files for further simulation and post
# processing 
# the content of files to be saved is:
# 1. positions of particles (the last step)
# 2. theta of particles (the last step)
# 3. absolute velocity of particles (the last step)
# 4. order parameter (list for all steps)
# 5. v_list (list for all steps)
# 6. theta_list (list for all steps)
# -----------------------------------------------------------------------------
f_config = open(saving_path + filename_config, 'w') # writing positions
for a in positions:
   f_config.write(str(a[0]) + ' ' + str(a[1]) + '\n')
f_config.close()

f_theta = open(saving_path + filename_theta, 'w') # writing theta
for b in theta:
   f_theta.write(str(b) + '\n')
f_theta.close()

f_v = open(saving_path + filename_v, 'w') # writing v
for c in v:
   f_v.write(str(c) + '\n')
f_v.close()

f_op = open(saving_path + filename_ord_para, 'a') # writing order parameter with additional mode
for op in order_para:
    f_op.write(str(op) + '\n')
f_op.close()

f_vl = open(saving_path + filename_v_list, 'a') # writing v_list with additional mode
for vl in v_list:
    for velo in vl:
        f_vl.write(str(velo) + ' ')
    f_vl.write('\n')
f_vl.close()

f_tl = open(saving_path + filename_theta_list, 'a') # writing theta_list with additional mode
for tl in theta_list:
    for tt in tl:
        f_tl.write(str(tt) + ' ')
    f_tl.write('\n')
f_tl.close()

f_fl = open(saving_path + filename_fln_list, 'a') # writing fln_num_list with additional mode
for fln in fln_num_list:
    for flnn in fln:
        f_fl.write(str(flnn) + ' ')
    f_fl.write('\n')
f_fl.close()

f_rs = open(saving_path + filename_running_steps, 'w') # writing the new total running steps into file
f_rs.write(str(running_steps + NSTEPS))
f_rs.close()

# print the ending time
print 'Ending: ' + time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())