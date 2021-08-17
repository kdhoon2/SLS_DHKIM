
import numpy as np
import matplotlib.pyplot as plt
import math
import copy
import time
import sys
import pandas as pd
import random
np.set_printoptions(threshold=sys.maxsize)
##########################################################################

def making_dict(xs,ys,outer_geo,*inner_geos):

	point_num=ys*xs.shape[0]+xs

	matrix_we_want=[outer_geo]

	for i in range(len(inner_geos)):  # mats is now TUPLE! Also, returned one will be LIST!
		matrix_we_want.append(~inner_geos[i])
        

	geo=np.all(matrix_we_want,axis=0)

	domain_dict={}

	for i in range(geo.shape[0]):
		for j in range(geo.shape[1]):
			if geo[i,j]==True:
				domain_dict[(xs[i,j],ys[i,j])]=point_num[i,j]
    

	Total_points_coord=list(domain_dict.keys())
	Total_points=list(domain_dict.values())

	total_num_points = len(Total_points)

	return geo, domain_dict, Total_points, Total_points_coord, total_num_points


def get_all_points(xs,ys,*ggs): # getting all points included in geometry (including boundary just right inside the geometry)
	
	gg_points_total=[]
	for i in range(len(ggs)):
		gg=ggs[i]
		gg_points=[]

		for i in range(gg.shape[0]):
			for j in range(gg.shape[1]):
				if gg[i,j]==True:
					gg_points.append((xs[i,j],ys[i,j]))

		gg_points_total.append(gg_points)

	return gg_points_total


def boundary_sort(all_points):   # all_points: e.g. [(0, 0), (1, 0), (2, 0), (0, 1), (2, 1), (0, 2), (1, 2), (2, 2)]


    BoundaryNodes=[]
    
    for i in range(len(all_points)):
        count=0
    
        tmpo1=(all_points[i][0]-1,all_points[i][1])
        if tmpo1 in all_points:
            count=count+1
        
        tmpo2=(all_points[i][0]+1,all_points[i][1])
        if tmpo2 in all_points:
            count=count+1
        
        tmpo3=(all_points[i][0],all_points[i][1]-1)
        if tmpo3 in all_points:
            count=count+1
        
        tmpo4=(all_points[i][0],all_points[i][1]+1)
        if tmpo4 in all_points:
            count=count+1
        
        tmpo5=(all_points[i][0]-1,all_points[i][1]+1)
        if tmpo5 in all_points:
            count=count+1
        
        tmpo6=(all_points[i][0]-1,all_points[i][1]-1)
        if tmpo6 in all_points:
            count=count+1
        
        tmpo7=(all_points[i][0]+1,all_points[i][1]+1)
        if tmpo7 in all_points:
            count=count+1
        
        tmpo8=(all_points[i][0]+1,all_points[i][1]-1)
        if tmpo8 in all_points:
            count=count+1
    
        if count!=8:
            BoundaryNodes.append(all_points[i])
    
    # sorting
    
    center_x=0
    center_y=0
    for i in range(len(BoundaryNodes)):
        center_x=center_x+BoundaryNodes[i][0]
        center_y=center_y+BoundaryNodes[i][1]
        
    center_x=center_x/len(BoundaryNodes)
    center_y=center_y/len(BoundaryNodes)
    center=(center_x,center_y)    

    angles=[]
    for i in range(len(BoundaryNodes)):
        rel_x=BoundaryNodes[i][0]-center[0]
        rel_y=BoundaryNodes[i][1]-center[1]
        angles.append(math.atan2(rel_y,rel_x))
        
    angles=np.array(angles)
    
    bd=np.array(BoundaryNodes)
    sorted_one=bd[np.argsort(angles)]
    
    sorted_bdry=[]
    for i in range(len(sorted_one)):
        sorted_bdry.append(tuple(sorted_one[i]))
    
    sorted_bdry=np.array(sorted_bdry)
    return sorted_bdry



def just_outer_nodes_sorted(inner_geometry_points):  # ONLY for inner geometries!
    # inner_geometry_points: e.g. geo2_points
    
    inner_boundary_nodes=boundary_sort(inner_geometry_points)
    inner_boundary_nodes=[tuple(inner_boundary_nodes[i]) for i in range(len(inner_boundary_nodes))]
    
    Nodes=[]
    
    for i in range(len(inner_boundary_nodes)):
        count=0
    
        tmpo1=(inner_boundary_nodes[i][0]-1,inner_boundary_nodes[i][1])
        if tmpo1 not in inner_geometry_points:
            Nodes.append(tmpo1)
        
        tmpo2=(inner_boundary_nodes[i][0]+1,inner_boundary_nodes[i][1])
        if tmpo2 not in inner_geometry_points:
            Nodes.append(tmpo2)
        
        tmpo3=(inner_boundary_nodes[i][0],inner_boundary_nodes[i][1]-1)
        if tmpo3 not in inner_geometry_points:
            Nodes.append(tmpo3)
        
        tmpo4=(inner_boundary_nodes[i][0],inner_boundary_nodes[i][1]+1)
        if tmpo4 not in inner_geometry_points:
            Nodes.append(tmpo4)
        
        tmpo5=(inner_boundary_nodes[i][0]-1,inner_boundary_nodes[i][1]+1)
        if tmpo5 not in inner_geometry_points:
            Nodes.append(tmpo5)
        
        tmpo6=(inner_boundary_nodes[i][0]-1,inner_boundary_nodes[i][1]-1)
        if tmpo6 not in inner_geometry_points:
            Nodes.append(tmpo6)
        
        tmpo7=(inner_boundary_nodes[i][0]+1,inner_boundary_nodes[i][1]+1)
        if tmpo7 not in inner_geometry_points:
            Nodes.append(tmpo7)
        
        tmpo8=(inner_boundary_nodes[i][0]+1,inner_boundary_nodes[i][1]-1)
        if tmpo8 not in inner_geometry_points:
            Nodes.append(tmpo8)

        
    Nodes=list(set(Nodes))
    
    # sorting
    
    center_x=0
    center_y=0
    for i in range(len(Nodes)):
        center_x=center_x+Nodes[i][0]
        center_y=center_y+Nodes[i][1]
        
    center_x=center_x/len(Nodes)
    center_y=center_y/len(Nodes)
    center=(center_x,center_y)    

    angles=[]
    for i in range(len(Nodes)):
        rel_x=Nodes[i][0]-center[0]
        rel_y=Nodes[i][1]-center[1]
        angles.append(math.atan2(rel_y,rel_x))
        
    angles=np.array(angles)
    
    bd=np.array(Nodes)
    sorted_one=bd[np.argsort(angles)]
    
    just_outer_nodes_of_inner_bdry=[]
    for i in range(len(sorted_one)):
        just_outer_nodes_of_inner_bdry.append(tuple(sorted_one[i]))
    
    just_outer_nodes_of_inner_bdry=np.array(just_outer_nodes_of_inner_bdry)
    return just_outer_nodes_of_inner_bdry


def screwed_up(flaglist):
	flaglist=flaglist.tolist()

	if flaglist==[0 for _ in range(len(flaglist))]:
		return False
    
	nz_indices=[i for i,j in enumerate(flaglist) if j!=0 ]
	if nz_indices==list(range(nz_indices[0],nz_indices[-1]+1)):
		return False
    
	elif flaglist[0] ==0:
		return True
    
	else:
		for i in range(len(flaglist)):
			if flaglist[i]==0:
				first_zero=i
				break
        
		shifted=flaglist[first_zero:]+flaglist[:first_zero]
        
		nz_indices2=[i for i,j in enumerate(shifted) if j!=0 ]
		if nz_indices2==list(range(nz_indices2[0],nz_indices2[-1]+1)):
			return False
		else:
			return True



def starting_points(points, dom_dict): # [(0, 0), (1, 0), (2, 0), (0, 1), (2, 1), (0, 2), (1, 2), (2, 2)]

    x_coord=np.array([points[i][0] for i in range(len(points))])
    y_coord=np.array([points[i][1] for i in range(len(points))])

    min_x=np.min(x_coord)
    max_x=np.max(x_coord)
    min_y=np.min(y_coord)
    max_y=np.max(y_coord)

    center_x=(min_x+max_x)/2
    center_y=(min_y+max_y)/2

    tuned_points=[(x_coord[i]-center_x,(y_coord[i]-center_y)) for i in range(len(points))]

    tp=np.array(tuned_points)

    xflipped=copy.deepcopy(tp)
    xflipped[:,0]=-xflipped[:,0]
    xflipped=[(xflipped[i][0],xflipped[i][1]) for i in range(len(xflipped))]

    yflipped=copy.deepcopy(tp)
    yflipped[:,1]=-yflipped[:,1]
    yflipped=[(yflipped[i][0],yflipped[i][1]) for i in range(len(yflipped))]

    # xy_flipped=copy.deepcopy(tp)
    # xy_flipped=-xy_flipped
    # xy_flipped=[(xy_flipped[i][0],xy_flipped[i][1]) for i in range(len(xy_flipped))] 

    diag_flipped=copy.deepcopy(tp)
    diag_flipped=np.flip(diag_flipped,axis=1)
    diag_minus_flipped=-diag_flipped

    diag_flipped=[(diag_flipped[i][0],diag_flipped[i][1]) for i in range(len(diag_flipped))]
    diag_minus_flipped=[(diag_minus_flipped[i][0],diag_minus_flipped[i][1]) for i in range(len(diag_minus_flipped))]

    ################################################################

    Yaxis_SYMMETRY= (set(tuned_points)==set(xflipped))
    Xaxis_SYMMETRY= (set(tuned_points)==set(yflipped))

    # Origin_SYMMETRY= (set(tuned_points)==set(xy_flipped))  # Origin point symmetry

    YeqX_SYMMETRY= (set(tuned_points)==set(diag_flipped))
    YeqminusX_SYMMETRY= (set(tuned_points)==set(diag_minus_flipped))  # Y=-X symmetry


    points_reduced=copy.deepcopy(tuned_points)

    ################################################################


    if Xaxis_SYMMETRY and Yaxis_SYMMETRY:  

        points_reduced=[(points_reduced[i][0], points_reduced[i][1]) \
                        for i in range(len(points_reduced)) if points_reduced[i][0]<=0 and points_reduced[i][1]<=0 ]

        if YeqX_SYMMETRY:

            points_reduced=[(points_reduced[i][0], points_reduced[i][1]) \
                            for i in range(len(points_reduced)) if points_reduced[i][1]<=points_reduced[i][0]]


        else:

            pass


    elif YeqX_SYMMETRY and YeqminusX_SYMMETRY:

        points_reduced=[(points_reduced[i][0], points_reduced[i][1]) \
                        for i in range(len(points_reduced)) \
                        if points_reduced[i][0]<=points_reduced[i][1] and points_reduced[i][1]<=-points_reduced[i][0] ]



    elif Xaxis_SYMMETRY:

        points_reduced=[(points_reduced[i][0], points_reduced[i][1]) \
         for i in range(len(points_reduced)) if points_reduced[i][1]<=0]


    elif Yaxis_SYMMETRY:

        points_reduced=[(points_reduced[i][0], points_reduced[i][1]) \
         for i in range(len(points_reduced)) if points_reduced[i][0]<=0]



    elif YeqX_SYMMETRY:

        points_reduced=[(points_reduced[i][0], points_reduced[i][1]) \
                        for i in range(len(points_reduced)) \
                        if points_reduced[i][0]<=points_reduced[i][1] ]


    elif YeqminusX_SYMMETRY:

        points_reduced=[(points_reduced[i][0], points_reduced[i][1]) \
                        for i in range(len(points_reduced)) \
                        if points_reduced[i][1]<=-points_reduced[i][0] ]


    else:

        pass


    points_reduced=[(int(points_reduced[i][0]+center_x),int(points_reduced[i][1]+center_y)) for i in range(len(points_reduced))]

    points_to_start=[dom_dict[points_reduced[i]] for i in range(len(points_reduced))]

    return points_to_start


def getting_Sorted_Boundaries(outer_geo_points, inner_geos_points, dom_dict):

    SB=[boundary_sort(outer_geo_points)]
    for i in range(len(inner_geos_points)):
        SB.append(just_outer_nodes_sorted(inner_geos_points[i]))

    Srtd_Bdries=[]
    for i in range(len(SB)):
        Srtd_Bdries.append( [dom_dict[tuple(SB[i][j])] for j in range(len(SB[i]))]  )

    return Srtd_Bdries



def getting_inner_nodes_of_material_grid(points):
# Filtering out the inner nodes of material grid points in geometry (for calculating gradient)
    
    inner_nodes=[]

    for i in range(len(points)):
    
        x_crd=points[i][0]
        y_crd=points[i][1]
        
        if (x_crd+1,y_crd) in points:
            pass
        else:
            continue
        
        if (x_crd-1,y_crd) in points:
            pass
        else:
            continue
            
        if (x_crd,y_crd+1) in points:
            pass
        else:
            continue
            
        if (x_crd,y_crd-1) in points:
            pass
        else:
            continue
            
        inner_nodes.append(points[i])

    x_arr=np.array(inner_nodes).T[0]
    y_arr=np.array(inner_nodes).T[1]

    return inner_nodes, x_arr, y_arr


print('All functions imported.')



