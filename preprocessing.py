from functions import *



N_Laser_grid = 10

N_Material_grid=3*N_Laser_grid+2 # Because we already configured laser-material grid system like this

Laser_Diameter=0.005 # meter # it is also laser grid gap size

dL=Laser_Diameter/3 # Because we already configured laser-material grid system like this
L=(N_Material_grid-1)*dL

SCALE_FACTOR_lasergrid = Laser_Diameter
SCALE_FACTOR_materialgrid = dL

lg_padding=2*dL   # Because we already configured laser-material grid system like this

###########################################################

pt=list(map(int, np.linspace(0,N_Laser_grid-1,N_Laser_grid)))  
xs, ys = np.meshgrid(pt, pt)

pt_mg=list(map(int, np.linspace(0,N_Material_grid-1,N_Material_grid)))  
xs_mg, ys_mg = np.meshgrid(pt_mg, pt_mg)

# ♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣ #
# geo_1 MUST BE the outer geometry!!!

def geometry1(x,y):  # outer_boundary

    # circle
    # xx=x-L/2
    # yy=y-L/2
    # rr=0.5*L
   
    # gg= xx**2+yy**2<=rr**2

    # square
    gg = ( (0<=x) & (x<=L) & (0<=y) & (y<=L))


    hole1 = ( (0<=x) & (x<=1*L/8) & (0<=y) & (y<=1*L/8))
    hole2 = ( (0<=x) & (x<=1*L/8) & (7*L/8<=y) & (y<=L))
    hole3 = ( (7*L/8<=x) & (x<=L) & (0<=y) & (y<=1*L/8))
    hole4 = ( (7*L/8<=x) & (x<=L) & (7*L/8<=y) & (y<=L))

    # hole1 = ( (4*L/10<=x) & (x<=6*L/10) & (0*L/10<=y) & (y<=2*L/10))
    # hole2 = ( (4*L/10<=x) & (x<=6*L/10) & (8*L/10<=y) & (y<=10*L/10))
    # hole3 = ( (0*L/10<=x) & (x<=2*L/10) & (4*L/10<=y) & (y<=6*L/10))
    # hole4 = ( (8*L/10<=x) & (x<=10*L/10) & (4*L/10<=y) & (y<=6*L/10))

    gg_final = gg & (~hole1) & (~hole2) & (~hole3) & (~hole4)
    # gg_final = gg

    return gg_final
   
    
def geometry2(x,y):  # hollow
    
    # square
    gg = ( (2*L/10<=x) & (x<=8*L/10) & (2*L/10<=y) & (y<=8*L/10))

    hole1 = ( (2*L/10<=x) & (x<=4*L/10) & (2*L/10<=y) & (y<=4*L/10))
    hole2 = ( (2*L/10<=x) & (x<=4*L/10) & (6*L/10<=y) & (y<=8*L/10))
    hole3 = ( (6*L/10<=x) & (x<=8*L/10) & (2*L/10<=y) & (y<=4*L/10))
    hole4 = ( (6*L/10<=x) & (x<=8*L/10) & (6*L/10<=y) & (y<=8*L/10))

    gg_final = gg & (~hole1) & (~hole2) & (~hole3) & (~hole4)

    return gg_final


######### Copy and paste the followings and change the name of geometry functions #########

geo_1=geometry1(SCALE_FACTOR_lasergrid*xs+lg_padding, SCALE_FACTOR_lasergrid*ys+lg_padding)  
geo_2=geometry2(SCALE_FACTOR_lasergrid*xs+lg_padding, SCALE_FACTOR_lasergrid*ys+lg_padding)

geo_mg1=geometry1(SCALE_FACTOR_materialgrid*xs_mg, SCALE_FACTOR_materialgrid*ys_mg)  
geo_mg2=geometry2(SCALE_FACTOR_materialgrid*xs_mg, SCALE_FACTOR_materialgrid*ys_mg)


## Put the arguments of the following functions having all geometries (geo_1 (outer geometry) and hollow ones) ##

geo, domain_dict, Total_points, Total_points_coord, total_num_points=making_dict(xs,ys,geo_1,geo_2)
outer_geo_points, *inner_geos_points = get_all_points(xs,ys,geo_1,geo_2)

geo_mg, domain_dict_mg, Total_points_mg, Total_points_mg_coord, total_num_points_mg = making_dict(xs_mg,ys_mg,geo_mg1,geo_mg2)

# ♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣♣ #

####################################################################################################################

# For calculating thermal gradient in the future
inner_nodes_mg, x_arr, y_arr = getting_inner_nodes_of_material_grid(Total_points_mg_coord) 
Sorted_Boundaries=getting_Sorted_Boundaries(outer_geo_points, inner_geos_points, domain_dict)
The_number_of_boundaries=len(Sorted_Boundaries)

POINTS_TO_START=starting_points(Total_points_coord,domain_dict)

inner_points_lg=copy.deepcopy(Total_points)

for kk in range(The_number_of_boundaries):
    inner_points_lg=set(inner_points_lg)-set(Sorted_Boundaries[kk])
    
inner_points_lg=list(inner_points_lg)


print(geo)
print()
print(Total_points)
print()

print(POINTS_TO_START)
print()
print(total_num_points)
print()



####################################################################################################################


blacklist_left=[N_Laser_grid*i for i in range(1,N_Laser_grid)]
blacklist_right=[N_Laser_grid*i-1 for i in range(1,N_Laser_grid)]


def NextAvailablePoints(CurrentPoint, NeverGonePoints):  

    # CurrentPoint: 7
    # NeverGonePoints: [8,9,10,13,14,15,16,19,20,21,22,25,26,27,28]
    
    coord_x=CurrentPoint%N_Laser_grid
    coord_y=CurrentPoint//N_Laser_grid

    
    if CurrentPoint in blacklist_left:
        NextPoints=[]

        p2=(coord_x+1) + N_Laser_grid*coord_y
        if p2 in NeverGonePoints:

            NextPoints.append(p2)

        p3=coord_x + N_Laser_grid*(coord_y-1)
        if p3 in NeverGonePoints:

            NextPoints.append(p3)

        p4=coord_x + N_Laser_grid*(coord_y+1)
        if p4 in NeverGonePoints:

            NextPoints.append(p4)
            
    elif CurrentPoint in blacklist_right:
        NextPoints=[]
        
        p1=(coord_x-1) + N_Laser_grid*coord_y
        if p1 in NeverGonePoints:

            NextPoints.append(p1)

        p3=coord_x + N_Laser_grid*(coord_y-1)
        if p3 in NeverGonePoints:

            NextPoints.append(p3)

        p4=coord_x + N_Laser_grid*(coord_y+1)
        if p4 in NeverGonePoints:

            NextPoints.append(p4)
            
    else:
        NextPoints=[]

        p1=(coord_x-1) + N_Laser_grid*coord_y
        if p1 in NeverGonePoints:

            NextPoints.append(p1)

        p2=(coord_x+1) + N_Laser_grid*coord_y
        if p2 in NeverGonePoints:

            NextPoints.append(p2)

        p3=coord_x + N_Laser_grid*(coord_y-1)
        if p3 in NeverGonePoints:

            NextPoints.append(p3)

        p4=coord_x + N_Laser_grid*(coord_y+1)
        if p4 in NeverGonePoints:

            NextPoints.append(p4)        
        
        
    return NextPoints





def flag_record_filtering(points_to_go, flag_rec):
	to_be_deleted_index=[]

	for kk in range(The_number_of_boundaries):
		sb=Sorted_Boundaries[kk]
		flr=flag_rec[kk]  # taking advantage of Python's property

		filtered_indices=[i for i,pt in enumerate(points_to_go) if pt in sb]

		for i in filtered_indices:
			flr[i][sb.index(points_to_go[i])]=1  # taking advantage of Python's property
			if screwed_up(flr[i]):
				to_be_deleted_index.append(i)

	return to_be_deleted_index



def self_collision_filtering(path_rec):
	to_be_deleted_index=[]

	filtered_indices=[i for i,a in enumerate(path_rec) if a[-1] in inner_points_lg]

	for i in filtered_indices:
		path_now=path_rec[i]

		path_end=path_now[-1]
		path_second_end=path_now[-2]

		direction=path_end-path_second_end
		forth_middle=path_end+direction

		letsgo=N_Laser_grid+1-abs(direction)
		seawall={forth_middle-letsgo,forth_middle,forth_middle+letsgo}
		fire_egg={path_end-letsgo,path_end+letsgo}

		spr=set(path_now)

		if seawall-spr != seawall and fire_egg-spr == fire_egg:
			to_be_deleted_index.append(i)

	return to_be_deleted_index



def assassinator(path_rec, points_to_go, flag_rec, depth):
	to_be_deleted_index=[]

	to_be_deleted_index_1=flag_record_filtering(points_to_go, flag_rec)
	to_be_deleted_index+=to_be_deleted_index_1

	if depth>6:
		to_be_deleted_index_2=self_collision_filtering(path_rec)
		to_be_deleted_index+=to_be_deleted_index_2

	points_to_go=np.delete(points_to_go,to_be_deleted_index)

	path_rec=np.delete(path_rec,to_be_deleted_index,axis=0)

	for kk in range(The_number_of_boundaries):
		flag_rec[kk]=np.delete(flag_rec[kk],to_be_deleted_index,axis=0)


	return path_rec, points_to_go, flag_rec



print('Preprocessing done.')