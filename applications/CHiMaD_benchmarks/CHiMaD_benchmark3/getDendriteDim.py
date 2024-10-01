import sys
import math

# ----------------------------------------------------------------------------------------
#
# Splits a string into a list of strings.
#
def SplitLines(s):
    lines = s.split("\n")
    return [line for line in lines if line != '']


#
# Extracts values from VisIt's pick output.
#
def ExtractPickValue(s):
    lines = SplitLines(s)
    retval = 0.
    #print lines
    if "<nodal>" in lines[5]:
        sum = 0.
        nSum = 0
        # Extract all of the values in the pick and average them. Not
        # exactly what we want but it's something.
        for line in lines[2:]:
            e = line.find("=")
            if e != -1:
                tok = line[e+2:]
                value = float(tok)
                sum += value
                nSum += 1
        retval = sum / float(nSum)
    elif "<zonal>" in lines[1]:
        e = lines[1].find("=")
        if e != -1:
            tok = lines[1][e+2:]
            retval = float(tok)
    return retval


# ----------------------------------------------------------------------------------------
#
def interfaceSearch(xIncrement, yIncrement, zIncrement, pickLocX, pickLocY, pickLocZ):
	intLoc = [0.0, 0.0];
	found_int = False

	NodePick((pickLocX,pickLocY,pickLocZ))
	pickVal = ExtractPickValue(GetPickOutput())

	while found_int != True:

		pickLocX += xIncrement
		pickLocY += yIncrement
		pickLocZ += zIncrement

		if (abs(xIncrement) > 0.0) and (yIncrement == 0.0) and (zIncrement == 0.0):
				if xIncrement > 0.0:
					if pickLocX > span[1]:
						pickLocX -= xIncrement
						return [pickLocX, pickLocX]
				else:
					if pickLocX < span[0]:
						pickLocX -= xIncrement
						return [pickLocX, pickLocX]
		elif (abs(yIncrement) > 0.0) and (xIncrement == 0.0) and (zIncrement == 0.0):
				if yIncrement > 0.0:
					if pickLocY > span[3]:
						pickLocY -= yIncrement
						return [pickLocY, pickLocY]
				else:
					if pickLocY < span[2]:
						pickLocY -= yIncrement
						return [pickLocY, pickLocY]
		elif (abs(zIncrement) > 0.0) and (xIncrement == 0.0) and (yIncrement == 0.0):
				if zIncrement > 0.0:
					if pickLocZ > span[5]:
						pickLocZ -= zIncrement
						return [pickLocZ, pickLocZ]
				else:
					if pickLocZ < span[0]:
						pickLocZ -= zIncrement
						return [pickLocZ, pickLocZ]
		else:
				print 'Error: step is not aligned with the mesh'
				found_int = True

		NodePick((pickLocX,pickLocY,pickLocZ))
		pickVal_next = ExtractPickValue(GetPickOutput())

		if (pickVal > 0.0) and (pickVal_next <= 0.0):
			if (abs(xIncrement) > 0.0) and (abs(yIncrement) < 1.0e-8) and (abs(zIncrement) < 1.0e-8):
				intLoc[0] = (0.0-pickVal)*(xIncrement)/(pickVal_next-pickVal) + (pickLocX-xIncrement)
				intLoc[1] = pickLocX-xIncrement
				found_int = True
			elif (abs(yIncrement) > 0.0) and (xIncrement == 0.0) and (zIncrement == 0.0):
				intLoc[0] = (0.0-pickVal)*(yIncrement)/(pickVal_next-pickVal) + (pickLocY-yIncrement)
				intLoc[1] = pickLocY-yIncrement
				found_int = True
			elif (abs(zIncrement) > 0.0) and (xIncrement == 0.0) and (yIncrement == 0.0):
				intLoc[0] = (0.0-pickVal)*(zIncrement)/(pickVal_next-pickVal) + (pickLocZ-zIncrement)
				intLoc[1] = pickLocZ-zIncrement
				found_int = True
			else:
				print 'Error: step is not aligned with the mesh'
				found_int = True

		pickVal = pickVal_next

	return intLoc

# ----------------------------------------------------------------------------------------
#
#  Searches for the interface, moving from the centroid outward
#
def findInterface2(xIncrement, yIncrement, zIncrement):
	pickLocX = centroid[0]
	pickLocY = centroid[1]
	pickLocZ = centroid[2]

	skipNodes = 8.0

	while skipNodes >= 1.0:

		intLoc = interfaceSearch(xIncrement*skipNodes, yIncrement*skipNodes, zIncrement*skipNodes, pickLocX, pickLocY, pickLocZ)

		if (abs(xIncrement) > 0.0) and (yIncrement == 0.0) and (zIncrement == 0.0):
				pickLocX = intLoc[1]
		elif (abs(yIncrement) > 0.0) and (xIncrement == 0.0) and (zIncrement == 0.0):
				pickLocY = intLoc[1]
		elif (abs(zIncrement) > 0.0) and (xIncrement == 0.0) and (yIncrement == 0.0):
				pickLocZ = intLoc[1]
		else:
				print 'Error: step is not aligned with the mesh'
				found_int = True

		skipNodes /= 2

	return intLoc[0]

# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
#  Searches for the interface, moving from the centroid outward
#
def findInterface3(xIncrement, yIncrement, zIncrement,starting_point,var):
	pickLocX = starting_point[0]
	pickLocY = starting_point[1]
	pickLocZ = starting_point[2]

	skipNodes = 8.0

	while skipNodes >= 1.0:
		(intLoc, last_point) = interfaceSearch3(xIncrement*skipNodes, yIncrement*skipNodes, zIncrement*skipNodes, pickLocX, pickLocY, pickLocZ, var)

		pickLocX = last_point[0]
		pickLocY = last_point[1]
		pickLocZ = last_point[2]


		skipNodes /= 2

	return intLoc

# ----------------------------------------------------------------------------------------
#
def interfaceSearch3(xIncrement, yIncrement, zIncrement, pickLocX, pickLocY, pickLocZ ,var):
	found_int = False

	pickVal = arbitraryPick((pickLocX,pickLocY,pickLocZ),var)

	while found_int != True:

		pickLocX += xIncrement
		pickLocY += yIncrement
		pickLocZ += zIncrement

		pickVal_next = arbitraryPick((pickLocX,pickLocY,pickLocZ),var)

		if (pickVal > 0.5) and (pickVal_next <= 0.5):
			step_size = math.sqrt(xIncrement*xIncrement + yIncrement*yIncrement + zIncrement*zIncrement)
			dist_to_int = (0.5-pickVal)*(step_size)/(pickVal_next-pickVal)

			intLoc = (pickLocX-xIncrement + xIncrement/step_size * dist_to_int, pickLocY-yIncrement + yIncrement/step_size * dist_to_int, pickLocZ-zIncrement + zIncrement/step_size * dist_to_int)
			found_int = True


		pickVal = pickVal_next

	return (intLoc, (pickLocX-xIncrement,pickLocY-yIncrement,pickLocZ-zIncrement))

# ----------------------------------------------------------------------------------------
#
def arbitraryPick(pick_coords,var):
# Trilinear interpolation algorithm for finding the value of a field at an arbitrary point

    ZonePick(pick_coords)
    temp = ExtractPickValue(GetPickOutput())
    temp_obj = GetPickOutputObject()
    node_list = temp_obj['incident_nodes']
    print temp_obj

    node_coords = []
    node_val = []

    for node in node_list:
        PickByNode(element=node)
        #PickByNode(element=node,domain=temp_obj["domain_id"])
        temp_obj = GetPickOutputObject()
        node_coords.append(temp_obj['point'])
        node_val.append(temp_obj[var])

	# Cycle through the nodes to determine the relationship between each node and the point of interest
    pX_pY_pZ_assigned = False
    pX_mY_pZ_assigned = False
    pX_pY_mZ_assigned = False
    pX_mY_mZ_assigned = False
    mX_pY_pZ_assigned = False
    mX_mY_pZ_assigned = False
    mX_pY_mZ_assigned = False
    mX_mY_mZ_assigned = False

    node_index_list = [0,1,2,3,4,5,6,7]
    for node_index in node_index_list:
        if ( abs(node_coords[node_index][0] - pick_coords[0]) < 1.0e-6):
            oX = True
            pX = False
            mX = False
        elif (node_coords[node_index][0] >= pick_coords[0]):
            oX = False
            pX = True
            mX = False
        else:
            oX = False
            pX = False
            mX = True

        if ( abs(node_coords[node_index][1] - pick_coords[1]) < 1.0e-6):
            oY = True
            pY = False
            mY = False
        elif node_coords[node_index][1] >= pick_coords[1]:
            oY = False
            pY = True
            mY = False
        else:
            oY = False
            pY = False
            mY = True

        if ( abs(node_coords[node_index][2] - pick_coords[2]) < 1.0e-6):
            oZ = True
            pZ = False
            mZ = False
        elif node_coords[node_index][2] >= pick_coords[2]:
            oZ = False
            pZ = True
            mZ = False
        else:
            oZ = False
            pZ = False
            mZ = True

        if ((pX and pY) and pZ):
            pX_pY_pZ_assigned = True
            node_coords_pX_pY_pZ = node_coords[node_index]
            node_val_pX_pY_pZ = node_val[node_index]

        if ((pX and mY) and pZ):
            pX_mY_pZ_assigned = True
            node_coords_pX_mY_pZ = node_coords[node_index]
            node_val_pX_mY_pZ = node_val[node_index]

        if ((pX and pY) and mZ):
            pX_pY_mZ_assigned = True
            node_coords_pX_pY_mZ = node_coords[node_index]
            node_val_pX_pY_mZ = node_val[node_index]

        if ((pX and mY) and mZ):
            pX_mY_mZ_assigned = True
            node_coords_pX_mY_mZ = node_coords[node_index]
            node_val_pX_mY_mZ = node_val[node_index]

        if ((mX and pY) and pZ):
            mX_pY_pZ_assigned = True
            node_coords_mX_pY_pZ = node_coords[node_index]
            node_val_mX_pY_pZ = node_val[node_index]

        if ((mX and mY) and pZ):
            mX_mY_pZ_assigned = True
            node_coords_mX_mY_pZ = node_coords[node_index]
            node_val_mX_mY_pZ = node_val[node_index]

        if ((mX and pY) and mZ):
            mX_pY_mZ_assigned = True
            node_coords_mX_pY_mZ = node_coords[node_index]
            node_val_mX_pY_mZ = node_val[node_index]

        if ((mX and mY) and mZ):
            mX_mY_mZ_assigned = True
            node_coords_mX_mY_mZ = node_coords[node_index]
            node_val_mX_mY_mZ = node_val[node_index]

	# Cycle through again to assign the nodes that share a coordinate with the pick coordinates
	for node_index in node_index_list:
		if ( abs(node_coords[node_index][0] - pick_coords[0]) < 1.0e-6):
			oX = True
			pX = False
			mX = False
		elif (node_coords[node_index][0] >= pick_coords[0]):
			oX = False
			pX = True
			mX = False
		else:
			oX = False
			pX = False
			mX = True

		if ( abs(node_coords[node_index][1] - pick_coords[1]) < 1.0e-6):
			oY = True
			pY = False
			mY = False
		elif node_coords[node_index][1] >= pick_coords[1]:
			oY = False
			pY = True
			mY = False
		else:
			oY = False
			pY = False
			mY = True

		if ( abs(node_coords[node_index][2] - pick_coords[2]) < 1.0e-6):
			oZ = True
			pZ = False
			mZ = False
		elif node_coords[node_index][2] >= pick_coords[2]:
			oZ = False
			pZ = True
			mZ = False
		else:
			oZ = False
			pZ = False
			mZ = True

		if (pX_pY_pZ_assigned == False):
			if ( ((pX == True) or (oX == True)) and ((pY == True) or (oY == True)) and ((pZ == True) or (oZ == True)) ):
				pX_pY_pZ_assigned = True
				node_coords_pX_pY_pZ = node_coords[node_index]
				node_val_pX_pY_pZ = node_val[node_index]
		if (pX_mY_pZ_assigned == False):
			if ( ((pX == True) or (oX == True)) and ((mY == True) or (oY == True)) and ((pZ == True) or (oZ == True)) ):
				pX_mY_pZ_assigned = True
				node_coords_pX_mY_pZ = node_coords[node_index]
				node_val_pX_mY_pZ = node_val[node_index]
		if (pX_pY_mZ_assigned == False):
			if ( ((pX == True) or (oX == True)) and ((pY == True) or (oY == True)) and ((mZ == True) or (oZ == True)) ):
				pX_pY_mZ_assigned = True
				node_coords_pX_pY_mZ = node_coords[node_index]
				node_val_pX_pY_mZ = node_val[node_index]
		if (pX_mY_mZ_assigned == False):
			if ( ((pX == True) or (oX == True)) and ((mY == True) or (oY == True)) and ((mZ == True) or (oZ == True)) ):
				pX_mY_mZ_assigned = True
				node_coords_pX_mY_mZ = node_coords[node_index]
				node_val_pX_mY_mZ = node_val[node_index]
		if (mX_pY_pZ_assigned == False):
			if ( ((mX == True) or (oX == True)) and ((pY == True) or (oY == True)) and ((pZ == True) or (oZ == True)) ):
				mX_pY_pZ_assigned = True
				node_coords_mX_pY_pZ = node_coords[node_index]
				node_val_mX_pY_pZ = node_val[node_index]
		if (mX_mY_pZ_assigned == False):
			if ( ((mX == True) or (oX == True)) and ((mY == True) or (oY == True)) and ((pZ == True) or (oZ == True)) ):
				mX_mY_pZ_assigned = True
				node_coords_mX_mY_pZ = node_coords[node_index]
				node_val_mX_mY_pZ = node_val[node_index]
		if (mX_pY_mZ_assigned == False):
			if ( ((mX == True) or (oX == True)) and ((pY == True) or (oY == True)) and ((mZ == True) or (oZ == True)) ):
				mX_pY_mZ_assigned = True
				node_coords_mX_pY_mZ = node_coords[node_index]
				node_val_mX_pY_mZ = node_val[node_index]
		if (mX_mY_mZ_assigned == False):
			if ( ((mX == True) or (oX == True)) and ((mY == True) or (oY == True)) and ((mZ == True) or (oZ == True)) ):
				mX_mY_mZ_assigned = True
				node_coords_mX_mY_mZ = node_coords[node_index]
				node_val_mX_mY_mZ = node_val[node_index]

	# Perform the interpolation (using the notation from the Wikipedia article)
	if (abs(node_coords_pX_pY_pZ[0]-node_coords_mX_pY_pZ[0]) > 1.0e-10):
		xd = (pick_coords[0]-node_coords_mX_pY_pZ[0])/(node_coords_pX_pY_pZ[0]-node_coords_mX_pY_pZ[0])
	else:
		xd = 1.0
	if (abs(node_coords_pX_pY_pZ[1]-node_coords_pX_mY_pZ[1]) > 1.0e-10):
		yd = (pick_coords[1]-node_coords_pX_mY_pZ[1])/(node_coords_pX_pY_pZ[1]-node_coords_pX_mY_pZ[1])
	else:
		yd = 1.0
	if (abs(node_coords_pX_pY_pZ[2]-node_coords_pX_pY_mZ[2]) > 1.0e-10):
		zd = (pick_coords[2]-node_coords_pX_pY_mZ[2])/(node_coords_pX_pY_pZ[2]-node_coords_pX_pY_mZ[2])
	else:
		zd = 1.0

	c00 = node_val_mX_mY_mZ * (1.0-xd) + node_val_pX_mY_mZ * xd
	c01 = node_val_mX_mY_pZ * (1.0-xd) + node_val_pX_mY_pZ * xd
	c10 = node_val_mX_pY_mZ * (1.0-xd) + node_val_pX_pY_mZ * xd
	c11 = node_val_mX_pY_pZ * (1.0-xd) + node_val_pX_pY_pZ * xd

	c0 = c00 * (1.0-yd) + c10*yd
	c1 = c01 * (1.0-yd) + c11*yd

	pick_val = c0 * (1.0-zd) + c1*zd

	return pick_val


# ----------------------------------------------------------------------------------------
# SCRIPT BODY
# ----------------------------------------------------------------------------------------


variable = 'phi'

output_num_list = raw_input("List of output numbers: ").split()

directory =  os.path.dirname(os.path.realpath(__file__)) + "/"


for output_num in output_num_list:

    filename = directory+"solution-"+output_num+".vtu"

    OpenDatabase(filename)

    AddPlot("Volume",variable)
    DrawPlots()

    if output_num == output_num_list[0]:
    	Query("SpatialExtents")
    	span = GetQueryOutputValue()
    	print "Span: ", span

    	print "Dendrite tip position, solid fraction:"

    dimension = len(span)/2

    centroid = [0,0,0]

    #numDivX = 90
    #numDivY = 90
    #numDivZ = 90

    numDivX = 128*3*3
    numDivY = 128*3*3
    numDivZ = 128*3*3

    h_x = span[1]/numDivX
    h_y = span[3]/numDivY

    if dimension > 2:
        h_z = span[5]/numDivZ

    starting_point = []
    starting_point.append(0)
    starting_point.append(0)
    starting_point.append(0)

    intLocX_plus = findInterface2(h_x, 0.0, 0.0)

    DefineScalarExpression("shifted_var", "(phi+1)/2")
    ChangeActivePlotsVar("shifted_var")

    Query("Weighted Variable Sum")
    int_var = GetQueryOutputValue()
    solid_frac = int_var/(span[1]*span[3])

    print intLocX_plus, solid_frac

sys.exit()
