def parse(fileName):
	"""Parse the information from a .vtk file and append the data to the object.
	   Does not check for the file consistency."""
	file = open(fileName, 'r')

	# Use local data first. 
	points = []
	cells = []
	cell_types = []
	lookup_table = []

	line = ""
	while(True):
		line = file.readline()
		if 'POINTS' in line:
			break

	print "Reading in points from file ", filename
	nPoints = 0

	if (line == ""):
		print "No POINTS defined in the .vtk file"
		return

	# Set the number of points
	nPoints = int(line.split()[1])

	# Append the numbers.
	#for i in xrange(nPoints):
	pointsRead = 0
	while(pointsRead < nPoints):
		#points.append(map(lambda x : float(x), file.readline().split()))
		line = file.readline()
		pointsRead = pointsRead + len(line.split())/3
		points.append(line)

	# Append cells
	line = ""

	while(True):
		line = file.readline()
		if 'CELLS' in line:
			break
			
	print "Reading in cells from file ", filename	
	
	if (line == ""):
		print "No CELLS defined in the .vtk file"
		return
	# Set the number of polygons.
	nCells = int(line.split()[1])
	nVerts = int(line.split()[2])
	
	print "Number of cells:", nCells, " and number of verticies: ", nVerts

	# Read the cells.
	for i in xrange(nCells):
		line = file.readline()
		vert_list = line.split()
		f_verts = map(float,vert_list)
		cells.append(f_verts)

	
	# Read the cell types
	while(True):
		line = file.readline()
		if 'CELL_TYPES' in line:
			break
			
	print "Reading in cell types from file ", filename
	
	nCellTypes = int(line.split()[1])
	
	cellType_size = 0
	while(cellType_size<nCellTypes):
	#for i in xrange(nCellTypes):
		line = file.readline()
		cellType_size = cellType_size + len(line.split())
		cell_types.append(line)
		
	# Read the point data
	while(True):
		line = file.readline()
		if 'POINT_DATA' in line:
			break
			
	print "Reading in points data from file ", filename
	
	nPointData = int(line.split()[1])
	
	# Read the variable name
	while(True):
		line = file.readline()
		if 'SCALARS' in line:
			break
	var_name = (line.split()[1])
	var_type = (line.split()[2])
	var_num = 1 #(line.split()[3])
	
	while(True):
		line = file.readline()
		if 'LOOKUP_TABLE' in line:
			break
			
	lut_size = 0		
	while(lut_size < nPointData):
		line = file.readline()
		lut_size = lut_size + len(line.split())
		for number in line:
			#print number
			lookup_table.append(number)
			#lookup_table.append(float(number))
		
	
	
	file.close()
	
	return [nPoints,points,nCells,nVerts,cells,nCellTypes,cell_types,nPointData,var_name,var_type,var_num,lookup_table]

def write(fileName, append, nPoints,points,nCells,nVerts,cells,nCellTypes,cell_types,nPointData,var_name,var_type,var_num,lookup_table):
	
	# Append or overwrite?
	if append:
		file = open(fileName, 'a')
	else:
		file = open(fileName, 'w')

	file.writelines("# vtk DataFile Version 2.0\n")
	file.writelines("appended vtk files\n")
	file.writelines("ASCII\n")
	file.writelines("DATASET UNSTRUCTURED_GRID\n\n")

	file.writelines("POINTS %d double \n" % nPoints)

	for point in points:
		#file.writelines("%.10f %.10f %.10f\n" % (point[0], point[1], point[2]))
		file.writelines(point)

	
	file.writelines("\nCELLS %d %d \n" % (nCells, nVerts))
	
	for cell in cells:
		#file.writelines(cell)
		for verts in cell:
			for vert in verts:
				file.write("%d " % (vert))
			file.write("\n")
			
		
	file.writelines("\nCELL_TYPES %d \n" % (nCellTypes))
	
	for cell_type in cell_types:
		file.writelines(cell_type)
		
	file.writelines("\nPOINT_DATA %d" % (nPointData))
	file.writelines("\nSCALARS %s %s %s" % (var_name, var_type, var_num))
	file.writelines("\nLOOKUP_TABLE default \n")
	
	for table_line in lookup_table:
		for number in table_line:
			#file.writelines("%s " % (number))
			file.writelines(number)

	# for cell in cells:
# 		file.writelines("%d " % len(polygon))
# 		for label in polygon:
# 			# I don't know how many labels are there in general.
# 			file.writelines("%d " % label)
# 		file.writelines("\n")

	file.close()


# ---------------------------------------------------------------------------------------
# Script body
# ---------------------------------------------------------------------------------------

#filenames = ["parse_test_0.vtk", "parse_test_1.vtk"]
#filenames = ["parse_test_1.vtk"]
#filenames = ["import_test.0.vtk", "import_test.1.vtk" , "import_test.2.vtk" , "import_test.3.vtk"]
filenames =[]
for i in xrange(8):
	filename = "import_test_vtk."+str(i)+".vtk"
	filenames.append(filename)

nPoints = 0
nCells = 0
nVerts = 0
nCellTypes = 0
nPointData = 0

points = []
cells = []
cell_types = []
lookup_table = []



for filename in filenames:
	[nPoints_temp,points_temp,nCells_temp,nVerts_temp,cells_temp,nCellTypes_temp,cell_types_temp,nPointData_temp,var_name,var_type,var_num,lookup_table_temp] = parse(filename)
	
	for i in xrange(len(cells_temp)):
		for j in xrange(len(cells_temp[i])):
			if j > 0:
				cells_temp[i][j] += (float(nCells)*8.0)
	
	nPoints += nPoints_temp
	nCells += nCells_temp
	nVerts += nVerts_temp
	nCellTypes += nCellTypes_temp
	nPointData += nPointData_temp
	points.append(points_temp)
	cells.append(cells_temp)
	cell_types.append(cell_types_temp)
	lookup_table.append(lookup_table_temp)
	
	
write("result.vtk",False,nPoints,points,nCells,nVerts,cells,nCellTypes,cell_types,nPointData,var_name,var_type,var_num,lookup_table)