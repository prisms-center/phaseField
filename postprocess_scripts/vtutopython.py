try:
    import sys,getopt,os,meshio,pickle
    import numpy as np
    import pandas as pd
except:
    print('meshio, pandas, pickle, and numpy packages are required.')

argumentList = sys.argv


# Options
options = "i:o:hs:v:d:"
 
# Long options
long_options = ["input","output","help","scalars", "vectors", "dimension"]
vectors = []
scalars = []
ndim = 2
grid = False
infile = ''
outfile = ''

try:
    # Parsing argument
    arguments, values = getopt.getopt(argumentList[1:], options, long_options)
    # checking each argument
    for currentArgument, currentValue in arguments:
 
        if currentArgument in ("-h", "--help"):
            print ("\
\nUsage: vtutopython.py -i <infile> -o <outfile> -s <scalars> -v <vectors> -d <dims> -g -h \n\n\
\
-i or --input: specify input file in vtu or vtk format\n\
-o or --output: specify output file in pkl format\n\
-s or --scalars: comma separated list of scalar fields (no spaces)\n\
-v or --vectors: comma separated list of vector fields (no spaces)\n\
-d or --dimension: number of dimensions, either 2 or 3 \n\
-h or --help: access this help screen\n\n\
\
For an easier experience, copy the core function, vtutopython, directly into your script.\n\n\
\
For this version, the input should be an evenly spaced rectilnear grid (no AMR, same spacing in every direction)\n\n\
\
To extract file, use:\n\
    with open('name_of_file.pkl', 'rb') as f:\n\
        loaded_file = pickle.load(f)\n")
            sys.exit(0)

        elif currentArgument in ("-i", "--input"):
            infile = currentValue
            if infile.split(sep='.')[1] not in ['vtu','vtk']:
                print('Input file should be in vtu or vtk format')
                sys.exit(1)
        elif currentArgument in ("-o", "--output"):
            outfile = currentValue
            if outfile.split(sep='.')[1] != 'pkl':
                print('Output file should be in pkl format')
                sys.exit(1)
        elif currentArgument in ("-v", "--vectors"):
            vectors = currentValue.replace(' ','').split(',')
        elif currentArgument in ("-s", "--scalars"):
            scalars = currentValue.replace(' ','').split(',')
        elif currentArgument in ("-d", "--dimension"):
            ndim = currentValue
             
except getopt.error as err:
    # output error, and return with an error code
    print (str(err))

if infile == '':
    print('Specify an input file.')
    sys.exit(1)

if outfile == '':
    print('Specify an output file.')
    sys.exit(1)

if scalars == [] and vectors == []:
    print('Specify some fields to output.')
    sys.exit(1)

def vtutopython(filename, scalars=[], vectors=[], dim=2):
    # for this version, the input should be an evenly spaced rectilnear grid (no AMR, same spacing in every direction)

    # filename is the name of the input vtu file (vtk probably works too) in the current directory 
    # scalars is a list of strings containing the names of scalar fields
    # vectors is a list of strings containing the names of vector fields
    # dim is the dimension of the problem (2 or 3)
 
    # returns a dictionary containing numpy arrays of the requested fields

    # required libraries

    vecindices = []
    indices = []
    outputs = {}

    mesh = meshio.read(os.path.join(os.getcwd(),filename))      # load the mesh

    if dim == 2:
        indexcount = 2

        # put the data into a pandas dataframe to process it (which requires it be loaded into a dictionary first)
        data = {'x': mesh.points[:,0], 'y': mesh.points[:,1]}       
        for field in vectors:
            data[field + 'x'] = mesh.point_data[field][:,0]
            data[field + 'y'] = mesh.point_data[field][:,1]
            vecindices.append((indexcount,field))                   # we will put the fields in a numpy array, so need to know what the indices will be
            indexcount += 2
        for field in scalars:
            data[field] = mesh.point_data[field]
            indices.append((indexcount,field))
            indexcount += 1
        df = pd.DataFrame(data=data)                                # initialize dataframe
        df = df.drop_duplicates()                                   # all values are repeated for each cell node, need to get rid of them
        df = df.sort_values(['x', 'y'], ascending=[True, True])     # sort by position so that we can easily reshape the array
        cleandata = df.to_numpy()                                   # we need to put things in a numpy array so they can be reshaped

        dh = cleandata[1,1]-cleandata[0,1]                          # find spacing between cells
        xsize = cleandata[-1,0]                                     # domain size
        ysize = cleandata[-1,1]
        xdiv = round(xsize/dh)+1                                    # number of points in each direction (number of cells + 1)
        ydiv = round(ysize/dh)+1

        # reshape the arrays and load them into a dictionary to be returned
        outputs['x'] = cleandata[:,0].reshape(xdiv,ydiv)
        outputs['y'] = cleandata[:,1].reshape(xdiv,ydiv)

        for i in indices:
            outputs[i[1]] = cleandata[:,i[0]].reshape(xdiv,ydiv)
        for i in vecindices:
            outputs[i[1]] = np.array([cleandata[:,i[0]].reshape(xdiv,ydiv),cleandata[:,i[0]+1].reshape(xdiv,ydiv)])

        return outputs

    elif dim == 3:
        indexcount = 3

        # put the data into a pandas dataframe to process it (which requires it be loaded into a dictionary first)
        data = {'x': mesh.points[:,0], 'y': mesh.points[:,1], 'z': mesh.points[:,2]}
        for field in vectors:
            data[field + 'x'] = mesh.point_data[field][:,0]
            data[field + 'y'] = mesh.point_data[field][:,1]
            data[field + 'z'] = mesh.point_data[field][:,2]
            vecindices.append((indexcount,field))                               # we will put the fields in a numpy array, so need to know what the indices will be
            indexcount += 3
        for field in scalars:
            data[field] = mesh.point_data[field]
            indices.append((indexcount,field))
            indexcount += 1
        df = pd.DataFrame(data=data)                                            # initialize dataframe
        df = df.drop_duplicates()                                               # all values are repeated for each cell node, need to get rid of them
        df = df.sort_values(['x', 'y', 'z'], ascending=[True, True, True])      # sort by position so that we can easily reshape the array
        cleandata = df.to_numpy()                                               # we need to put things in a numpy array so they can be reshaped

        dh = cleandata[1,2]-cleandata[0,2]                                      # find spacing between cells
        xsize = cleandata[-1,0]                                                 # domain size
        ysize = cleandata[-1,1]
        zsize = cleandata[-1,2]
        xdiv = round(xsize/dh)+1                                                # number of points in each direction (number of cells + 1)
        ydiv = round(ysize/dh)+1
        zdiv = round(zsize/dh)+1

        # reshape the arrays and load them into a dictionary to be returned
        outputs['x'] = cleandata[:,0].reshape(xdiv,ydiv,zdiv)
        outputs['y'] = cleandata[:,1].reshape(xdiv,ydiv,zdiv)
        outputs['z'] = cleandata[:,2].reshape(xdiv,ydiv,zdiv)

        for i in indices:
            outputs[i[1]] = cleandata[:,i[0]].reshape(xdiv,ydiv,zdiv)
        for i in vecindices:
            outputs[i[1]] = np.array([cleandata[:,i[0]].reshape(xdiv,ydiv,zdiv),cleandata[:,i[0]+1].reshape(xdiv,ydiv,zdiv),cleandata[:,i[0]+2].reshape(xdiv,ydiv,zdiv)])

        return outputs
    
with open(outfile, 'wb') as f:
    pickle.dump(vtutopython(filename=infile,scalars=scalars,vectors=vectors,dim=ndim), f)
