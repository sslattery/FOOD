# pytaps script to tag the fine mesh with data and setup a tag without data 
# on the coarse mesh to interpolate onto.

import math
from itaps import iBase, iMesh

mesh = iMesh.Mesh()
mesh.load("fine_tet_part.vtk")

fine_tag = mesh.createTag("domain", 1, float)

num_vert = mesh.getNumOfType( iBase.Type.vertex )
vertices = mesh.getEntities( iBase.Type.vertex )
coords = mesh.getVtxCoords(vertices)

tag_data = []

for i in xrange(num_vert):
    tag_data += [ math.sqrt( coords[i][0]*coords[i][0] + \
                             coords[i][1]*coords[i][1] + \
                             coords[i][2]*coords[i][2] ) ]


fine_tag[vertices] = tag_data

mesh.save("tagged_fine.vtk")
