# pytaps script to tag the fine mesh with data and setup a tag without data 
# on the coarse mesh to interpolate onto.

import math
from itaps import iBase, iMesh

##---------------------------------------------------------------------------##
## Large quadratic hex mesh - tagged with a scalar domain of coordinate 
## distance from origin
mesh = iMesh.Mesh()
mesh.load("hex_mesh.vtk")

tag = mesh.createTag("domain", 1, float)

num_vert = mesh.getNumOfType( iBase.Type.vertex )
vertices = mesh.getEntities( iBase.Type.vertex )
coords = mesh.getVtxCoords(vertices)

tag_data = []

for i in xrange(num_vert):
    tag_data += [ math.cos(coords[i][0]/100.0) * \
                  math.cos(coords[i][1]/100.0) * \
                  math.cos(coords[i][2]/100.0) ]

tag[vertices] = tag_data

mesh.save("hex_domain.vtk")

del mesh
del tag
del vertices
del coords
del tag_data
