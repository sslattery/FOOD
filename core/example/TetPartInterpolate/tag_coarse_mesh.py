# pytaps script to tag the fine mesh with data and setup a tag without data 
# on the coarse mesh to interpolate onto.

import math
from itaps import iBase, iMesh

mesh = iMesh.Mesh()
mesh.load("coarse_tet_part.vtk")

fine_tag = mesh.createTag("range", 1, float)

vertices = mesh.getEntities( iBase.Type.vertex )
coords = mesh.getVtxCoords(vertices)

tag_data = []

for i in xrange len(coords):
    
    tag_data += [ 0.0 ]


fine_tag[vertices] = tag_data
