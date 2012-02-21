# pytaps script to tag the fine mesh with data and setup a tag without data 
# on the coarse mesh to interpolate onto.

import math
from itaps import iBase, iMesh

##---------------------------------------------------------------------------##
## Large linear tet mesh - tagged with a 3-vector domain of vertex coordinates.
mesh = iMesh.Mesh()
mesh.load("tet_mesh.vtk")

tag = mesh.createTag("domain", 3, float)

num_vert = mesh.getNumOfType( iBase.Type.vertex )
vertices = mesh.getEntities( iBase.Type.vertex )
coords = mesh.getVtxCoords(vertices)

tag_data = []

for i in xrange(num_vert):
    tag_data += [ [coords[i][0]*coords[i][0], \
                   coords[i][1]*coords[i][1], \
                   coords[i][2]*coords[i][2]] ]

tag[vertices] = tag_data

mesh.save("tet_domain.vtk")

del mesh
del tag
del vertices
del coords
del tag_data

##---------------------------------------------------------------------------##
## hex mesh - tagged with a 3-vector range.
mesh = iMesh.Mesh()
mesh.load("hex_mesh.vtk")

tag = mesh.createTag("range", 3, float)
tag_grad = mesh.createTag("grad_range", 9, float)

num_vert = mesh.getNumOfType( iBase.Type.vertex )
vertices = mesh.getEntities( iBase.Type.vertex )
coords = mesh.getVtxCoords(vertices)

tag_data = []
tag_grad_data = []

for i in xrange(num_vert):
    tag_data += [ [0.0, 0.0, 0.0] ]
    tag_grad_data += [ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ]

tag[vertices] = tag_data
tag_grad[vertices] = tag_grad_data

mesh.save("hex_range.vtk")

del mesh
del tag
del vertices
del coords
del tag_data

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

##---------------------------------------------------------------------------##
## Small linear tet mesh - tagged with a scalar range.
## from origin
mesh = iMesh.Mesh()
mesh.load("tet_mesh.vtk")

tag = mesh.createTag("range", 1, float)
tag_grad = mesh.createTag("grad_range", 3, float)

num_vert = mesh.getNumOfType( iBase.Type.vertex )
vertices = mesh.getEntities( iBase.Type.vertex )
coords = mesh.getVtxCoords(vertices)

tag_data = []
tag_grad_data = []

for i in xrange(num_vert):
    tag_data += [ 0.0 ];
    tag_grad_data += [ [0.0, 0.0, 0.0] ];

tag[vertices] = tag_data
tag_grad[vertices] = tag_grad_data

mesh.save("tet_range.vtk")

del mesh
del tag
del vertices
del coords
del tag_data
