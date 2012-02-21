# pytaps script to tag the fine mesh with data and setup a tag without data 
# on the coarse mesh to interpolate onto.

import math
from itaps import iBase, iMesh

##---------------------------------------------------------------------------##
## Quad flat surface - tagged with a scalar domain of coordinate 
## distance from origin
mesh = iMesh.Mesh()
mesh.load("quad_flat_surf.vtk")

tag = mesh.createTag("domain", 1, float)

num_vert = mesh.getNumOfType( iBase.Type.vertex )
vertices = mesh.getEntities( iBase.Type.vertex )
coords = mesh.getVtxCoords(vertices)

tag_data = []

for i in xrange(num_vert):
    tag_data += [ coords[i][0]*coords[i][0] + \
                  coords[i][1]*coords[i][1] + \
                  coords[i][2]*coords[i][2] ]

tag[vertices] = tag_data

mesh.save("tagged_quad_flat_surf.vtk")

del mesh
del tag
del vertices
del coords
del tag_data

##---------------------------------------------------------------------------##
## Delaunay flat surface - tagged with a scalar domain of coordinate 
## distance from origin
mesh = iMesh.Mesh()
mesh.load("delaunay_flat_surf.vtk")

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

mesh.save("tagged_delaunay_flat_surf.vtk")

del mesh
del tag
del vertices
del coords
del tag_data

##---------------------------------------------------------------------------##
## Quad curved surface - tagged with a scalar domain of coordinate 
## distance from origin
mesh = iMesh.Mesh()
mesh.load("quad_curved_surf.vtk")

tag = mesh.createTag("domain", 1, float)

num_vert = mesh.getNumOfType( iBase.Type.vertex )
vertices = mesh.getEntities( iBase.Type.vertex )
coords = mesh.getVtxCoords(vertices)

tag_data = []

for i in xrange(num_vert):
    tag_data += [ coords[i][0]*coords[i][0] + \
                  coords[i][1]*coords[i][1] + \
                  coords[i][2]*coords[i][2] ]

tag[vertices] = tag_data

mesh.save("tagged_quad_curved_surf.vtk")

del mesh
del tag
del vertices
del coords
del tag_data

##---------------------------------------------------------------------------##
## Delaunay curved surface - tagged with a scalar domain of coordinate 
## distance from origin
mesh = iMesh.Mesh()
mesh.load("delaunay_curved_surf.vtk")

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

mesh.save("tagged_delaunay_curved_surf.vtk")

del mesh
del tag
del vertices
del coords
del tag_data
