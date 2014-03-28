from pyne.mesh import Mesh, IMeshTag
m = Mesh(structured=True, structured_coords=[[0,1,2],[0,1],[0,1]])
m.q = IMeshTag(1, float)
m.q_bias = IMeshTag(1, float)
m.q[:] = [1, 3]
m.q_bias[:] = [1, 1]

m.mesh.save("simpliest.h5m")
