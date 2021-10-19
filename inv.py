import numpy as np
from scipy.linalg import null_space
import igl
import polyscope as ps

v0, _, n0, f0, _, _ = igl.read_obj("./data/output.obj")

from pylab import cross,dot,inv

def rot(U,V):
    W=cross(U,V)
    A=np.array([U,W,cross(U,W)]).T
    B=np.array([V,W,cross(V,W)]).T
    return dot(B,inv(A))

def get_R(a):
  b = np.array([[0, 0, 1]])

  frame_a = np.concatenate([a, null_space(a).T], 0)
  frame_b = np.concatenate([b, null_space(b).T], 0)
  
  R = frame_b @ (np.linalg.inv(frame_a))

  print(R @ frame_a)
  print(R)
  print(frame_a)

  return R


n0 = np.zeros([f0.shape[0], 3])
n0 = igl.per_face_normals(v0, f0, n0)
vec_n0 = n0.mean(0)
vec_n0 /= np.linalg.norm(vec_n0)

# R = get_R(vec_n0.reshape(1, 3))
R = rot(vec_n0.reshape(-1), np.array([0, 0, 1]))

v05_norm = v0[5] + vec_n0 * 10
# v0 = np.vstack([v0, v05_norm])
# v0 = np.vstack([v0, np.array([0, 0, 0])])
# v0 = np.vstack([v0, vec_n0 * 10])

v0_ = R @ v0.T

e = np.array([[1629, 1630]])

v0_ = R @ v0.T
v0_ = v0_.T
v0_[:, 2] = 0
# v0_ = v0_[:1628, :]

igl.write_obj("rotated.obj", v0_, f0)

# ps.init()

# lc = ps.register_curve_network('c0', v0, e)
# pc = ps.register_point_cloud('p0', v0)

# lc = ps.register_curve_network('c1', v0_.T, e)
# pc = ps.register_point_cloud('p1', v0_.T)

# ps.show()