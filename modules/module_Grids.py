#================================================================================
#
#			       Module Grids
#			===========================
#
#	The various classes implemented in this module are used to compute
#	sums on the first Brillouin zone by creating grids and tesselations.
#================================================================================

from module_Constants import *


# utility function
def get_index(i,j):
	I = j*(j+1)/2+i
	return I


class Wedge:
	"""
	This is a container class which stores the 
	k points and their connection topology to form
	an irreducible wedge of the first Brillouin zone.
	"""

	def __init__(self,list_k,triangles_indices):
		self.list_k            = list_k
		self.triangles_indices = triangles_indices

		# create cross product for every triangle

		self.cross_products = []
		self.Jacobians      = []

		for t in self.triangles_indices:
			k1 = self.list_k[t[0]]
			k2 = self.list_k[t[1]]
			k3 = self.list_k[t[2]]

			self.cross_products.append( N.linalg.norm(N.cross(k2-k1,k3-k1)) )

                        J = complex(1.,0.)*N.array([ [k2[0]-k1[0],k3[0]-k1[0]] ,\
                                                     [k2[1]-k1[1],k3[1]-k1[1]] ])
		
	
                        self.Jacobians.append(J)


                self.cross_products = N.array(self.cross_products)
                self.Jacobians      = N.array(self.Jacobians)


class TesselationGrid:
	"""
	This class generates and stores a k-grid which is consistent
	with the hexagonal symmetry of graphene.
	"""

	def __init__(self,nx_gridsize):

		self.get_symmetries()
		irreducible_wedge = self.generate_irreducible_Wedge(nx_gridsize)
		self.generate_list_wedges()

	def generate_list_wedges(self):

		self.list_wedges = []

		self.irreducible_wedge.triangles_indices
		for R in self.list_D6:
			list_k = []
			for k in self.irreducible_wedge.list_k:
				list_k.append(N.dot(R,k))
			list_k = N.array(list_k)



                        # Make sure the Jacobian is always positive!
                        triangles_indices = deepcopy(self.irreducible_wedge.triangles_indices)
                        if N.linalg.det(R) < 0.:
                                triangles_indices[:,1] = self.irreducible_wedge.triangles_indices[:,2]
                                triangles_indices[:,2] = self.irreducible_wedge.triangles_indices[:,1]

			wedge = Wedge(list_k,triangles_indices)

			self.list_wedges.append(wedge)



	def generate_irreducible_Wedge(self,nx_gridsize):
		"""
		Generate the points inside a regular square
		"""

                delta = 0.001
		fy = (N.sqrt(3.)/3.-delta)*twopia # make the size just a tiny bit smaller, to avoid edge of 1BZ
		fx = fy/N.sqrt(3.)

		integers = N.arange(nx_gridsize+1)

		list_triangles = []

		ki = []
		kj = []

		triangles_indices  = []

		# create point array
		for j in integers:
        		for i in integers[:j+1]:
				k = get_index(i,j)
	
				ki.append(i)
				kj.append(j)

		list_k = N.array([N.array(ki)*fx/nx_gridsize, N.array(kj)*fy/nx_gridsize]).transpose()

                """
                tol = 1e-10
                list_norm_k = N.sqrt(N.sum(list_k**2,axis=1))
                I           = N.argwhere(N.abs(list_norm_k - norm_K_point) < tol)

                if len(I) != 0: 
                        iK = I[0]
                        I = N.argsort(N.abs(list_norm_k - norm_K_point))

                        # Move the K point just a little bit off the actual K point

                        K0    = list_k[iK]
                        k_nn1 = list_k[I[1]]
                        k_nn2 = list_k[I[2]]

                        K_avg    = k_nn1 + 0.999*(K0-k_nn1)

                        list_k[iK] = K_avg
                """


		# create connections
		for i in integers[:-1]:

			j  = i
			I1 = get_index(i,j+1)
			I2 = get_index(i,j)
			I3 = get_index(i+1,j+1)
			triangles_indices.append([I1,I2,I3])

        		for j in integers[i+1:-1]:

				I1 = get_index(i,j+1)
				I2 = get_index(i,j)
				I3 = get_index(i+1,j+1)
				triangles_indices.append([I1,I2,I3])

				I1 = get_index(i+1,j)
				I2 = get_index(i+1,j+1)
				I3 = get_index(i,j)

				triangles_indices.append([I1,I2,I3])



		triangles_indices = N.array(triangles_indices )
		self.irreducible_wedge = Wedge(list_k=list_k,triangles_indices=triangles_indices)

		return 


	def get_symmetries(self):
		theta = N.pi/3.
		# Identity
		E =  N.array([[  1.,  0.],
			      [  0.,  1.]])


		# pi/3 rotation left
		C6_1  = N.array([[ N.cos(theta), N.sin(theta)],
				 [-N.sin(theta), N.cos(theta)] ])

		# pi/3 rotation right
		C6_2  = N.array([[ N.cos(theta),-N.sin(theta)],
				 [ N.sin(theta), N.cos(theta)] ])

		# 2pi/3 rotation left
		C3_1  = N.dot(C6_1,C6_1)

		# 2pi/3 rotation right
		C3_2  = N.dot(C6_2,C6_2)

		# pi rotation 
		C2    = N.dot(C6_1,C3_1)

		# vertical mirror plane, and rotations by 2pi/3
		# these mirror planes cross the face of the hexagon
		C2p_1  =  N.array([[ -1.,  0.],
				   [  0.,  1.]])
		C2p_2  =  N.dot(C3_1,C2p_1)
		C2p_3  =  N.dot(C3_2,C2p_1)

		# mirror planes
		# these mirror planes cross the corners of the hexagon
		C2pp_1 =  N.dot(C6_1,C2p_1)
		C2pp_2 =  N.dot(C3_1,C2pp_1)
		C2pp_3 =  N.dot(C3_2,C2pp_1)


		self.list_D6 = N.array([ E, C6_1, C6_2, C3_1, C3_2, C2, C2p_1, C2p_2, C2p_3, C2pp_1, C2pp_2, C2pp_3])


