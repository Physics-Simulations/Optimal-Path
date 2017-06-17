from __future__ import division
import numpy as np

TOLERANCE = 5e-7
EPSILON = 1e-6

h = 0.1
DIRECTIONS = [[h, 0, 0.0], [-h, 0, 0.0], [0, h, 0.0], [0, -h, 0.0]]

ALPHA = 1

def hash_key(r):
	return '%.4f;%.4f' % (r[0], r[1])

def dehash(r):
	s = r.split(';')
	return float(s[0]), float(s[1])

def compare(u, v, h):
	return all(np.abs(u - v) < h)

def save_path(path, U, file_name):
	with open(file_name, 'w') as f:
		for r in path:
			f.write("%.8f %.8f %.8f\n" % (r[0], r[1], U(r)))

def U1(r):
	return np.sin(r[0]) + np.sin(r[1])

def U2(r):
	return np.sinh(r[0]) + np.sinh(r[1])

def U3(r):
	return -10 / np.sqrt(r[0]*r[0]+r[1]*r[1])

def U4(r):
	return 0.1 * (r[0]*r[0]+r[1]*r[1]) / 2

def dist(p, q, k):
	return np.power(p[0]-q[0],2) + np.power(p[1]-q[1],2) + k * np.power(p[2]-q[2], 2)

def weight_v1(r, dr):
	return dr[2]

def weight_v2(r, dr):
	return dr[2] - r[2]

def calc_energy(path):
	return sum(abs(P(path[k]) - P(path[k - 1])) for k in xrange(1, len(path))) * h

	return path

def method_1(r0, rf, U, weight, k):
	r = r0
	d = dist(r, rf, k)

	path = [r]

	inf = float('inf')

	while not compare(r, rf, h):
		next_r, min_U = None, inf
		for direc in DIRECTIONS:

			dr = r + direc
			dr[2] = U(dr)

			Udr = weight(r, dr)

			if dist(dr, rf, k) <= d and Udr < min_U:
				next_r = dr
				min_U = Udr

		if next_r is None:
			raise ValueError("No path was found")

		r = next_r
		d = dist(r, rf, k)
		path.append(r)

	return path

def method_2(r0, rf, U, weight, k):
	# dictionary of precessors
	# each item contains a list of all the points that converg to it
	G = { }
	# dictionary of distances (energy)
	D = { } 

	D[hash_key(r0)] = 0

	def _relax(u, v, w):
		if D[v] > D[u] + w:
			D[v] = D[u] + w
			G[v] = u

	def loop(r):
		hr = hash_key(r)
		d = dist(r, rf, k)
		if compare(r, rf, h):
			return

		for direc in DIRECTIONS:
			dr = r + direc

			hdr = hash_key(dr)
			dr[2] = U(dr)

			if dist(dr, rf, k) <= d:
				if hdr not in G:
					G[hdr] = hr
					D[hdr] = D[hr] + weight(r, dr)
					loop(dr)
				_relax(hr, hdr, weight(r, dr))

	def find_path(start, end):
		"""
		Finds the shortest path between start and end.

		Returns the path with the energy spend.
		"""

		if end not in G:
			raise ValueError("No path was found")

		path = []
		while True:
			path.append(dehash(end))
			if end in start:
				break
			end = G[end]

		return path[::-1]

	loop(r0)
	return find_path(hash_key(r0), hash_key(rf))

def bellman_ford(r0, rf, weight, U, k):
	def build_graph():
		"""
		Builds the graph with the vertex that verify the condition
				d(a, rf) < d(b, rf)
		where a is one of the four neighbours of b
		"""
		G = { }

		def loop(r):
			hr = hash_key(r)
			if hr in G:
				return
			else:
				G[hr] = { }

			d = dist(r, rf, k)
			for direc in DIRECTIONS:
				dr = r + direc
				dr[2] = U(dr)

				dist_dr = dist(dr, rf, k)

				if dist_dr <= d:
					hdr = hash_key(dr)
					G[hr][hdr] = weight(r, dr)
					loop(dr)
		loop(r0)
		return G

	def algo(graph, start, end):
		"""
		Implementation of Bellman-Ford's shortest paths algorithm.
		It can be used for graphs with negative weight edges and if
		there is a neg-cycle it will find it.

		Complexity: O(V*E)
				given that E ~ O(V^2) this runs in cubic time
		"""
		def _initialise(graph, source):
			inf = float('inf')

			D = { } # dictionary of distances
			P = { } # dictionary of predecessors

			for u in graph:
				D[u] = inf
				P[u] = None
			D[source] = 0
			return D, P

		def _relax(graph, u, v, D, P):
			if D[v] > D[u] + graph[u][v]:
				D[v] = D[u] + graph[u][v]
				P[v] = u

		D, P = _initialise(graph, start)
		for _ in xrange(len(graph)):
			for u in graph:
				for v in graph[u]:
					_relax(graph, u, v, D, P)

		for u in graph:
			for v in graph[u]:
				if D[v] > D[u] + graph[u][v]:
					raise "Negative weight cycle found"

		return (D, P)

	def find_path(D, P, start, end):
		"""
		Finds the shortest path between start and end.

		returns the path with the weight.
		"""
		weight = D[end]
		path = []
		while True:
			path.append(dehash(end))
			if end == start:
				break
			end = P[end]

		return path[::-1], weight

	start, end = hash_key(r0), hash_key(rf)
	G = build_graph()
	D, P = algo(G, start, end)

	return find_path(D, P, start, end)


def dijkstra(r0, rf, weight, U, k):

	def build_graph():
		"""
		Builds the graph with the vertex that verify the condition
				d(a, rf) < d(b, rf)
		where a is one of the four neighbours of b
		"""
		G = { }
		def loop(r):
			hr = hash_key(r)
			if hr in G:
				return
			else:
				G[hr] = { }

			d = dist(r, rf, k)
			for direc in DIRECTIONS:
				dr = r + direc
				dr[2] = U(dr)

				dist_dr = dist(dr, rf, k)

				if dist_dr <= d:
					hdr = hash_key(dr)
					G[hr][hdr] = abs(weight(r, dr))
					loop(dr)

		loop(r0)
		return G

	def algo(graph, start, end=None):
		"""
		Implementation of Dijkstra shortest paths algorithm. This will find
		the shortest path to all elements from start (if no end is specified).

		In this case each vertex in the graph has a NON-NEGATIVE weight
		associated with it, represented in the graph as a dictionary in which
		each element u: graph[u] = {(vertex, weight) for each neighbourg of u}

		We here use a priority queue (binary heap) to find the nearest
		element to the start one, this is the same as doing a min operation
		over all the elements in the adjacency dictionary, and takes O(lg V) time.

		Complexity: O(V*lg(V) + E*lg(V))
		This can be improved by using a fibonacci heap instead of a binary heap.
		"""
		D, P = {}, {}
		Q = PriorityQueue() # estimated distances
		Q[start] = 0

		for u in Q:
			# finds de nearest element
			D[u] = Q[u]
			if u == end:
				break

			for v in graph[u]: # relax distances from u
				#relax(graph, u, v, D, P)
				uv_length = D[u] + graph[u][v]
				if v in D:
					if uv_length < D[v]:
						raise "Already found a better path to go."
				elif v not in Q or uv_length < Q[v]:
					Q[v] = uv_length
					P[v] = u

		return (D, P)

	def find_path(D, P, start, end):
		"""
		Finds the shortest path between start and end.

		returns the path with the weight.
		"""
		weight = D[end]
		path = []
		while True:
			path.append(dehash(end))
			if end == start:
				break
			end = P[end]
		return path[::-1], weight

	start, end = hash_key(r0), hash_key(rf)
	G = build_graph()
	D, P = algo(G, start, end)

	return find_path(D, P, start, end)