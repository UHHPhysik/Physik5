import sys, json, math


_EPSILON = 1e-12


def _signed_sqrt(value):
	result = abs(value)**0.5
	if result < _EPSILON:
		return 0.
	return math.copysign(result, value)


def _calc_eta(pz, pt):
	if pt < _EPSILON:
		return math.copysign(float('inf'), pz)
	return math.asinh(pz / pt)


def _calc_y(pz, E):
	arg = min(1, abs(pz / E))
	if arg == 1:
		return math.copysign(float('inf'), pz)
	return math.atanh(pz / E)


class LorentzVector(object):
	def __init__(self, pt=None, eta=None, phi=None, m=None, E=None, px=None, py=None, pz=None, q=0):
		self.q = q
		if E is None:
			self.E = (pt**2 * math.cosh(eta)**2 + m**2)**0.5
		else:
			self.E = float(E)
		if px is None:
			self.px = pt * math.cos(phi)
		else:
			self.px = float(px)
		if py is None:
			self.py = pt * math.sin(phi)
		else:
			self.py = float(py)
		if pz is None:
			self.pz = pt * math.sinh(eta)
		else:
			self.pz = float(pz)
		if abs(self.px / self.E) < _EPSILON:
			self.px = 0
		if abs(self.py / self.E) < _EPSILON:
			self.py = 0
		if abs(self.pz / self.E) < _EPSILON:
			self.pz = 0

	def __repr__(self):
		return '(pt=%.2f, eta=%.2f, phi=%.2f, m=%.2f, q=%d)' % (self.pt, self.eta, self.phi, self.m, self.q)

	def __str__(self):
		return '(E=%.2f, px=%.2f, py=%.2f, pz=%.2f, m=%.2f)' % (self.E, self.px, self.py, self.pz, self.m)

	beta = property(lambda self: self.p / self.E)
	charge = property(lambda self: self.q)
	eta = property(lambda self: _calc_eta(self.pz, self.pt))
	Et = property(lambda self: self.E * self.pt / max(_EPSILON, self.p))
	gamma = property(lambda self: self.E / self.m)
	m = property(lambda self: _signed_sqrt(self.E**2 - self.p2))
	mt = property(lambda self: _signed_sqrt(self.E**2 - self.pz**2))
	p2 = property(lambda self: (self.px**2 + self.py**2 + self.pz**2))
	phi = property(lambda self: math.atan2(self.py, self.px))
	p = property(lambda self: _signed_sqrt(self.p2))
	pt = property(lambda self: _signed_sqrt(self.px**2 + self.py**2))
	theta = property(lambda self: math.atan2(self.pt, self.pz))
	y = property(lambda self: _calc_y(self.pz, self.E))

	def __add__(self, other):
		return LorentzVector(q=self.q + other.q, E=self.E + other.E,
			px=self.px + other.px, py=self.py + other.py, pz=self.pz + other.pz)

	def __mul__(self, other):
		if isinstance(other, LorentzVector):
			return (self.E * other.E - self.px * other.px - self.py * other.py - self.pz * other.pz)
		return LorentzVector(E=self.E * other, px=self.px * other, py=self.py * other, pz=self.pz * other, q=self.q)

	__rmul__ = __mul__

	def __sub__(self, other):
		return LorentzVector(q=self.q - other.q, E=self.E - other.E,
			px=self.px - other.px, py=self.py - other.py, pz=self.pz - other.pz)



def read_events(fn_list, hlt_name = 'HLT_Mu30'):
	for fn in fn_list:
		with open(fn) as fp:
			sys.stdout.write('Loading ... %s\n' % fn)
			collision_data = json.load(fp)
			sys.stdout.write('Available Trigger: %s\n' % str.join(', ', collision_data['metadata']['trigger']))
			hlt_mask = (1 << collision_data['metadata']['trigger'].index(hlt_name))
			for event in collision_data['events']:
				(hlt_bits, muon_list) = event
				if hlt_bits & hlt_mask:
					yield [LorentzVector(pt=pt, eta=eta, phi=phi, m=0.105658, charge=charge) for (charge,pt,eta,phi) in muon_list]
	sys.stdout.write('Finished processing!')

	
if __name__ == '__main__':
	import matplotlib.pyplot as plt
	data_pt = []
	for muon_list in read_events(['mu_data_nmin2_nmax4_part1.json'], 'HLT_Mu30'):
		data_pt.append(muon_list[0].pt)
	plt.hist(data_pt, bins=50, range=(0, 170), log=True)
	plt.savefig('mu_spec.png')
	plt.show()
