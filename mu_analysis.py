import sys, json, math


class LorentzVector(object):
	def __init__(self, pt=None, eta=None, phi=None, m=None, E=None, px=None, py=None, pz=None, charge=0):
		(self.E, self.px, self.py, self.pz, self.charge) = (E, px, py, pz, charge)
		if E is None:
			self.E = (pt**2 * math.cosh(eta)**2 + m**2)**0.5
		if px is None:
			self.px = pt * math.cos(phi)
		if py is None:
			self.py = pt * math.sin(phi)
		if pz is None:
			self.pz = pt * math.sinh(eta)

	def _calc_m(self):
		m2 = self.E**2 - self.p2
		if m2 >= 0:
			return m2**0.5
		return - (-m2)**0.5

	m = property(_calc_m)
	y = property(lambda self: 0.5*math.log((self.E + self.pz)/(self.E - self.pz)))
	pt = property(lambda self: (self.px**2 + self.py**2)**0.5)
	p = property(lambda self: self.p2**0.5)
	p2 = property(lambda self: (self.px**2 + self.py**2 + self.pz**2))
	eta = property(lambda self: math.asinh(self.pz / self.pt))
	phi = property(lambda self: math.atan(self.py / self.px))

	def __add__(self, other):
		return LorentzVector(charge=self.charge + other.charge, E=self.E + other.E,
			px=self.px + other.px, py=self.py + other.py, pz=self.pz + other.pz)

	def __sub__(self, other):
		return LorentzVector(charge=self.charge - other.charge, E=self.E - other.E,
			px=self.px - other.px, py=self.py - other.py, pz=self.pz - other.pz)

	def __mul__(self, other):
		if isinstance(other, LorentzVector):
			return (self.E * other.E - self.px * other.px - self.py * other.py - self.pz * other.pz)
		return LorentzVector(E=self.E * other, px=self.px * other, py=self.py * other, pz=self.pz * other, charge=self.charge)

	__rmul__ = __mul__

	def __repr__(self):
		try:
			return '(pt=%.2f, eta=%.2f, phi=%.2f, m=%.2f, Q=%d)' % (self.pt, self.eta, self.phi, self.m, self.charge)
		except:
			return '(pt=%.2f, eta=%.2f, phi=%.2f, m=%.2f, Q=%d)' % (self.pt, self.eta, self.phi, self.m, self.charge)


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
