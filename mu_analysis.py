import sys, json, math


class LorentzVector(object):
	def __init__(self, pt=None, eta=None, phi=None, m=None, E=None, px=None, py=None, pz=None, charge=0):
		(self._E, self._px, self._py, self._pz, self.charge) = (E, px, py, pz, charge)
		if E is None:
			self._E = (pt**2 * math.cosh(eta)**2 + m**2)**0.5
		if px is None:
			self._px = pt * math.cos(phi)
		if py is None:
			self._py = pt * math.sin(phi)
		if pz is None:
			self._pz = pt * math.sinh(eta)
	m = property(lambda self: (self._E**2 - self.p2)**0.5)
	pt = property(lambda self: (self._px**2 + self._py**2)**0.5)
	p2 = property(lambda self: (self._px**2 + self._py**2 + self._pz**2))
	eta = property(lambda self: math.asinh(self._pz / self.pt))
	phi = property(lambda self: math.atan(self._py / self._px))

	def __add__(self, other):
		return LorentzVector(charge=self.charge + other.charge, E=self._E + other._E,
			px=self._px + other._px, py=self._py + other._py, pz=self._pz + other._pz)

	def __mul__(self, other):
		if isinstance(other, LorentzVector):
			return (self._E * other._E - self._px * other._px - self._py * other._py - self._pz * other._pz)
		return LorentzVector(E=self._E * other, px=self._px * other, py=self._py * other, pz=self._pz * other, charge=self.charge)

	def __repr__(self):
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


if __name__ == '__main__':
	import matplotlib.pyplot as plt
	data_pt = []
	for muon_list in read_events(['mu_data_nmin2_nmax4_part1.json'], 'HLT_Mu30'):
		data_pt.append(muon_list[0].pt)
	plt.hist(data_pt, bins=50, range=(0, 170), log=True)
	plt.savefig('mu_spec.png')
	plt.show()
