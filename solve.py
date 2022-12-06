import sympy

def solve():
	a = 0.2
	m = 1.4
	n = 3.4
	deltaQ = 3.0
	ls = 10.0
	lm = 20.0


	p = [100.0]
	q = []
	pm = [300.0]
	ps = [175.0]
	c = []
	deltaC = []

	def solveQ(i):
		qi = (-1.0 * (m**n)*((p[i] + a * pm[i])**(n-1))* (p[i] - pm[i] - ps[i])) ** (1.0/n)
		q.append(qi)

	def solveDeltaC(i):
		A = n * (q[i] ** (n-1))
		B = (m ** n) * (p[i] + a * pm[i]) ** (n-2) * (n * p[i] - (n-1-a) * pm[i] - (n-1) * ps[i])
		C = (m ** n) * ((p[i] + a * pm[i]) ** (n-2)) * ((a * (n-1) - 1) * p[i] - n * a * pm[i] - a * (n-1) * ps[i])
		D = -1.0 * (m ** n) * ((p[i] + a * pm[i]) ** (n-1))

		R2 = q[0] * q[0] + p[0] * p[0] 
		deltaP = -1.0 * q[i] * ((R2 - q[i] * q[i]) ** 0.5) * deltaQ
		
		deltaCi = (A * deltaQ + B * deltaP) / (C * lm * pm[i] + D * ls * ps[i])
		deltaC.append(deltaCi)

	def solveC(i):
		if i == 0:
			c.append(deltaC[0])
		else:
			c.append(c[i-1] + deltaC[i])

	def solvePMPS(i):
		pmi = pm[i-1] - lm * pm[i-1]*deltaC[i-1]
		psi = ps[i-1] - ls * ps[i-1]*deltaC[i-1]
		pm.append(pmi)
		ps.append(psi)
		if pm[i] < 0:
			pm[i] = 0.0

	def solvePQ(i):
		print("solving PQ", i)
		R2 = q[0] * q[0] + p[0] * p[0] 
		pi, qi = sympy.symbols('pi, qi')
		f1 = qi ** 2 + pi ** 2 - R2
		f2 = qi ** 2 + (m**n) * (pi + a * pm[i]) ** (n-1) * (pi - pm[i] - ps[i])
		solution = sympy.nsolve((f1, f2), (pi, qi), (1, 1))
		sol = list(solution)
		pii = sol[0]
		qii = sol[1]
		p.append(pii)
		q.append(qii)

	def printAll():
		print('p:')
		print(p)
		print('q:')
		print(q)
		print('end')
			
	MAX = 100
	for i in range(0, 100):
		solveQ(i)
		solveDeltaC(i)
		solveC(i)
		solvePMPS(i + 1)
		solvePQ(i)

		if q[i] * 1.0 / p[i] < m or c[i] > 0.2 or pm[i] <= 0:
			break
	
	printAll()

if __name__ == '__main__':
	solve()