#########################################
#										#
#	CENG798 - Quantum Computation		#
#										#
#	Shor's Algorithm Implementation		#
#										#
#	Muhammed Furkan YAĞBASAN			#
#	2099505								#
#										#
#########################################

import random
import sys
import numpy as np
from math import pi, ceil, log, floor, sqrt
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, execute, Aer

#############################################################################

def swap_registers(circuit, reg, n):
	for i in range(n//2):
		circuit.swap(reg[i], reg[n-i-1])

# applies qft to first n qubits of reg
def qft(circuit, reg, n, swaps):
	for j in range(n):
		circuit.h(reg[n-j-1])
		for m in range(n-j-1):
			circuit.cu1(pi/float(2**(n-j-1-m)), reg[m], reg[n-j-1])

	if(swaps):
		swap_registers(circuit, reg, n)

# applies inverse qft to first n qubits of reg
def qft_dagger(circuit, reg, n, swaps):
	if(swaps):
		swap_registers(circuit, reg, n)

	for j in range(n):
		for m in range(j):
			circuit.cu1(-pi/float(2**(j-m)), reg[m], reg[j])
		circuit.h(reg[j])

##############################################################################################

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

"""Function that calculates the array of angles to be used in the addition in Fourier Space"""
def getAngles(a,N):
    s=bin(int(a))[2:].zfill(N) 
    angles=np.zeros([N])
    for i in range(0, N):
        for j in range(i,N):
            if s[j]=='1':
                angles[N-i-1]+=pow(2, -(j-i))
        angles[N-i-1]*=np.pi
    return angles

"""Creation of a doubly controlled phase gate"""
def ccphase(circuit,angle,ctl1,ctl2,tgt):
    circuit.cu1(angle/2,ctl1,tgt)
    circuit.cx(ctl2,ctl1)
    circuit.cu1(-angle/2,ctl1,tgt)
    circuit.cx(ctl2,ctl1)
    circuit.cu1(angle/2,ctl2,tgt)

"""Creation of the circuit that performs addition by a in Fourier Space"""
"""Can also be used for subtraction by setting the parameter inv to a value different from 0"""
def phiADD(circuit,q,a,N,inv):
    angle=getAngles(a,N)
    for i in range(0,N):
        if inv==0:
            circuit.u1(angle[i],q[i])
        else:
            circuit.u1(-angle[i],q[i])

"""Single controlled version of the phiADD circuit"""
def cphiADD(circuit,q,ctl,a,n,inv):
    angle=getAngles(a,n)
    for i in range(0,n):
        if inv==0:
            circuit.cu1(angle[i],ctl,q[i])
        else:
            circuit.cu1(-angle[i],ctl,q[i])

"""Doubly controlled version of the phiADD circuit"""      
def ccphiADD(circuit,q,ctl1,ctl2,a,n,inv):
    angle=getAngles(a,n)
    for i in range(0,n):
        if inv==0:
            ccphase(circuit,angle[i],ctl1,ctl2,q[i])
        else:
            ccphase(circuit,-angle[i],ctl1,ctl2,q[i])
        
"""Circuit that implements doubly controlled modular addition by a"""
def ccphiADDmodN(circuit, q, ctl1, ctl2, aux, a, N, n):
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)
    phiADD(circuit, q, N, n, 1)
    qft_dagger(circuit, q, n, 0)
    circuit.cx(q[n-1],aux)
    qft(circuit,q,n,0)
    cphiADD(circuit, q, aux, N, n, 0)
    
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)
    qft_dagger(circuit, q, n, 0)
    circuit.x(q[n-1])
    circuit.cx(q[n-1], aux)
    circuit.x(q[n-1])
    qft(circuit,q,n,0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)

"""Circuit that implements the inverse of doubly controlled modular addition by a"""
def ccphiADDmodN_inv(circuit, q, ctl1, ctl2, aux, a, N, n):
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)
    qft_dagger(circuit, q, n, 0)
    circuit.x(q[n-1])
    circuit.cx(q[n-1],aux)
    circuit.x(q[n-1])
    qft(circuit, q, n, 0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)
    cphiADD(circuit, q, aux, N, n, 1)
    qft_dagger(circuit, q, n, 0)
    circuit.cx(q[n-1], aux)
    qft(circuit, q, n, 0)
    phiADD(circuit, q, N, n, 0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)

"""Circuit that implements single controlled modular multiplication by a"""
def cMULTmodN(circuit, ctl, q, aux, a, N, n):
    qft(circuit,aux,n+1,0)
    for i in range(0, n):
        ccphiADDmodN(circuit, aux, q[i], ctl, aux[n+1], (2**i)*a % N, N, n+1)
    qft_dagger(circuit, aux, n+1, 0)

    for i in range(0, n):
        circuit.cswap(ctl,q[i],aux[i])

    a_inv = modinv(a, N)
    qft(circuit, aux, n+1, 0)
    i = n-1
    while i >= 0:
        ccphiADDmodN_inv(circuit, aux, q[i], ctl, aux[n+1], pow(2,i)*a_inv % N, N, n+1)
        i -= 1
    qft_dagger(circuit, aux, n+1, 0)


# calculate period of f(x)=a^x mod N
def periodFinding(N, a):
	n = ceil(log(N, 2))
	if(n<4):
		n=4
	m = 2*n

	qregUp 		= QuantumRegister(m)
	qregDown 	= QuantumRegister(n)
	aux 		= QuantumRegister(n+2)
	creg 		= ClassicalRegister(m)

	qc = QuantumCircuit(qregUp, qregDown, aux, creg)

	qft(qc, qregUp, m, 1)
	qc.x(qregDown[0])

	for i in range(0, m):
		cMULTmodN(qc, qregUp[i], qregDown, aux, int(pow(a, pow(2, i))), N, n)

	qft_dagger(qc, qregUp, m, 1)
	qc.measure(qregUp, creg)


	shots = 1000
	job = execute(qc, Aer.get_backend('qasm_simulator'), shots=shots)
	counts = job.result().get_counts(qc)

	### process data start
	i=0
	avg = 0
	while i < len(counts):
		val = list(counts.values())[i]
		output_desired = list(counts.keys())[i]
		x_value = int(output_desired, 2)
		print("value:{0} \t times:{1}".format(x_value, val))

		avg += (val*val)/shots

		i += 1

	i=0
	Variance = 0
	while i < len(counts):
		val = list(counts.values())[i]

		Variance += (((val-avg)**2)*val)/shots

		i += 1

	stdDev = sqrt(Variance)
	treshold = avg - stdDev

	if(stdDev<avg*0.3):
		treshold = 0

	print("---------------------------")
	print("avg     : {0}".format(avg))
	print("stdDev  : {0}".format(stdDev))
	print("treshold: {0}".format(treshold))
	print("---------------------------")
	### process data end

	i=0
	MOverR = 2**m
	while i < len(counts):
		val = list(counts.values())[i]
		if(val<treshold):
			i+=1
			continue
		output_desired = list(counts.keys())[i]
		x_value = int(output_desired, 2)

		print("value:{0} \t times:{1}".format(x_value, val))

		if(x_value<MOverR and x_value!=0):
			MOverR = x_value
			
		i+=1;

	print("---------------------------")
	return floor(((2**m)/MOverR)+0.5) #contFrac(MOverR, m)

def gcd(a, b):
	if(b==0):
		return a
	else:
		return gcd(b, a%b)

# select a random number a
def Step1(N):
	a = random.randint(2, N-1)
	print("Selected a: {0}".format(a))

	gcdAN = gcd(a,N)
	if(gcdAN!=1):
		p = gcdAN
		q = N/gcdAN
		return (p, q)
	else:
		return (a, 0)

# find the period and check if it satisfies the conditions 
def Step2(N, a):
	print("executing the quantum circuit...")
	r = periodFinding(N, a)
	print("found r:{0}".format(r))
	if(r%2==1):
		return 0
	else:
		if((a**(r/2)+1)%N==0):
			print("a**(r/2)+1) mod N = 0")
			return 0
		else:
			return r


def MyShor(N):
	r = 0
	a = 0
	p = 0
	q = 0
	while(True):
		print("-----------Step1--------------")
		(p,q) = Step1(N)
		if(q!=0):
			print("Lucky selection for a")
			return (p,q)
		else:
			a = p
			print("-----------Step2--------------")
			r = Step2(N,a)
			if(r!=0):
				p = gcd((a**(r/2)-1),N)
				q = gcd((a**(r/2)+1),N)
				break

			print("\n+++++++++++++++++++++++++++++++++++++++++++++\n")

	return (p,q)


def main():
	if(len(sys.argv) != 2):
		print("usage: python HW1.py num_to_factor")
	else:
		factors = MyShor(int(sys.argv[1]))
		print("factors: {0} and {1}".format(factors[0], factors[1]))

if __name__ == "__main__":
	main()


