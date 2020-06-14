# Shor's Algorithm Implementation

## Step1

- Pick any int a<N. See if a and N are relatively prime.
- if gcd(a, N) != 1, then that must be a common factor of N. No need to process further.
- if gcd(a, N) == 1, then go to Step2

## Step2

- Compute the period(r) of f(x) = a^x mod N.
- if r is odd, go to Step1
- if r is even and a^(r/2)+1=0 mod N, go to Step1
- otherwise go to step 3

## Step3

- We know that a^r-1 = k*N
- and we know that r is even
- the (a^(r/2)-1)*(a^(r/2)+1) = k*p*q
- then
	p = gcd(a^(r/2)-1, N)
	q = gcd(a^(r/2)+1, N)

## Execute

```bash
python HW1.py num_to_factorize
```

## Video Explanation

- [video](https://www.youtube.com/watch?v=fzNOtX6CLWk&feature=youtu.be)


# References
- [qiskit textbook](https://qiskit.org/textbook/ch-algorithms/shor.html)
- [example project](https://github.com/ttlion/ShorAlgQiskit/blob/master/Shor_Normal_QFT.py)
- [Period finding in Shor’s Algorithm](https://medium.com/@jonathan_hui/qc-period-finding-in-shors-algorithm-7eb0c22e8202)
- [Circuit for Shor’s algorithm using 2n+3 qubits](https://arxiv.org/pdf/quant-ph/0205095.pdf)
- [IBM explanation of Shor's Algorithm](https://quantum-computing.ibm.com/docs/guide/q-algos/shor-s-algorithm)
- [Ahmet Çevik Quantum Computation Course Lecture Notes](http://ahmetcevik.com/wp-content/uploads/2020/04/Quantum-computing-lecture-notes.pdf)

Muhammed Furkan YAĞBASAN
