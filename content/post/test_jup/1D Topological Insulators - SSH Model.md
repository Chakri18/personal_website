
We are interested to understand some non-trivial topological phases in 1 dimension system. Our model system is a tight binding model with two atoms per unit cell with $t_1$ and $t_2$ as hopping parameters and the on-site energy is assumed to be zero. The Hamiltonian is given by,

$$ \hat{\mathcal{H}} = \sum_n \left[ t_1|n_a\rangle \langle n_b| + t_2|n_b\rangle \langle (n+1)_a| + h.c. \right]$$
where, $t_1 - $ Intra cell hopping amplitude; $t_2 - $ Inter cell hopping amplitude

### Importing the required libraries


```python
import numpy as np
from bqplot import pyplot as plt
import ipywidgets as ipy
%matplotlib notebook
```

### Functions required for the task


```python
def SSH_model(N = 25, t1 = 2.0, t2 = 2.0):
    Hmat = np.zeros([2*N,2*N])
    for i in range(len(Hmat)-1):
        Hmat[i][i+1] = (1-i%2)*t1 + (i%2)*t2 #0.5*np.random.random()
        Hmat[i+1][i] = Hmat[i][i+1]
        #Hmat[i][i] = np.random.random()
    return Hmat

def compute(Hmat):
    eigenenergies, eigenstate = np.linalg.eigh(Hmat)
    return eigenenergies, eigenstate

def update(t1 = 2.0):
    N = 100; t2 = 2.0;
    ee,es = compute(SSH_model(N, t1,t2))
    fig1 = plt.figure(title="Eigenenergies", padding_y = 0, 
                      max_aspect_ratio =1, min_aspect_ratio =1)
    plt.scatter(range(2*N),ee, colors = ['red'])
    fig2 = plt.figure(title="Eigenstates", padding_y = 0)
    plt.plot(range(2*N),np.absolute(es[:,0]))
    plt.plot(range(2*N),np.absolute(es[:,N]))
    return ipy.HBox([fig1,fig2])
```

### Visualization of the Tight-Binding Matrix


```python
plt.figure(title="The tight-binding matrix", padding_y = 0, 
           max_aspect_ratio =1, min_aspect_ratio =1)
Hmat = SSH_model(100,1,2)
plt.gridheatmap(Hmat[:10,:10])
plt.show()
```

### Energy Spectrum and the Eigenstates


```python
ipy.interact(update, t1 = (0,4,0.05))
N = 15
ih = np.linspace(0,4,101); ee1 = np.zeros([len(ih),2*N])
for j in range(len(ih)):
    ee1[j],_ = compute(SSH_model(N,ih[j],2))
plt.figure(title="Energy Spectrum versus intra cell hopping(N = 15)")
for k in range(2*N):
    plt.plot(ih,ee1[:,k])
plt.ylabel('Energy')
plt.xlabel('Intra cell hopping(t1)')
plt.show()
N = 15
ih = np.linspace(0,4,101); ee1 = np.zeros([len(ih),2*N])
for j in range(len(ih)):
    ee1[j],_ = compute(SSH_model(N,2,ih[j]))
plt.figure(title="Energy Spectrum versus inter cell hopping(N = 15)")
for k in range(2*N):
    plt.plot(ih,ee1[:,k])
plt.ylabel('Energy')
plt.xlabel('Inter cell hopping(t2)')
plt.show()
```

### Let us add some random noises to our Hamiltonian

### 1. Noise term in the hopping terms


```python
def SSH_model(N = 25, t1 = 2.0, t2 = 2.0):
    Hmat = np.zeros([2*N,2*N])
    for i in range(len(Hmat)-1):
        Hmat[i][i+1] = (1-i%2)*t1 + (i%2)*t2 + np.random.random()
        Hmat[i+1][i] = Hmat[i][i+1]
        #Hmat[i][i] = np.random.random()
    return Hmat

def compute(Hmat):
    eigenenergies, eigenstate = np.linalg.eigh(Hmat)
    return eigenenergies, eigenstate
```


```python
N = 15
ih = np.linspace(0,4,101); ee1 = np.zeros([len(ih),2*N])
for j in range(len(ih)):
    ee1[j],_ = compute(SSH_model(N,ih[j],2))
plt.figure(title="Energy Spectrum versus intra cell hopping(N = 15)")
for k in range(2*N):
    plt.plot(ih,ee1[:,k])
plt.ylabel('Energy')
plt.xlabel('Intra cell hopping(t1)')
plt.show()
```

### 2. Noise term in the onsite term


```python
def SSH_model(N = 25, t1 = 2.0, t2 = 2.0):
    Hmat = np.zeros([2*N,2*N])
    for i in range(len(Hmat)-1):
        Hmat[i][i+1] = (1-i%2)*t1 + (i%2)*t2 
        Hmat[i+1][i] = Hmat[i][i+1]
        Hmat[i][i] = np.random.random()#(1-i%2)*2 + (i%2)*(-1*2) ; 
        Hmat[-1][-1] = np.random.random()
    return Hmat

def compute(Hmat):
    eigenenergies, eigenstate = np.linalg.eigh(Hmat)
    return eigenenergies, eigenstate
```


```python
N = 15
ih = np.linspace(0,4,101); ee1 = np.zeros([len(ih),2*N])
for j in range(len(ih)):
    ee1[j],_ = compute(SSH_model(N,ih[j],2))
plt.figure(title="Energy Spectrum versus intra cell hopping(N = 15)")
for k in range(2*N):
    plt.plot(ih,ee1[:,k])
plt.ylabel('Energy')
plt.xlabel('Intra cell hopping(t1)')
plt.show()
```
