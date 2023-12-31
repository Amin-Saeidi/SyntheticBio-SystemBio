import numpy as np    
class ForwardEuler:
    def __init__(self, f):
        # Wrap user’s f in a new function that always
        # converts list/tuple to array (or let array be array)
        self.f = lambda u, t: np.asarray(f(u, t), float)
    def set_initial_condition(self, U0):
        if isinstance(U0, (float,int)): # scalar ODE
            self.neq = 1 # no of equations
            U0 = float(U0)
        else: # system of ODEs
            U0 = np.asarray(U0)
            self.neq = U0.size # no of equations
        self.U0 = U0
    def solve(self, time_points):
        self.t = np.asarray(time_points)
        N = len(self.t)
        if self.neq == 1: # scalar ODEs
            self.u = np.zeros(N)
        else: # systems of ODEs
            self.u = np.zeros((N,self.neq))
        # Assume that self.t[0] corresponds to self.U0
        self.u[0] = self.U0
        # Time loop
        for n in range(N-1):
            self.n = n
            self.u[n+1] = self.advance()
        return self.u, self.t
    def advance(self):
        # Create local variables to get rid of "self." in
        # the numerical formula
        u, f, n, t = self.u, self.f, self.n, self.t
        #dt is not necessarily constant:
        dt = t[n+1]-t[n]
        unew = u[n] + dt*f(u[n], t[n])
        return unew