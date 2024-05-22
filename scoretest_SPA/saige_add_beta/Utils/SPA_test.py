import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.optimize import root_scalar
from scipy import stats

# Setting the seed
np.random.seed(1234)

class SPA():
    def __init__(self, G, mu):
        self.G = G
        self.mu = mu
        self.n = mu.shape[0] 

    def K(self, t):
        tem1 = np.sum(np.log(np.ones(self.n) - self.mu + self.mu * np.exp(self.G * np.ones(self.n) * t)))
        tem2 = t * np.sum(self.G * self.mu)
        return tem1 - tem2
    
    def K1(self, t):
        tem0 = (np.ones(self.n) - self.mu) * np.exp( - self.G * np.ones(self.n) * t) + self.mu
        tem1 = np.sum(self.G * self.mu / tem0)
        tem2 = np.sum(self.G * self.mu)
        return tem1 - tem2
    
    def K2(self, t):
        tem1 = (np.ones(self.n) - self.mu) * self.mu * np.power(self.G, 2) * np.exp( - self.G * np.ones(self.n) * t)
        tem2 = np.power((np.ones(self.n) - self.mu) * np.exp( - self.G * np.ones(self.n) * t) + self.mu, 2)
        return np.sum(tem1 / tem2)
    
    def saddle_assumption_func(self, t):
         return self.K1(t) - self.x

    def solve_t(self, x):
        '''solve the assumption to get t: K'(t) =  x_t'''
        self.x = x
        sol = root_scalar(self.saddle_assumption_func, bracket=[-100, 100], method='brentq')
        print(sol.root)
        return sol.root

    def calculate_F(self, x):
        t_hat = self.solve_t(x)
        w = np.sign(t_hat) * np.sqrt(2 * (t_hat * x - self.K(t_hat)))
        v = t_hat * np.sqrt(self.K2(t_hat))
        input = w + 1/w * np.log(v/w)
        P = stats.norm.pdf(loc = 0, scale = 1, x = input)
        # density_value = (1 / np.sqrt(2 * np.pi * self.K2(t))) * np.exp(self.K(t) - t * x)
        return P

if __name__ == "__main__":

    spa_stdnormal = SPA()


    # Bootstrapping
    # boots = np.array([np.mean(np.random.choice(x, size=len(x), replace=True)) for _ in range(100)])

    # Generating random exponential data
    boots = np.random.normal(0, 1, 1000)
    # boot = [np.random.choice(sample, ) for i in range(100)]

    # Plotting the histogram
    plt.hist(boots, bins=30, density=True, alpha=0.6, color='g')

    # Adding the saddlepoint approximation plot
    fhat_values = np.linspace(-3, 3, 100)  # Adjust range and number of points as needed
    density_values = [spa_stdnormal.f_hat(i) for i in fhat_values]
    plt.plot(fhat_values, np.array(density_values), color='red')

    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.title('Bootstrapped Sample Means with Saddlepoint Approximation')
    plt.savefig('spa.png', bbox_inches='tight')

