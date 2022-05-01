import numpy as np
import matplotlib.pyplot as plt


class SIRModel:
    def __init__(self, N, beta=0.2, gamma=1 / 10):
        self.N = N
        self.beta = beta
        self.gamma = gamma
        self.init_infected = 1
        self.init_recovered = 0
        self.init_susceptible = N - 1

    def dsdt(self, curr_susceptible, curr_infected):
        return -(-self.beta * curr_susceptible * curr_infected) / self.N

    def didt(self, curr_susceptible, curr_infected):
        return ((self.beta * curr_susceptible * curr_infected) / self.N) - (self.gamma * it)

    def drdt(self, curr_infected):
        return self.gamma * curr_infected

    def simulate(self, time_steps):
        # return 3 buffers of length time_steps+1 for the infected/recovered/susceptible values
        infected = np.zeros(time_steps + 1)
        recovered = np.zeros(time_steps + 1)
        susceptible = np.zeros(time_steps + 1)
        infected[0] = self.init_infected
        recovered[0] = self.init_recovered
        susceptible[0] = self.init_susceptible
        
        for i in range(1, time_steps+1):
            next_values = self.simulate_step(
                curr_infected=infected[i-1],
                curr_recovered=recovered[i-1],
                curr_susceptible=susceptible[i-1],
            )
            infected[i] = next_values[0]
            recovered[i] = next_values[1]
            susceptible[i] = next_values[2]
            
        return infected, recovered, susceptible

    def simulate_step(self, curr_infected, curr_recovered, curr_susceptible):
        # Compute next values of infected/recovered/susceptible
        next_infected, next_recovered, next_susceptible = curr_infected, curr_recovered, curr_susceptible
        return next_infected, next_recovered, next_susceptible


if __name__ == '__main__':
    model = SIRModel(N=1000)
    infected, recovered, susceptible = model.simulate(100)
