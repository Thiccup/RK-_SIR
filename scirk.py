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
        # return 3 buffers of length time_steps for the infected/recovered/susceptible values
        return None, None, None


model = SIRModel(N=1000)
infected, recovered, susceptible = model.simulate(100)
