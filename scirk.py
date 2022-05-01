import numpy as np
import matplotlib.pyplot as plt


class SIRModel:
    def __init__(self, N, beta=0.2, gamma=1 / 10):
        self.N = N
        self.beta = beta
        self.gamma = gamma
        self.init_susceptible = N - 1
        self.init_infected = 1
        self.init_recovered = 0

    def dsdt(self, curr_susceptible, curr_infected, curr_recovered):
        return -(self.beta * curr_susceptible * curr_infected) / self.N

    def didt(self, curr_susceptible, curr_infected, curr_recovered):
        return ((self.beta * curr_susceptible * curr_infected) / self.N) - (self.gamma * curr_infected)

    def drdt(self, curr_susceptible, curr_infected, curr_recovered):
        return self.gamma * curr_infected

    def simulate(self, time_steps):
        # return 3 buffers of length time_steps+1 for the infected/recovered/susceptible values
        susceptible = np.zeros(time_steps + 1)
        infected = np.zeros(time_steps + 1)
        recovered = np.zeros(time_steps + 1)
        susceptible[0] = self.init_susceptible
        infected[0] = self.init_infected
        recovered[0] = self.init_recovered

        for i in range(1, time_steps + 1):
            next_values = self.simulate_step(
                curr_susceptible=susceptible[i - 1],
                curr_infected=infected[i - 1],
                curr_recovered=recovered[i - 1],
            )
            susceptible[i] = next_values[0]
            infected[i] = next_values[1]
            recovered[i] = next_values[2]

        return susceptible, infected, recovered

    @staticmethod
    def RK4(f, y_n):
        k1 = f(y_n)
        k2 = f(y_n + k1 / 2)
        k3 = f(y_n + k2 / 2)
        k4 = f(y_n + k3)
        return y_n + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    def simulate_step(self, curr_susceptible, curr_infected, curr_recovered):
        return self.RK4(
            lambda SIR: np.array([self.dsdt(*SIR), self.didt(*SIR), self.drdt(*SIR)]),
            np.array([curr_susceptible, curr_infected, curr_recovered]))

    def simulate_and_plot(self, time_steps):
        susceptible, infected, recovered = self.simulate(time_steps)
        time_axis = [x for x in range(time_steps + 1)]
        plt.plot(time_axis, infected, 'b', label="Infected")
        plt.plot(time_axis, recovered, 'g', label="Recovered")
        plt.plot(time_axis, susceptible, 'r', label="Susceptible")
        plt.xlabel("Time (days)")
        plt.ylabel("Number of People")
        plt.title(f'SIR Model with γ=1/10 and β=0.2')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    SIRModel(N=1000, beta=0.2).simulate_and_plot(160)
