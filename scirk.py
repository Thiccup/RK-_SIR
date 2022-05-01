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

    def dsdt(self, curr_infected, curr_susceptible):
        return -(self.beta * curr_susceptible * curr_infected) / self.N

    def didt(self, curr_infected, curr_susceptible):
        return ((self.beta * curr_susceptible * curr_infected) / self.N) - (self.gamma * curr_infected)

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

        for i in range(1, time_steps + 1):
            next_values = self.simulate_step(
                curr_infected=infected[i - 1],
                curr_recovered=recovered[i - 1],
                curr_susceptible=susceptible[i - 1],
            )
            infected[i] = next_values[0]
            recovered[i] = next_values[1]
            susceptible[i] = next_values[2]

        return infected, recovered, susceptible

    def simulate_step(self, curr_infected, curr_recovered, curr_susceptible):
        k11 = self.dsdt(curr_infected, curr_susceptible)
        k12 = self.didt(curr_infected, curr_susceptible)
        k13 = self.drdt(curr_infected)
        k21 = self.dsdt(curr_infected + k12/2, curr_susceptible + k11/2)
        k22 = self.didt(curr_infected + k12/2, curr_susceptible + k11/2)
        k23 = self.drdt(curr_infected + k12/2)
        k31 = self.dsdt(curr_infected + k22/2, curr_susceptible + k21/2)
        k32 = self.didt(curr_infected + k22/2, curr_susceptible + k21/2)
        k33 = self.drdt(curr_infected + k22/2)
        k41 = self.dsdt(curr_infected + k32, curr_susceptible + k31)
        k42 = self.didt(curr_infected + k32, curr_susceptible + k31)
        k43 = self.drdt(curr_infected + k32)
        next_susceptible = curr_susceptible + (1/6) * (k11 + (2 * k21) + (2 * k31) + k41)
        next_infected = curr_infected + (1/6) * (k12 + (2 * k22) + (2 * k32) + k42)
        next_recovered = curr_recovered + (1/6) * (k13 + (2 * k23) + (2 * k33) + k43)
        return next_infected, next_recovered, next_susceptible


if __name__ == '__main__':
    model = SIRModel(N=1000)
    infected, recovered, susceptible = model.simulate(200)
    time_axis = [x for x in range(0, 201)]
    plt.plot(time_axis, infected, 'b', label="Infected")
    plt.plot(time_axis, recovered, 'g', label="Recovered")
    plt.plot(time_axis, susceptible, 'r', label="Susceptible")
    plt.xlabel("Time (days)")
    plt.ylabel("Number of People")
    plt.title(f'SIR Model with γ=1/10 and β=0.2')
    plt.legend()
    plt.show()
