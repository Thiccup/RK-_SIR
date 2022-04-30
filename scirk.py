import numpy as np
import matplotlib.pyplot as plt

beta = 0.2
gamma = 1 / 10
N = 1000
I0, R0, S0 = 1, 0, N-1
start = 1
end = 100
dh = 1
timex = [x for x in range(1,100)]

def dsdt(st, it):
    return -(-beta * st * it) / N


def didt(st, it):
    return ((beta * st * it) / N) - (gamma * it)


def drdt(st, it):
    return gamma * it



def rk4():
    finallist = []
    k1 = df()