# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 15:20:23 2023

@author: jaidev
"""


import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as ex
import plotly.graph_objects as go


def lu_decomposition(matrix):
    n = len(matrix)
    lower = np.zeros((n, n))
    upper = np.zeros((n, n))

    for i in range(n):
        # Upper Triangular
        for k in range(i, n):
            sum = 0
            for j in range(i):
                sum += (lower[i][j] * upper[j][k])
            upper[i][k] = matrix[i][k] - sum

        # Lower Triangular
        for k in range(i, n):
            if i == k:
                lower[i][i] = 1
            else:
                sum = 0
                for j in range(i):
                    sum += (lower[k][j] * upper[j][i])
                lower[k][i] = (matrix[k][i] - sum) / upper[i][i]

    return lower, upper

def forward_substitution(L, b):
    n = len(L)
    y = np.zeros(n)

    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i][j] * y[j]
        y[i] /= L[i][i]

    return y

def backward_substitution(U, y):
    n = len(U)
    x = np.zeros(n)

    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= U[i][j] * x[j]
        x[i] /= U[i][i]

    return x

A = np.array([[3, 3, -8, -2], [-8, -4, 8, 4], [-6, 2, 9, -1], [-2, -7, 1, -10]])
b = np.array([-60, 36, 42, -33])

L, U = lu_decomposition(A)
y = forward_substitution(L, b)
x = backward_substitution(U, y)

print("Solution:", x)





st.header("MA 203 Project Group 1")
st.title("Modelling the heating in an IC with multiple heat transfer mechanisms using Numerical Methods")
st.sidebar.title("User Inputs")
liquid_input = st.sidebar.slider("Input temperature of the Coolant", 0, 100, 50)
liquid_output = st.sidebar.slider("Vaapour Temperature of the Coolant", 0, 100, 75)
n = st.sidebar.radio("Number of Transistors in a row", [2,3,4,5,6])
t = st.slider("Select the time duration", 0, 100)

st.write(f"The value of time is {t},  the value of number of transistors is {n},  liquid_temp is {liquid_input} and {liquid_output}")
