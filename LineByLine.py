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





## Solving the System using LU Decomposition
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


### Code to generate temperature matrix A, B
def generate_temperature_matrix(n, T0, T1, T):
    # Create the A matrix
    A = np.zeros((n**2, n**2))
    # Create the B matrix
    B = np.zeros(n**2)

    # Iterate over each point in the square plate
    for i in range(n):
        for j in range(n):
            index = i * n + j

            # Set the diagonal element in A matrix
            A[index, index] = 4

            # Set the elements for the surrounding points
            if i > 0:
                A[index, (i-1) * n + j] = -1
            if i < n-1:
                A[index, (i+1) * n + j] = -1
            if j > 0:
                A[index, i * n + (j-1)] = -1
            if j < n-1:
                A[index, i * n + (j+1)] = -1

            # Set the corresponding element in B matrix
            if i == 0:
                B[index] += T0
            if i == n-1:
                B[index] += T1
            if j == 0:
                B[index] += T[index]
            if j == n-1:
                B[index] += T[index]

    # Solve the system of equations Ax = B

    return A, B

## Shooting method for generating the temp distri in surr

def shooting_for_temp_distri_in_heatpipe(h, Ta, T0, Tf, x0, xfinal, delta_x, z01, z02,n):

    # Number of position steps
    num_steps = int((xfinal - x0) / delta_x + 1)

    # Initialize arrays to store results
    T_final = [0, 0]
    z_final = [0, 0]


    z = [0] * num_steps
    T = [0] * num_steps
    z[0] = z01
    T[0] = T0
    x = [0] * num_steps
    # Midpoint method
    for i in range(1, num_steps):
        # Update position
        x[i] = x[i-1] + delta_x
        # Predict using Euler's method
        z_pred = z[i-1] + delta_x * (h*(T[i-1]-Ta))
        T_pred = T[i-1] + delta_x * (z[i-1])

        # Correct using Heun's method
        z[i] = z[i-1] + 0.5 * delta_x * (h*(T[i-1]-Ta) + h*(T_pred-Ta))
        T[i] = T[i-1] + 0.5 * delta_x * (z[i-1]+z_pred)
    T_final[0] = T[-1]
    z_final[0] = z01
    true_guess = z_final[0] + ((z_final[1] - z_final[0])/(T_final[1]- T_final[0]))*(Tf - T_final[0])

    z = [0] * num_steps
    T = [0] * num_steps
    z[0] = z02
    T[0] = T0
    x = [0] * num_steps
    # Midpoint method
    for i in range(1, num_steps):
        # Update position
        x[i] = x[i-1] + delta_x
        # Predict using Euler's method
        z_pred = z[i-1] + delta_x * (h*(T[i-1]-Ta))
        T_pred = T[i-1] + delta_x * (z[i-1])

        # Correct using Heun's method
        z[i] = z[i-1] + 0.5 * delta_x * (h*(T[i-1]-Ta) + h*(T_pred-Ta))
        T[i] = T[i-1] + 0.5 * delta_x * (z[i-1]+z_pred)
    T_final[1] = T[-1]
    z_final[1] = z02


    z = [0] * num_steps
    z[0] = true_guess
    T = [0] * num_steps
    T[0] = T0
    x = [0] * num_steps

    # print('actual initial value of z, hopefully =', z[0])
    for i in range(1, num_steps):
        # Update position
        x[i] = x[i-1] + delta_x
        # Predict using Euler's method
        z_pred = z[i-1] + delta_x * (h*(T[i-1]-Ta))
        T_pred = T[i-1] + delta_x * (z[i-1])

        # Correct using Heun's method
        z[i] = z[i-1] + 0.5 * delta_x * (h*(T[i-1]-Ta) + h*(T_pred-Ta))
        T[i] = T[i-1] + 0.5 * delta_x * (z[i-1]+z_pred)
    # print(T[-1])

    if n<=len(T):
      step = len(T) // (n - 1)
      sampled_T = [T[i] for i in range(0, len(T), step)]
      # print(sampled_T)

    sampled_T = np.array(sampled_T)
    return sampled_T


#########################################################################################################################################################

st.header("MA 203 Project Group 1")
st.title("Modelling the heating in an IC with multiple heat transfer mechanisms using Numerical Methods (Liebmann's Method)")
st.link_button("Go to LU Decomposition Method", "https://ma203-project-group1-49zu9pkqddsqvymnpceqrj.streamlit.app/")
st.link_button("Go to the Project Report", "")
st.sidebar.title("User Inputs")
liquid_input = st.sidebar.slider("Input temperature of the Coolant", 0, 100, 50)
liquid_output = st.sidebar.slider("Vaapour Temperature of the Coolant", 0, 100, 75)


vistype = st.radio("Select the type of Visualisation", ["2D Heatmaps", "3D Surface Plots"])
# Parameters
i = st.sidebar.slider("Select the current Value", 1, 10)
r = st.sidebar.slider("Select the Resistance Value", 1, 100)
c = st.sidebar.slider("Select the C Value for the Power source", 1, 100)
t = st.slider("Select the time duration", 1, 100)
h = 0.01 # For water flowing in metal tubes (coefficient of heat transfer)
Ta = st.sidebar.slider("Select the Surrounding Temperature", 15, 35)

# Initial temp. conditions
T0 = st.sidebar.slider("Select the Heatsink Temperature", 10, 35)
Tf = i**2*r*t/c # At x = xfinal
T1 = Tf

# Pipe length and dist. intervals
x0 = 0
xfinal = 10
delta_x = 0.1

z01 = 10
z02 = 20
n = st.slider("Select the Number of transistors in a row", 3, 13, 3)

###########################################################################################################################################################

TTT = shooting_for_temp_distri_in_heatpipe(h, Ta, T0, Tf, x0, xfinal, delta_x, z01, z02, n)
print(TTT)

lis = []

for i in range(len(TTT)):
  for j in range(len(TTT)):
    lis.append(TTT[i])
T=lis

A, B = generate_temperature_matrix(n, T0, T1, T)


#############################################################################################################################################################



## Line by line solver, liebmann method

import numpy as np

def liebmann_method(A, B, initial_guess, max_iterations=100, tolerance=1e-6):
    n = len(A)
    X = initial_guess.copy()

    for _ in range(max_iterations):
        X_prev = X.copy()

        for i in range(n):
            X[i] = (B[i] - np.dot(A[i, :i], X[:i]) - np.dot(A[i, i+1:], X_prev[i+1:])) / A[i, i]

        if np.linalg.norm(X - X_prev) < tolerance:
            break

    return X

# Example usage

initial_guess = np.zeros_like(B)

solution = liebmann_method(A, B, initial_guess)


x = list(solution)
X = np.array(x).reshape((n,n))

masterlist=[]
lis = [T0]*(n+2)
masterlist.append(np.array(lis))
for i in range(len(TTT)):
  a=TTT[i]
  lis=list(X[i])
  lis.insert(0, a)
  lis.insert(len(TTT)+1, a)
  masterlist.append(np.array(lis))

st.write("Line by Line solver, using Gauss Siedel approach")
lis = [Tf]*(n+2)
masterlist.append(np.array(lis))
# Example 2D NumPy array
data = np.array(masterlist)

if vistype == '2D Heatmaps':
    heatmap_trace = go.Heatmap(z=data, colorscale='Hot')
    layout = go.Layout(title='Heatmap Plot', xaxis=dict(title='X-axis'), yaxis=dict(title='Y-axis'), width=800, height=800)
    fig = go.Figure(data=[heatmap_trace], layout=layout)
    st.plotly_chart(fig)

else:
    sh_0, sh_1 = data.shape
    x, y = np.linspace(0, 1, sh_0), np.linspace(0, 1, sh_1)
    fig = go.Figure(data=[go.Surface(z=data, x=x, y=y, colorscale='Temps')])
    fig.update_layout(scene_aspectmode='manual',
                      scene_aspectratio=dict(x=1, y=1, z=0.3))
    fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                      highlightcolor="limegreen", project_z=True))
    st.plotly_chart(fig)

