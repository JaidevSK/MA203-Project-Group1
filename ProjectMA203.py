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

st.header("MA 203 Project Group 1")
st.title("Modelling the heating in an IC with multiple heat transfer mechanisms using Numerical Methods")
st.sidebar.title("User Inputs")
liquid_input = st.sidebar.slider("Input temperature of the Coolant", 0, 100, 50)
liquid_output = st.sidebar.slider("Vaapour Temperature of the Coolant", 0, 100, 75)
n = st.sidebar.radio("Number of Transistors in a row", [2,3,4,5,6])
t = st.slider("Select the time duration", 0, 100)

st.write(f"The value of time is {t},  the value of number of transistors is {n},  liquid_temp is {liquid_input} and {liquid_output}")