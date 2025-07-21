#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import copy
import math

# Compute sigmoid for logistic regression
def sigmoid(z):
    """
    Compute the sigmoid of z
    Args:
        z (ndarray): A scalar, numpy array of any size.
    Returns:
        g (ndarray): sigmoid(z), with the same shape as z
    """
    g = 1 / (1 + np.exp(-z))
    return g


def compute_gradient(X, y, w, b, *argv): 
    """
    Computes the gradient for logistic regression
    Args:
        X : (ndarray Shape (m,n)) data, m examples by n features
        y : (ndarray Shape (m,))  target value 
        w : (ndarray Shape (n,))  values of parameters of the model      
        b : (scalar)              value of bias parameter of the model
        *argv : unused, for compatibility with regularized version below
    Returns
        dj_dw : (ndarray Shape (n,)) The gradient of the cost w.r.t. the parameters w. 
        dj_db : (scalar)             The gradient of the cost w.r.t. the parameter b. 
    """
    m, n = X.shape
    dj_dw = np.zeros(w.shape)
    dj_db = 0.
    for i in range(m):
        z_wb = 0
        for j in range(n):
            z_wb += w[j] * X[i][j]
        z_wb += b
        f_wb = sigmoid(z_wb)
        dj_db += f_wb - y[i]
        for j in range(n):
            dj_dw[j] += (f_wb - y[i]) * X[i][j]
    dj_dw = dj_dw / m
    dj_db = dj_db / m
    return dj_db, dj_dw


def gradient_descent(X, y, w_in, b_in, cost_function, gradient_function, alpha, num_iters, lambda_): 
    """
    Performs batch gradient descent to learn theta. Updates theta by taking 
    num_iters gradient steps with learning rate alpha
    Args:
        X :    (ndarray Shape (m, n) data, m examples by n features
        y :    (ndarray Shape (m,))  target value 
        w_in : (ndarray Shape (n,))  Initial values of parameters of the model
        b_in : (scalar)              Initial value of parameter of the model
        cost_function :              function to compute cost
        gradient_function :          function to compute gradient
        alpha : (float)              Learning rate
        num_iters : (int)            number of iterations to run gradient descent
        lambda_ : (scalar, float)    regularization constant
    Returns:
        w : (ndarray Shape (n,)) Updated values of parameters of the model after
            running gradient descent
        b : (scalar)             Updated value of parameter of the model after
            running gradient descent
    """
    m = len(X)
    # An array to store cost J and w's at each iteration
    J_history = []
    w_history = []
    for i in range(num_iters):
        # Calculate the gradient and update the parameters
        dj_db, dj_dw = gradient_function(X, y, w_in, b_in, lambda_)
        # Update Parameters using w, b, alpha and gradient
        w_in = w_in - alpha * dj_dw
        b_in = b_in - alpha * dj_db
        # Save cost J at each iteration
        if i < 100000: # prevent resource exhaustion 
            cost = cost_function(X, y, w_in, b_in, lambda_)
            J_history.append(cost)
        # Print cost every at intervals 10 times or as many iterations if < 10
        if i% math.ceil(num_iters/10) == 0 or i == (num_iters-1):
            w_history.append(w_in)
            print(f"Iteration {i:4}: Cost {float(J_history[-1]):8.2f}   ")
    return w_in, b_in, J_history, w_history


def predict(X, w, b): 
    """
    Predict whether the label is 0 or 1 using learned logistic
    regression parameters w
    Args:
        X : (ndarray Shape (m,n)) data, m examples by n features
        w : (ndarray Shape (n,))  values of parameters of the model      
        b : (scalar)              value of bias parameter of the model

    Returns:
        p : (ndarray (m,)) The predictions for X using a threshold at 0.5
    """
    m, n = X.shape
    p = np.zeros(m)
    # Loop over each example
    for i in range(m):
        z_wb = 0
        # Loop over each feature
        for j in range(n):
            # Add the corresponding term to z_wb
            z_wb += X[i, j] * w[j]
        # Add bias term
        z_wb += b
        # Calculate the prediction for this example
        f_wb = sigmoid(z_wb)
        # Apply the threshold
        p[i] = f_wb >= 0.5
    return p


def compute_cost_reg(X, y, w, b, lambda_ = 1):
    """
    Computes the cost over all examples
    Args:
        X : (ndarray Shape (m,n)) data, m examples by n features
        y : (ndarray Shape (m,))  target value
        w : (ndarray Shape (n,))  values of parameters of the model
        b : (scalar)              value of bias parameter of the model
        lambda_ : (scalar, float) Controls amount of regularization
    Returns:
        total_cost : (scalar)     cost
    """
    m, n = X.shape
    cost_without_reg = compute_cost(X, y, w, b) 
    reg_cost = 0.
    for j in range(n):
        reg_cost += w[j]**2
    reg_cost = lambda_ / (2 * m) * reg_cost
    total_cost = cost_without_reg + reg_cost
    return total_cost


def compute_gradient_reg(X, y, w, b, lambda_ = 1): 
    """
    Computes the gradient for logistic regression with regularization
    Args:
        X : (ndarray Shape (m,n)) data, m examples by n features
        y : (ndarray Shape (m,))  target value
        w : (ndarray Shape (n,))  values of parameters of the model
        b : (scalar)              value of bias parameter of the model
        lambda_ : (scalar,float)  regularization constant
    Returns
        dj_db : (scalar)             The gradient of the cost w.r.t. the parameter b.
        dj_dw : (ndarray Shape (n,)) The gradient of the cost w.r.t. the parameters w.
    """
    m, n = X.shape
    dj_db, dj_dw = compute_gradient(X, y, w, b)
    for j in range(n):
        dj_dw_j_reg = (lambda_ / m) * w[j]
        dj_dw[j] = dj_dw[j] + dj_dw_j_reg
    return dj_db, dj_dw
