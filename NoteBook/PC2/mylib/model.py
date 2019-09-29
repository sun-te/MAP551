import numpy as np

class frictionless_pendulum_model(object):

    def __init__(self):
        pass

    def fcn(self, t, y):
        y1, y2 = y
        y1_dot = y2
        y2_dot = -np.sin(y1)
        return (y1_dot, y2_dot)

class friction_pendulum_model(object):

    def __init__(self, alpha):
        self.alpha = alpha

    def fcn(self, t, y):
        alpha = self.alpha
        y1, y2 = y
        y1_dot = y2
        y2_dot = -np.sin(y1) -alpha*y2
        return (y1_dot, y2_dot)

class lokta_model(object):

    def __init__(self, k):
        self.k = k

    def fcn(self, t, u):
        u1, u2 = u
        k = self.k
        u1_dot = u1 * (1 - u2)
        u2_dot = u2 * (-k + u1)
        return (u1_dot, u2_dot)

class lokta_competitive_model(object):

    def __init__(self, alpha, beta):
        self.alpha = alpha
        self.beta = beta

    def fcn(self, t, u):
        u1, u2 = u
        alpha = self.alpha
        beta = self.beta
        u1_dot = u1 * (1 - u1 - u2)
        u2_dot = beta * (u1 - alpha) * u2
        return (u1_dot, u2_dot)

class rosenzweig_model(object):

    def __init__(self, alpha, beta, gamma):
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    def fcn(self, t, u):
        u1, u2 = u
        alpha = self.alpha
        beta = self.beta
        gamma = self.gamma
        u1_dot = u1 * (1 - (u1/gamma)) - ((u1*u2)/(1+u1))
        u2_dot = beta * u2 * ((u1/(1+u1)) - alpha)
        return (u1_dot, u2_dot)
