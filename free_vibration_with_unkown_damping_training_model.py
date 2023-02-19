import deepxde as dde
import numpy as np

from deepxde.backend import tf


def pde(x, y):
    dy_dt = dde.grad.jacobian(y, x, i=0, j=1)
    d2y_dt2 = dde.grad.hessian(y, x, i=0, j=1)
    
    return (
        d2y_dt2 +x[:, 0:1]*dy_dt + y
        )


    


def func(x):
    return np.sin(np.pi * x[:, 0:1]) * np.exp(-x[:, 1:])


geom = dde.geometry.Interval(0, 1)
timedomain = dde.geometry.TimeDomain(0, 20)
geomtime = dde.geometry.GeometryXTime(geom, timedomain)


ic_1 = dde.icbc.IC(geomtime, lambda x: -1, lambda _, on_initial: on_initial)


#//////////////////////////////////////////////////////////////
ic_2 = dde.icbc.OperatorBC(
    geomtime,
    lambda x, y, _: dde.grad.jacobian(y, x, i=0, j=1)-2,
    lambda _, on_initial: on_initial,
)


#//////////////////////////////////////////////////////////////

data = dde.data.TimePDE(
    geomtime,
    pde,
    [ic_1, ic_2],
    num_domain=400,
    num_boundary=20,
    num_initial=10,
    solution=func,
)

layer_size = [2] + [80] * 5 + [1]
activation = "tanh"
initializer = "Glorot uniform"
net = dde.nn.FNN(layer_size, activation, initializer)

model = dde.Model(data, net)

model.compile("adam", lr=0.001, metrics=["l2 relative error"])
losshistory, train_state = model.train(iterations=20000)

model.compile('L-BFGS', metrics=['l2 relative error'])
losshistory, train_state = model.train()

dde.saveplot(losshistory, train_state, issave=True, isplot=True)


