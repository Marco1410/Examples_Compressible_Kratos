import random
import numpy as np
import matplotlib.pyplot as plt


class RBF_Interpolator(object):
    def __init__(self, eps, dim, function='gauss'):
        self.eps = eps
        self.dim = dim
        self.function = function
    
    def fit(self, xk, yk):
        self.xk_ = xk
        if self.dim == 1:
            if self.function == 'linear':
                transformation = euclidean_distance(xk, xk)
            elif self.function == 'gauss':
                transformation = gauss_rbf(euclidean_distance(xk, xk), self.eps)
            elif self.function == 'inverse_quadratic':
                transformation = inverse_quadratic_rbf(euclidean_distance(xk, xk), self.eps)
            elif self.function == 'multiquadratic':
                transformation = multiquadratic_rbf(euclidean_distance(xk, xk), self.eps)
        elif self.dim == 2:
            if self.function == 'linear':
                transformation = euclidean_2d_distance(xk, xk)
            elif self.function == 'gauss':
                transformation = gauss_rbf(euclidean_2d_distance(xk, xk), self.eps)
            elif self.function == 'inverse_quadratic':
                transformation = inverse_quadratic_rbf(euclidean_2d_distance(xk, xk), self.eps)
            elif self.function == 'multiquadratic':
                transformation = multiquadratic_rbf(euclidean_2d_distance(xk, xk), self.eps)
        self.w_ = np.linalg.solve(transformation, yk)

    def __call__(self, xn):
        if self.dim == 1:
            if self.function == 'linear':
                transformation = euclidean_distance(xn, self.xk_)
            elif self.function == 'gauss':
                transformation = gauss_rbf(euclidean_distance(xn, self.xk_), self.eps)
            elif self.function == 'inverse_quadratic':
                transformation = inverse_quadratic_rbf(euclidean_distance(xn, self.xk_), self.eps)
            elif self.function == 'multiquadratic':
                transformation = multiquadratic_rbf(euclidean_distance(xn, self.xk_), self.eps)
        elif self.dim == 2:
            if self.function == 'linear':
                transformation = euclidean_2d_distance(xn, self.xk_)
            elif self.function == 'gauss':
                transformation = gauss_rbf(euclidean_2d_distance(xn, self.xk_), self.eps)
            elif self.function == 'inverse_quadratic':
                transformation = inverse_quadratic_rbf(euclidean_2d_distance(xn, self.xk_), self.eps)
            elif self.function == 'multiquadratic':
                transformation = multiquadratic_rbf(euclidean_2d_distance(xn, self.xk_), self.eps)
        return transformation.dot(self.w_)
    
#########################################################################################
#########################################################################################
############ Kernels
#########################################################################################
#########################################################################################
def gauss_rbf(radius, eps): return np.exp(-(eps*radius)**2)

def inverse_quadratic_rbf(radius, eps): return 1/(1+(eps*radius)**2)

def multiquadratic_rbf(radius, eps): return np.sqrt(1+(eps*radius)**2)

#########################################################################################
#########################################################################################
############ 1D Functions
#########################################################################################
#########################################################################################
def true_fn(x): return x**2-x-np.cos(np.pi*x)

def euclidean_distance(x, xk):
    return np.sqrt(((x.reshape(-1,1))-xk.reshape(1,-1))**2)

#########################################################################################
#########################################################################################
############ 2D Functions
#########################################################################################
#########################################################################################
def func(x1, x2): return x1**2*np.cos(2*np.pi*x2)

def generate_data():
    x1 = np.linspace(0, 1, 10)
    x2 = np.linspace(0, 1, 10)
    xx, xy = np.meshgrid(x1, x2)
    y = func(xx, xy)
    return np.concatenate([xx.reshape(-1,1), xy.reshape(-1,1)], axis=1), y.reshape(-1,1)

def plot_real_function(func, resolution=100):
    x1 = np.linspace(0, 1, resolution)
    x2 = np.linspace(0, 1, resolution)
    xx, xy = np.meshgrid(x1, x2)
    y = func(xx, xy)
    cs = plt.contour(xx, xy, y)
    plt.clabel(cs)

def plot_interp(func, resolution=100):
    x1 = np.linspace(0, 1, resolution)
    x2 = np.linspace(0, 1, resolution)
    xx, xy = np.meshgrid(x1, x2)
    data = np.concatenate([xx.reshape(-1,1), xy.reshape(-1,1)],axis=1)
    y = func(data)
    cs = plt.contour(xx, xy, y.reshape(resolution, resolution))
    plt.clabel(cs)

def euclidean_2d_distance(x, xk):
    x = x[np.newaxis,...].swapaxes(0,2)
    xk = xk.T[np.newaxis,...].swapaxes(0,1)
    return np.sqrt(((x-xk)**2).sum(axis=0))


if '__main__':

    ########################################################
    ########################################################
    ###########  1D RBF Interpolator
    ########################################################
    ########################################################
    # xk = np.linspace(0,1,5)
    # yk = true_fn(xk)

    # x = np.linspace(0,1,100)

    # linear_interp = RBF_Interpolator(eps=0, dim=1, function='linear')
    # linear_interp.fit(xk, yk)

    # gauss_interp = RBF_Interpolator(eps=2, dim=1, function='gauss')
    # gauss_interp.fit(xk, yk)

    # inverse_quadratic_interp = RBF_Interpolator(eps=3, dim=1, function='inverse_quadratic')
    # inverse_quadratic_interp.fit(xk, yk)

    # multiquadratic_interp = RBF_Interpolator(eps=1, dim=1, function='multiquadratic')
    # multiquadratic_interp.fit(xk, yk)

    # plt.figure(figsize=(9,6))
    # plt.plot(xk, yk        , 'x'  , label='Points')
    # plt.plot(x , linear_interp(x) , '--', label='Linear interpolation')
    # plt.plot(x , gauss_interp(x) , '--', label='Gaussian interpolation')
    # plt.plot(x , inverse_quadratic_interp(x) , '--', label='Inverse quadratic interpolation')
    # plt.plot(x , multiquadratic_interp(x) , '--', label='Multiquadratic interpolation')
    # plt.plot(x , true_fn(x), 'b'  , label='Original function')
    # plt.legend()
    # plt.title('1D RBF Interpolator')
    # plt.show()


    ########################################################
    ########################################################
    x_gtz      = np.loadtxt('2.0, 0.75.dat', usecols=(0,))[79:][::-1]
    cp_val_gtz = -np.loadtxt('2.0, 0.75.dat', usecols=(1,))[79:][::-1]
    x_lez      = np.loadtxt('2.0, 0.75.dat', usecols=(0,))[:78]
    cp_val_lez = -np.loadtxt('2.0, 0.75.dat', usecols=(1,))[:78]

    xk = x_gtz
    yk = cp_val_gtz

    ordered_indexes = sorted(range(len(xk)), key=lambda i: xk[i])
    xk = np.array([xk[i] for i in ordered_indexes])
    yk = np.array([yk[i] for i in ordered_indexes])

    # indexes = random.sample(range(len(xk)), 25)
    # indexes = sorted(indexes)
    # print(indexes)
    
    indexes = [0,5,10,15,20,22,25,27,30,32,35,37,40,42,45,47,50,52,55,57,60,62,65,67,70,72,75,77,80,82,84,89,92,95,100,112]
    xk = xk[indexes]
    yk = yk[indexes]

    x = np.linspace(0,1,100)

    linear_interp = RBF_Interpolator(eps=1, dim=1, function='linear')
    linear_interp.fit(xk, yk)

    gauss_interp = RBF_Interpolator(eps=15, dim=1, function='gauss')
    gauss_interp.fit(xk, yk)

    # inverse_quadratic_interp = RBF_Interpolator(eps=7, dim=1, function='inverse_quadratic')
    # inverse_quadratic_interp.fit(xk, yk)

    # multiquadratic_interp = RBF_Interpolator(eps=7, dim=1, function='multiquadratic')
    # multiquadratic_interp.fit(xk, yk)

    plt.figure(figsize=(9,6))
    plt.plot(x , linear_interp(x)            , '--', label='Linear interpolation e=1')
    plt.plot(x , gauss_interp(x)             , '--', label='Gaussian interpolation e=15')
    # plt.plot(x , inverse_quadratic_interp(x) , '--', label='Inverse quadratic interpolation')
    # plt.plot(x , multiquadratic_interp(x)    , '--', label='Multiquadratic interpolation')
    plt.plot(xk, yk                          , 'x' , label='Original points')
    plt.legend()
    plt.title('1D RBF Interpolator')
    plt.show()


    ########################################################
    ########################################################
    ###########  2D RBF Interpolator
    ########################################################
    ########################################################
    # data, ytest = generate_data()

    # plt.figure(figsize=(12,8))
    # plot_real_function(func)
    # plt.scatter(data[:,0], data[:,1])

    # interp = RBF_Interpolator(eps=3, dim=2)
    # interp.fit(data, ytest)

    # plot_interp(interp)
    # plt.title('2D RBF Gaussian Interpolator')
    # plt.show()


    ########################################################
    ########################################################
    # kratos_skin_data_filename = f"Onera_Cp_data.dat"
    # x  = np.loadtxt(kratos_skin_data_filename, usecols=(0,))
    # y  = np.loadtxt(kratos_skin_data_filename, usecols=(1,))
    # z  = np.loadtxt(kratos_skin_data_filename, usecols=(2,))
    # cp = np.loadtxt(kratos_skin_data_filename, usecols=(3,))

    # indexes  = z > 0
    # xx = x[indexes]
    # xy = y[indexes]
    # xz = cp[indexes]

    # ordered_indexes = sorted(range(len(xx)), key=lambda i: xx[i])
    # xx = np.array([xx[i] for i in ordered_indexes])
    # xy = np.array([xy[i] for i in ordered_indexes])
    # xz = np.array([xz[i] for i in ordered_indexes])

    # indexes = random.sample(range(len(xx)), 50)
    # indexes = sorted(indexes)
    # # print(indexes)

    # xf  = xx[indexes]
    # yf  = xy[indexes]
    # cpf = xz[indexes]

    # xf, yf = np.meshgrid(xf, yf)

    # data = np.concatenate([xf.reshape(-1,1), yf.reshape(-1,1)], axis=1)
    # ytest = cpf.reshape(-1,1)

    # plt.figure(figsize=(12,8))
    # cs = plt.contour(xf, yf, xz)
    # plt.clabel(cs)
    # plt.show()