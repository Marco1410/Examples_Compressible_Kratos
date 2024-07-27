import numpy as np


def E_p(u, c):
    """
    c: direction vector onto which to project.
    u: vector or collection of column vectors to project onto the direction of c.
    """
    c = c.reshape(-1, 1)
    if len(u.shape) == 1:
        u = u.reshape(-1, 1)
    projection_coefficients = (u.T @ c) / (c.T @ c)
    projection_of_u_onto_c = projection_coefficients.T * c
    projection_error = np.linalg.norm(u - projection_of_u_onto_c, axis=0) / np.linalg.norm(u, axis=0)
    return projection_error

def ObtainBasis(Snapshots, truncation_tolerance=0):
        U,_,_= truncated_svd(Snapshots,truncation_tolerance)
        return U

def truncated_svd(Matrix, epsilon=0):
    M,N=np.shape(Matrix)
    dimMATRIX = max(M,N)
    U, s, V = np.linalg.svd(Matrix, full_matrices=False) #U --> M xN, V --> N x N
    V = V.T
    tol = dimMATRIX*np.finfo(float).eps*max(s)/2
    R = np.sum(s > tol)  # Definition of numerical rank
    if epsilon == 0:
        K = R
    else:
        SingVsq = np.multiply(s,s)
        SingVsq.sort()
        normEf2 = np.sqrt(np.cumsum(SingVsq))
        epsilon = epsilon*normEf2[-1] #relative tolerance
        T = (sum(normEf2<epsilon))
        K = len(s)-T
    K = min(R,K)
    return U[:, :K], s[:K], V[:, :K]

class PEBL:
    def __init__(self, bisection_tolerance=0.15,  POD_tolerance=1e-3):
        self.bisection_tolerance = bisection_tolerance
        self.POD_tolerance = POD_tolerance

    def fit(self, Snapshots):
        self.Snapshots = Snapshots
        self.Stage_1()
        self.Stage_2()
        return self

    def Stage_1(self):
        #stage 1, generation of bisection tree with accuracy 'bisection_tolerance'
        max_index = np.argmax( np.linalg.norm(self.Snapshots, axis=0) )
        first_snapshot = self.Snapshots[:,max_index]
        self.Tree = Node([first_snapshot,np.arange(0,self.Snapshots.shape[1], 1, dtype=int)])
        bisect_flag = True
        while bisect_flag == True:
            bisect_flag = False
            for leaf in self.Tree.leaves():
                errors = E_p(self.Snapshots[:,leaf.val[1]], leaf.val[0])
                max_error = max(errors)
                print(f'PEBL clustering max local error {max_error}')
                if max_error > self.bisection_tolerance:
                    bisect_flag = True
                    #find next anchor point
                    max_index = np.argmax(errors)
                    c_new = self.Snapshots[:,leaf.val[1]][:,max_index]
                    new_errors = E_p(self.Snapshots[:,leaf.val[1]], c_new)
                    indexes_left = np.where( errors <= new_errors)
                    indexes_right = np.where( errors > new_errors)
                    #divide the snapshots among the two children
                    leaf.left =  Node([leaf.val[0], leaf.val[1][indexes_left[0]]])
                    leaf.right = Node([c_new, leaf.val[1][indexes_right[0]]])
                    leaf.val[1] = None

    def Stage_2(self):
        #stage 2, generation of local POD bases'
        index = 0
        for leaf in self.Tree.leaves():
            Phi_i = ObtainBasis(self.Snapshots[:,leaf.val[1]], self.POD_tolerance)
            leaf.val.append(Phi_i)
            leaf.val.append(index)
            index+=1

    def predict(self, u):
        current_node = self.Tree
        while not current_node.is_leaf():
            if E_p(u, current_node.left.val[0]) < E_p(u, current_node.right.val[0]):
                current_node = current_node.left
            else:
                current_node = current_node.right
        return current_node.val[3]

class Node:
    def __init__(self, val):
        self.left = None
        self.right = None
        self.val = val

    def leaves(self):
        current_nodes = [self]
        leaves = []

        while len(current_nodes) > 0:
            next_nodes = []
            for node in current_nodes:
                if node.left is None and node.right is None:
                    leaves.append(node)
                    continue
                if node.left is not None:
                    next_nodes.append(node.left)
                if node.right is not None:
                    next_nodes.append(node.right)
            current_nodes = next_nodes
        return leaves

    def is_leaf(self):
        return self.left is None and self.right is None