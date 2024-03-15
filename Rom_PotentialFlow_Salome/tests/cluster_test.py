import numpy as np
import pickle
from matplotlib import pyplot as plt

#import kmeans
from sklearn.cluster import KMeans
######

def kmeans_test(test_data, n_clusters,mu):
    kmeans_object = KMeans(n_clusters=n_clusters, random_state=0).fit(test_data)
    fig = plt.figure()
    fig.set_figwidth(12.0)
    fig.set_figheight(8.0)
    for j in range(n_clusters):
        fig = plt.scatter(mu[:, 1][kmeans_object.labels_==j], mu[:, 0][kmeans_object.labels_==j], label="Cluster "+str(j))

    fig = plt.title("k-means clustering")
    fig = plt.ylabel('Alpha')
    fig = plt.xlabel('Mach')
    fig = plt.minorticks_on()
    fig = plt.grid(which='both')
    fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
    # fig = plt.show()
    fig = plt.savefig('k-means clustering.pdf',bbox_inches='tight' )
    fig = plt.close('all')
    
    # print("Cluster for '1.448, 0.732.npy': " + str(kmeans_object.predict(np.load("1.448, 0.732.npy").T)))

def load_mu_parameters(name):  
    archivo = open(name, 'rb')
    mu = pickle.load(archivo)
    archivo.close()
    return mu


def E_p(u, c):
    """
    c: direction vector onto which to project
    u: vector or colection of column vectors to project onto the direction of c
    """
    c = c.reshape(-1,1)
    if len(u.shape)==1:
        u = u.reshape(-1,1)
    projection_of_u_onto_c = ((c@c.T) / (c.T@c)) @ u
    projection_error = np.linalg.norm(u - projection_of_u_onto_c, axis=0) / np.linalg.norm(u,axis=0)
    return projection_error

def PEBL(Snapshots, bisection_tolerance=0.15,  POD_tolerance=1e-3):
    #stage 1, generation of bisection tree with accuracy 'bisection_tolerance'
    max_index = np.argmax( np.linalg.norm(Snapshots, axis=0) )
    first_snapshot = Snapshots[:,max_index]
    Tree = Node([first_snapshot, Snapshots])
    bisect_flag = True
    while bisect_flag == True:
        bisect_flag = False
        for leaf in Tree.leaves():
            errors = E_p(leaf.val[1], leaf.val[0])
            max_error = max(errors)
            print(max_error)
            if max_error > bisection_tolerance:
                bisect_flag = True
                #find next anchor point
                max_index = np.argmax(errors)
                c_new = leaf.val[1][:,max_index]
                new_errors = E_p(leaf.val[1], c_new)
                indexes_left = np.where( errors <= new_errors)
                indexes_right = np.where( errors > new_errors)
                #divide the snapshots among the two children
                leaf.left =  Node([leaf.val[0], leaf.val[1][:,indexes_left[0]]])
                leaf.right = Node([c_new, leaf.val[1][:,indexes_right[0]]])
                leaf.val[1] = None
    #stage 2, generation of local POD bases'
    for leaf in Tree.leaves():
        Phi_i = ObtainBasis(leaf.val[1], POD_tolerance)
        leaf.val.append(Phi_i)
    return Tree


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

def pebl_test(test_data,  bisection_tolerance=0.15,  POD_tolerance=1e-3, list=False):
    Tree = PEBL(test_data.T,  bisection_tolerance=bisection_tolerance,  POD_tolerance=POD_tolerance)
    fig = plt.figure()
    fig.set_figwidth(12.0)
    fig.set_figheight(8.0)
    print("leaves: ", len(Tree.leaves()))
    leaf_counter = 0
    for leaf in Tree.leaves():
        leaf_counter+=1
        #save the centroid data (snapshotmatrix)
        # np.save("initial_solution_cluster_"+str(leaf_counter),leaf.val[0])
        print("Number of cases in leaf "+str(leaf_counter)+": ", len(leaf.val[1][1]))
        angle = np.zeros(len(leaf.val[1][1]))
        mach  = np.zeros(len(leaf.val[1][1]))
        for j in range(len(leaf.val[1][1])): 
            for i in range(test_data.T.shape[1]):
                if np.array_equal(test_data.T[:,i], leaf.val[1][:,j]):
                    angle[j] = list[i][0]
                    mach[j]  = list[i][1]
                    break
        fig = plt.scatter(mach, angle, label="Cluster "+str(leaf_counter))
    fig = plt.title("PEBL clustering")
    fig = plt.ylabel('Alpha')
    fig = plt.xlabel('Mach')
    fig = plt.minorticks_on()
    fig = plt.grid(which='both')
    fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
    fig = plt.savefig('PEBL clustering.pdf',bbox_inches='tight' )
    # fig = plt.show()
    fig = plt.close('all')

if __name__=='__main__':

    #load snapshotsmatrix
    snapshotsmatrix = np.load("M_simple.npy")

    #load mu data
    mu      = np.array(np.double(load_mu_parameters('mu_train.dat')))
    mu_list = load_mu_parameters('mu_train.dat')

    #launch tests
    kmeans_test(snapshotsmatrix.T,5,mu)

    pebl_test(snapshotsmatrix.T, bisection_tolerance=0.15,  POD_tolerance=1e-12, list=mu_list)