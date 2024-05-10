import os
import numpy as np
from ezyrb import Database
from ezyrb import ReducedOrderModel as ROM
from ezyrb import RBF, POD
from matplotlib import pyplot as plt

def load_mu_parameters():
    if os.path.exists("../mu_parameters/mu_train.npy") and os.path.exists("../mu_parameters/mu_test.npy"):
        archivo = '../mu_parameters/mu_train.npy'
        mu_train = np.load(archivo)
        archivo = '../mu_parameters/mu_test.npy'
        mu_test = np.load(archivo)
    elif os.path.exists("../mu_parameters/mu_train.npy"):
        archivo = '../mu_parameters/mu_train.npy'
        mu_train = np.load(archivo)
        mu_test = []
    elif os.path.exists("../mu_parameters/mu_test.npy"):
        archivo = '../mu_parameters/mu_test.npy'
        mu_test = np.load(archivo)
        mu_train = []
    return np.array(mu_train), np.array(mu_test)

def error_estimation(mu_list, mesh_database):
    approximation_error = 0.0
    FOM_model = []; RBF_model = []
    for mu in mu_list:
        FOM_model.append(np.load(f'../DataBase/{mesh_database}/Snapshots/{mu[0]}, {mu[1]}.npy'))
        RBF_model.append(np.load(f"{mesh_database}/Interpolated_data/{mu[0]}, {mu[1]}.npy").T)
    FOM_model = np.block(FOM_model)
    RBF_model = np.block(RBF_model)
    print(FOM_model.shape, RBF_model.shape)
    approximation_error = np.linalg.norm(FOM_model - RBF_model)/np.linalg.norm(FOM_model)
    print(f'Clustering approximation error: {approximation_error*100} %')

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
                print(f'Clustering local max error {max_error}')
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

if __name__ == "__main__":

    with os.scandir(f'../DataBase') as folders:
        folders = [dir.name for dir in folders if dir.is_dir()]

    for mesh_database in folders:

        if not os.path.exists(f"{mesh_database}"):
            os.mkdir(f"{mesh_database}")

        bisection_tolerance      = 0.05
        POD_tolerance            = 1e-12

        mu_train, mu_test = load_mu_parameters()

        plot_custering_info      = True
        plot_clusters            = True
        plot_clusters_prediction = True
        
        #### DATA
        #################################################################################################
        parameters = np.array(mu_train)
        snapshots = np.load(f'../DataBase/{mesh_database}/SnapshotsMatrix_conv.npy')

        if plot_custering_info:print(f'Snapshots Matrix shape: {snapshots.shape}')
        if plot_custering_info:print(f'Parameters shape: {parameters.shape}')

        #### RBF TRAINING
        #################################################################################################
        db = Database(parameters, snapshots.T)
        pod = POD()
        rbf = RBF()
        rom = ROM(db, pod, rbf).fit()

        #### PEBL CLUSTERING 
        #################################################################################################
        if not os.path.exists(f"{mesh_database}/Centroids"):
            os.mkdir(f"{mesh_database}/Centroids")
        pebl_object = PEBL(bisection_tolerance=bisection_tolerance,  POD_tolerance=POD_tolerance).fit(snapshots)
        leaves = pebl_object.Tree.leaves()
        mu_train_aux = np.array(mu_train)
        mu_train_with_leaves = []
        leaf_id = 0
        for leaf in leaves:
            for i in range(len(leaf.val[1])):
                mu_train_with_leaves.append([mu_train_aux[leaf.val[1]][i,0], mu_train_aux[leaf.val[1]][i,1], leaf_id])
            np.save(f'{mesh_database}/Centroids/{leaf_id}', leaf.val[0])
            leaf_id += 1
        np.save(f'{mesh_database}/mu_train_with_indexes', mu_train_with_leaves)

        if plot_custering_info: print(f'Number of leaves (clusters): {len(leaves)}')
        
        if plot_clusters:
            #### PLOT LEAVES / CLUSTERS
            #################################################################################################
            fig = plt.figure()
            fig.set_figwidth(24.0)
            fig.set_figheight(16.0)
            leaf_id = 0
            for leaf in leaves:
                if plot_custering_info: print("Number of cases in leaf "+str(leaf_id)+": ", len(leaf.val[1]))
                fig = plt.scatter(mu_train_aux[leaf.val[1]][:,1], mu_train_aux[leaf.val[1]][:,0], label="Cluster "+str(leaf_id))
                fig = plt.scatter(np.mean(mu_train_aux[leaf.val[1]][:,1]), np.mean(mu_train_aux[leaf.val[1]][:,0]), c='k',marker=f'${leaf_id}$', s= 150)
                leaf_id += 1
            fig = plt.title("PEBL clustering")
            fig = plt.ylabel('Alpha')
            fig = plt.xlabel('Mach')
            fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
            fig = plt.savefig(f'{mesh_database}/PEBL clustering.png',bbox_inches='tight' )
            fig = plt.close('all')

        if len(mu_test) > 0:
            #### SOLUTION AND CLUSTER PREDICTION TEST
            #################################################################################################
            if not os.path.exists(f"{mesh_database}/Interpolated_data"):
                os.mkdir(f"{mesh_database}/Interpolated_data")
            interpolated_solutions_list = [rom.predict([element]).snapshots_matrix for element in mu_test]
            predicted_indexes_list = [pebl_object.predict(solution.T) for solution in interpolated_solutions_list]
            for i, solution in enumerate(interpolated_solutions_list):
                np.save(f"{mesh_database}/Interpolated_data/{mu_test[i][0]}, {mu_test[i][1]}.npy", solution)
            np.save(f'{mesh_database}/predicted_indexes', predicted_indexes_list)

            if plot_custering_info: error_estimation(mu_test, mesh_database)

            if plot_clusters_prediction:
                #### PLOT LEAVES / CLUSTERS PREDICTION
                #################################################################################################
                if not os.path.exists(f"{mesh_database}/Clustering_mu_test_plots"):
                    os.mkdir(f"{mesh_database}/Clustering_mu_test_plots")
                for i, predicted_index in enumerate(predicted_indexes_list):
                    if plot_custering_info: print(f'Test {i}: {mu_test[i][1]} {mu_test[i][0]}, predicted index: {predicted_index}')
                    #### PLOT CLUSTERS PREDICTION
                    fig = plt.figure()
                    fig.set_figwidth(24.0)
                    fig.set_figheight(16.0)
                    leaf_id = 0
                    for leaf in leaves:
                        if leaf_id == predicted_index:
                            fig = plt.scatter(mu_train_aux[leaf.val[1]][:,1], mu_train_aux[leaf.val[1]][:,0], label="Cluster "+str(leaf_id))
                            fig = plt.scatter(np.mean(mu_train_aux[leaf.val[1]][:,1]), np.mean(mu_train_aux[leaf.val[1]][:,0]), c='k',marker=f'${leaf_id}$', s= 150)
                        else:
                            fig = plt.scatter(mu_train_aux[leaf.val[1]][:,1], mu_train_aux[leaf.val[1]][:,0], label="Cluster "+str(leaf_id), alpha=0.2)
                            fig = plt.scatter(np.mean(mu_train_aux[leaf.val[1]][:,1]), np.mean(mu_train_aux[leaf.val[1]][:,0]), c='k',marker=f'${leaf_id}$', alpha=0.2, s= 150)
                        leaf_id += 1
                    fig = plt.scatter(mu_test[i][1], mu_test[i][0], c='k',marker="s", s=150)
                    fig = plt.text(mu_test[i][1], mu_test[i][0]-0.05, 
                                   f'Mu test {i}', ha='center',  
                                   bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),))
                    fig = plt.title("PEBL online prediction")
                    fig = plt.ylabel('Alpha')
                    fig = plt.xlabel('Mach')
                    fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
                    fig = plt.savefig(f'{mesh_database}/Clustering_mu_test_plots/PEBL prediction {i} {mu_test[i][1]} {mu_test[i][0]}.png',bbox_inches='tight' )
                    fig = plt.close('all')   
