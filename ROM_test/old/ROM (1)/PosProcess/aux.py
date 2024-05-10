
def plot_mu_values_with_errors():
    mu_train, mu_test = load_mu_parameters(with_errors=True)
    markercolor = ["ob","xr","+g","sc","*m","Dy"]
    tolerancias = [1.0, 1e-3]
    fig = plt.figure()
    fig.set_figwidth(12.0)
    fig.set_figheight(8.0)
    if len(mu_train)>0:
        items_filtrados_1 = [item for item in mu_train if item[2] >= tolerancias[0]]
        items_filtrados_11 = [item for item in mu_train if item[2] < tolerancias[0]]
        items_filtrados_2 = [item for item in items_filtrados_11 if item[2] >= tolerancias[1]]
        items_filtrados_3 = [item for item in mu_train if item[2] < tolerancias[1]]

        train_a1 = [item[0] for item in items_filtrados_1]
        train_a2 = [item[0] for item in items_filtrados_2]
        train_a3 = [item[0] for item in items_filtrados_3]

        train_m1 = [item[1] for item in items_filtrados_1]
        train_m2 = [item[1] for item in items_filtrados_2]
        train_m3 = [item[1] for item in items_filtrados_3]

        if len(train_a1) > 0:
            fig = plt.plot(train_m1, train_a1, markercolor[0], label="train error > 1%")
        if len(train_a2) > 0:
            fig = plt.plot(train_m2, train_a2, markercolor[1], label="train 1% > error > 1e-3%")
        if len(train_a3) > 0:
            fig = plt.plot(train_m3, train_a3, markercolor[2], label="train error < 1e-3%")

    if len(mu_test)>0:
        items_filtrados_1 = [item for item in mu_test if item[2] >= tolerancias[0]]
        items_filtrados_11 = [item for item in mu_test if item[2] < tolerancias[0]]
        items_filtrados_2 = [item for item in items_filtrados_11 if item[2] >= tolerancias[1]]
        items_filtrados_3 = [item for item in mu_test if item[2] < tolerancias[1]]

        test_a1 = [item[0] for item in items_filtrados_1]
        test_a2 = [item[0] for item in items_filtrados_2]
        test_a3 = [item[0] for item in items_filtrados_3]

        test_m1 = [item[1] for item in items_filtrados_1]
        test_m2 = [item[1] for item in items_filtrados_2]
        test_m3 = [item[1] for item in items_filtrados_3]

        if len(test_a1) > 0:
            fig = plt.plot(test_m1, test_a1, markercolor[3], label="test error > 1%")
        if len(test_a2) > 0:
           fig = plt.plot(test_m2, test_a2, markercolor[4], label="test 1% > error > 1e-3%")
        if len(test_a3) > 0:
            fig = plt.plot(test_m3, test_a3, markercolor[5], label="test error < 1e-3%")

    fig = plt.title('Mu Values and approximation errors FOM vs ROM')
    fig = plt.ylabel('Alpha')
    fig = plt.xlabel('Mach')
    fig = plt.grid(True)
    fig = plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
    fig = plt.savefig("MuValuesWithErrors.png")
    fig = plt.close('all')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# plot
#

def plot_Cps():
    mu_train, mu_test = load_mu_parameters()
    case_names = ["FOM","ROM","HROM"]
    markercolor = ["ob","xr","+g","sc","*m","Dy"]
    for j in range(len(mu_train)):
        if os.path.exists("DataBase/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat"):
            cp_min = 0
            cp_max = 0
            fig = plt.figure()
            fig.set_figwidth(12.0)
            fig.set_figheight(8.0)
            for n, name in enumerate(case_names):
                casename = "Fit"
                if name == "FOM":
                    x  = np.loadtxt("DataBase/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(0,))
                    cp = np.loadtxt("DataBase/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(3,))
                    fig = plt.plot(x, cp, markercolor[n], markersize = 2.0, label = name)
                    if np.min(cp) < cp_min:
                        cp_min = np.min(cp)
                    if np.max(cp) > cp_max:
                        cp_max = np.max(cp)
                else:
                    if os.path.exists("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat"):
                        x  = np.loadtxt("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(0,))
                        cp = np.loadtxt("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(3,))
                        fig = plt.plot(x, cp, markercolor[n], markersize = 2.0, label = name)
                        if np.min(cp) < cp_min:
                            cp_min = np.min(cp)
                        if np.max(cp) > cp_max:
                            cp_max = np.max(cp)
            fig = plt.title('Cp vs x - ' + casename + " " + str(mu_train[j][0]) + ", " + str(mu_train[j][1]))
            fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
            fig = plt.ylabel('Cp')
            fig = plt.xlabel('x')
            fig = plt.grid()
            fig = plt.legend()
            fig = plt.tight_layout()
            fig = plt.savefig("Captures/" + casename + str(j) + ".png")
            fig = plt.close('all')
        else:
            print("The file DataBase/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat doesn't exist.")

    for j in range(len(mu_test)):
        if os.path.exists("DataBase/Data/" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat"):
            cp_min = 0
            cp_max = 0
            fig = plt.figure()
            fig.set_figwidth(12.0)
            fig.set_figheight(8.0)
            for n, name in enumerate(case_names):
                casename = "Test"
                if name == "FOM":
                    x  = np.loadtxt("DataBase/Data/" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(0,))
                    cp = np.loadtxt("DataBase/Data/" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(3,))
                    fig = plt.plot(x, cp, markercolor[n], markersize = 2.0, label = name)
                    if np.min(cp) < cp_min:
                        cp_min = np.min(cp)
                    if np.max(cp) > cp_max:
                        cp_max = np.max(cp)
                else:
                    if os.path.exists("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat"):
                        x  = np.loadtxt("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(0,))
                        cp = np.loadtxt("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(3,))
                        fig = plt.plot(x, cp, markercolor[n], markersize = 2.0, label = name)
                        if np.min(cp) < cp_min:
                            cp_min = np.min(cp)
                        if np.max(cp) > cp_max:
                            cp_max = np.max(cp)
            fig = plt.title('Cp vs x - ' + casename + " " + str(mu_test[j][0]) + ", " + str(mu_test[j][1]))
            fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
            fig = plt.ylabel('Cp')
            fig = plt.xlabel('x')
            fig = plt.grid()
            fig = plt.legend()
            fig = plt.tight_layout()
            fig = plt.savefig("Captures/" + casename + str(j) + ".png")
            fig = plt.close('all')
        else:
            print("The file DataBase/Data/" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat doesn't exist.")

def PlotClusterListWithRomModes(self,n_clusters,mu_list,cluster_list,rom_modes_list,clustering_method):
    fig = plt.figure()
    fig.set_figwidth(15.0)
    fig.set_figheight(10.0)
    for i in range(n_clusters):
        mu_aux = [mu for mu, cluster in zip(mu_list, cluster_list[:][0]) if cluster == i]
        mu_train_a = np.zeros(len(mu_aux))
        mu_train_m = np.zeros(len(mu_aux))
        for j in range(len(mu_aux)):
            mu_train_a[j] = mu_aux[j][0]
            mu_train_m[j] = mu_aux[j][1]
        fig = plt.scatter(mu_train_m, mu_train_a, label="Cluster " + str(i) + ": Rom modes: " + str(np.int0(rom_modes_list[i])))
        fig = plt.scatter(self.centroids_to_plot[i,0], self.centroids_to_plot[i,1], c='k',marker="s", s= 150)
        fig = plt.text(self.centroids_to_plot[i,0], self.centroids_to_plot[i,1]-0.05, "Centroid C"+str(i), ha='center')
    fig = plt.title(clustering_method+" clustering")
    fig = plt.ylabel('Alpha')
    fig = plt.xlabel('Mach')
    fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
    fig = plt.savefig(clustering_method+' clustering.pdf',bbox_inches='tight' )
    fig = plt.close('all')