
import os
import numpy as np
import pandas as pd
from scipy.spatial import distance
from itertools import combinations
from scipy.stats import hypergeom,gaussian_kde,norm
import gc
from statsmodels.stats.multitest import multipletests
import scanpy as sc


class RedeTIL_Features:
    def __init__(self, data, target, perturbation, combo_target = None,
                 T_cells=None, Cancer_cells=None,
                 outdir='./results',
                 genesets = ['CD4','CD8A',"TNFRSF9","CTLA4","TIGIT","CXCL13","TCF7","HAVCR2","LAYN","CXCR5","ENTPD1","LAG3","GZMB","BATF","CCL5",'CD69',"CD3E","CD274","PDCD1"]):
        """
        Initialize the RedeTIL_Features class.

        :param data: AnnData object containing single-cell data, data.obs must contain the 'label' column, where the cell annotation are listed.
        :param target: Target gene.
        :param perturbation: Perturbation type, 'block' or 'depletion'.
        :param combo_target: Combo target gene (optional).
        :param T_cells: Label for T cells (optional), if not provided, Dynamic_features are not available.
        :param Cancer_cells: Label for cancer cells (optional), if not provided, Dynamic_features are not available.
        :param outdir: Output directory.
        :param genesets: List of gene+ T cell subsets of interest.
        """
        self.TPM = pd.DataFrame(data.X.toarray(),index = data.obs.index, columns = data.var.index).T
        self.label = pd.DataFrame(data.obs['labels'].reset_index().values, columns = ['cells','labels'])
        self.target = target
        if perturbation not in ["block", "depletion"]:
            raise ValueError(f"perturbation must be one of 'block' or 'depletion', but got '{perturbation}'")
        self.combo_target = combo_target
        self.perturbation = perturbation
        self.T_cells = T_cells
        self.Cancer_cells = Cancer_cells
        self.outdir = outdir
        self.genesets = genesets
        os.makedirs(self.outdir, exist_ok=True)
        
    def getSignificance(self, coordinates, labels,  k=3):
        """
        Calculate cell connections.

        :param coordinates: Cell coordinates.
        :param labels: Cell labels.
        :param k: Dimensions of the projection space, only support k=3 for now.
        :return: Results of connections, including cell-cell distance, the number cell-cell interaction pairs.
        """
        
        # Check Input 
        assert isinstance(coordinates, np.ndarray), "coordinates is not an array"
        assert 2 <= coordinates.shape[1] <=3, "Unsupported dimensional number"

        # Data preprocessing
        standards = np.unique(labels['labels'])
        cellCounts=pd.Series(labels['labels']).value_counts().to_dict()

        # Calculate distances
        dist_mt = distance.cdist(coordinates, coordinates, 'euclidean')
        euclidean_distance = pd.DataFrame(dist_mt,index= labels['labels'],columns = labels['labels'])
        euclidean_distance = euclidean_distance.groupby(euclidean_distance.index).mean()
        euclidean_distance = euclidean_distance.T.groupby(euclidean_distance.T.index).mean().T #Mean distance

        # Mean rank distance
        rank_distance = pd.DataFrame(dist_mt).rank()-1
        rank_distance.index,rank_distance.columns = ['To-' + x for x in labels['labels']],labels['labels']
        rank_distance = rank_distance.groupby(rank_distance.index).mean()
        rank_distance = rank_distance.T.groupby(rank_distance.T.index).mean().T #Mean rank distance

        # Connections
        diagonal_mask = np.eye(len(dist_mt), dtype=bool)
        dist_mt[diagonal_mask] = np.nan
        topK = np.zeros(dist_mt.shape[0])
        for i in range(dist_mt.shape[0]):
            topK[i]=sorted(dist_mt[i][dist_mt[i]>0])[k-1]
        topK = np.median(topK)
        
        connects_mt = pd.DataFrame(dist_mt<=topK)
        connects_mt.index,connects_mt.columns = labels['labels'],labels['labels']
        counts = connects_mt.groupby(connects_mt.index).sum()
        counts = counts.T.groupby(counts.T.index).sum()
        diagonal_mask = np.eye(len(counts), dtype=bool)
        counts[diagonal_mask] /= 2
        diag = np.diag(counts)

        # Calculate connection details
        #detailed_connections = {}
        detailed_connections = pd.DataFrame()
        clusterPairs2run = np.vstack((np.asarray(list(combinations(standards, 2))),
                                    np.tile(standards, (2, 1)).T))
        for i_row in range(clusterPairs2run.shape[0]):
            cluster1 = clusterPairs2run[i_row, 0]
            cluster2 = clusterPairs2run[i_row, 1]
            cellsfrom1 = [idx for idx in range(len(labels)) if labels['labels'][idx] == cluster1]
            cellsfrom2 = [idx for idx in range(len(labels)) if labels['labels'][idx] == cluster2]
            sub_connects_mt = np.array(connects_mt)[cellsfrom1, :][:, cellsfrom2]

            if cluster1 == cluster2:
                sub_connects_mt[np.tril_indices(sub_connects_mt.shape[0], k=-1)] = False
            indices = np.where(sub_connects_mt)
            cell1 = np.array([cellsfrom1[index[0]] for index in list(zip(*indices))], dtype=int)
            cell2 = np.array([cellsfrom2[index[1]] for index in list(zip(*indices))], dtype=int)
            cell1 = np.array(labels.loc[cell1,'cells'])
            cell2 = np.array(labels.loc[cell2,'cells'] )
            detailed_connections = pd.concat([detailed_connections,pd.DataFrame({"cell1":cell1,
                                                "cell2":cell2,
                                                "Type":cluster1 + "---" + cluster2})])

        # Calculate p-values and adjusted q-values
        K = (counts.sum().sum() + np.sum(diag))/2
        N = sum(cellCounts.values()) * (sum(cellCounts.values()) - 1) / 2
        p_value = np.zeros((len(standards), len(standards)))
        for i in range(p_value.shape[0]):
            for j in range(p_value.shape[1]):
                if i == j:
                    M = int(cellCounts[standards[i]] * (cellCounts[standards[j]]) / 2)
                else:
                    M = cellCounts[standards[i]] * cellCounts[standards[j]]
                p_value[i, j] = 1-hypergeom.cdf(counts.iloc[i, j], N, M, K) 

        p_value = pd.DataFrame(p_value)
        p_value.index,p_value.columns = counts.index,counts.columns

        # Post-processing
        clusters = pd.Series(labels['labels']).unique()
        row_col_idx = np.triu_indices(len(clusters), k=0)
        q_value = np.array(p_value)
        tmp = multipletests(q_value[row_col_idx], method='fdr_bh')[1]
        q_value[row_col_idx] = tmp
        q_value = q_value.T
        q_value[row_col_idx] = tmp
        q_value = pd.DataFrame(q_value)
        q_value.index,q_value.columns = counts.index,counts.columns
        
        result = {}
        result["connections"] = counts
        result["pvalue"] = p_value
        result['qvalue'] = q_value
        result['detailed_connections'] = detailed_connections
        result['topK'] = topK
        result['euclidean_distance'] = euclidean_distance
        result['rank_distance'] = rank_distance/len(coordinates)
        return result

    def optimization(self, affinityMat, initial_config=None, k=3, max_iter=1000, min_cost=0, condition="tight"):
        """
        Embedding cells into three dimensional space.

        :param affinityMat: Affinity matrix.
        :param initial_config: Initial coordinates.
        :param k: Dimensions of the projection space, only support k=3 for now.
        :param max_iter: Maximum number of iterations.
        :param min_cost: Minimum cost.
        :param condition: Condition for optimization.
        :return: Optimized coordinates.
        """
        
        n = affinityMat.shape[0]
        momentum = 0.5  # initial momentum
        final_momentum = 0.8  # final momentum
        mom_switch_iter = 250  # value to which momentum is changed
        epsilon = 1000  # initial learning rate
        min_gain = 0.01  # minimum gain for delta-bar-delta
        eps = 2 ** (-52)
        epoch = 100
        
        if initial_config is not None and isinstance(initial_config, np.ndarray):
            # if initial_config.shape[0] != n or initial_config.shape[1] != k:
            if  initial_config.shape[1] != k:
                raise ValueError("initial_config argument does not match necessary configuration for X")
            ydata = np.random.normal(size=(n, k))+50
            ydata[:len(initial_config)] = initial_config
        else:
            ydata = np.random.normal(size=(n, k))
        
        P = 0.5 * (affinityMat + affinityMat.T)
        P[P < eps] = eps
        P /= P.sum()
        
        grads = np.zeros_like(ydata)
        incs = np.zeros_like(ydata)
        gains = np.ones_like(ydata)
        
        for iter in range(1,max_iter+1):
            sum_ydata = np.apply_along_axis(lambda row: np.sum(row ** 2), axis=1, arr=ydata)
            d =  sum_ydata + (-2 * np.dot(ydata, ydata.T)) + sum_ydata.reshape((-1, 1))

            num = 1 / (1 + d)
            np.fill_diagonal(num, 0)
            Q = num / num.sum()
            
            if np.any(np.isnan(num)):
                print("NaN in grad. descent")
                
            Q[Q < eps] = eps
            P_Q = P - Q
            P_Q[(P_Q > 0) & (d <= 0.01)] = -0.01  # distance restraining
            #stiffnesses = 4 * (P_Q) * num # distance restraining
            stiffnesses = 4 * (P - Q) * num # distance restraining
            
            for i in range(n):
                grads[i, :] = np.sum((-ydata + ydata[i, :]) * stiffnesses[:, i].reshape(-1, 1), axis=0)
            
            gains = ((gains + 0.2) * np.abs(np.sign(grads) != np.sign(incs)) + 
                    gains * 0.8 * np.abs(np.sign(grads) == np.sign(incs)))
            gains[gains < min_gain] = min_gain
            
            incs = momentum * incs - epsilon * (gains * grads)
            ydata += incs
            ydata = ydata - ydata.mean(axis=0, keepdims=True)
            
            if iter == mom_switch_iter:
                momentum = final_momentum
            
            if iter % epoch == 0:
                cost = np.sum(P * np.log((P + eps) / (Q + eps)))
                print(f"Iteration #{iter} loss function cost is: {cost}")
                
                if cost < min_cost:
                    break
                    
            range_ = np.max(np.abs(ydata))
            if condition == "tight":
                if range_ > 50 and iter % 10 == 0:
                    ydata *= 50 / range_
            else:
                if range_ > 50 and iter % max_iter == 0:
                    ydata *= 50 / range_
        
        return ydata

    def getContribution(self, TPM, LR, detailed_connections, verbose=True):
        """
        Calculate the contribution of ligand-receptor pairs.

        :param TPM: TPM matrix.
        :param LR: Ligand-receptor pairs.
        :param detailed_connections: Detailed connections between cells.
        :param verbose: Print verbose information.
        :return: Contribution list of ligand-receptor pairs.
        """
        LR.loc[:, [0,1]] =  LR.loc[:, [0,1]].astype(str)
        if verbose:
            print("Extracting data matrix")
        ligands_existed = np.intersect1d(TPM.index, LR.loc[:,0])
        receptors_existed = np.intersect1d(TPM.index, LR.loc[:,1])
        
        flt_LR = LR.loc[(np.isin(LR.loc[:, 0], ligands_existed)) & (np.isin(LR.loc[:, 1], receptors_existed)), :]
        
        if LR.shape[1] > 2:
            LRscores = flt_LR.loc[:, 2]
        else:
            LRscores = np.repeat(1, LR.shape[0])
        
        if verbose:
            print(f"Calculate contribution of {len(np.unique(detailed_connections['Type']))} cluster pairs.")
        
        LR_contri_lst = pd.DataFrame()
        
        for target_clusterPair in np.unique(detailed_connections['Type']):
            if detailed_connections[detailed_connections['Type']==target_clusterPair].shape[0] < 3:
                print(f"Number of connected cells in {target_clusterPair} is lower than 3.")
                continue
            
            L1S = TPM.loc[flt_LR.loc[:, 0], detailed_connections[detailed_connections['Type']==target_clusterPair].iloc[:, 0]]
            R1S = TPM.loc[flt_LR.loc[:, 1], detailed_connections[detailed_connections['Type']==target_clusterPair].iloc[:, 0]]
            L1R = TPM.loc[flt_LR.loc[:, 0], detailed_connections[detailed_connections['Type']==target_clusterPair].iloc[:, 1]]
            R1R = TPM.loc[flt_LR.loc[:, 1], detailed_connections[detailed_connections['Type']==target_clusterPair].iloc[:, 1]]
            
            all_intensity = np.array(L1S)  * np.array(R1R) * LRscores.values[:, np.newaxis].repeat(L1S.shape[1], axis=1) + np.array(R1S) * np.array(L1R) * LRscores.values[:, np.newaxis].repeat(R1S.shape[1], axis=1)
            all_intensity = pd.DataFrame(all_intensity)
            all_intensity.index = [f"{flt_LR.iloc[i, 0]}---{flt_LR.iloc[i, 1]}" for i in range(len(flt_LR))]
            contribution_mt = all_intensity.div(all_intensity.sum(axis=0), axis=1)
            
            contribution_forCluster = pd.DataFrame(contribution_mt.sum(axis=1).div(contribution_mt.shape[1]).sort_values(ascending=False),columns=['Contribution'])
            contribution_forCluster['Type'] = target_clusterPair
            LR_contri_lst = pd.concat([LR_contri_lst,contribution_forCluster])

            del L1S, L1R, R1S, R1R, all_intensity, contribution_mt, contribution_forCluster
            gc.collect()

        return LR_contri_lst

    def getDensity3D(self, x, y, z, n=100, **kwargs):
        """
        Calculate 3D density.

        :param x: X coordinates.
        :param y: Y coordinates.
        :param z: Z coordinates.
        :param n: Number of grid points.
        :return: Density values.
        """
        try:
            #dens = gaussian_kde(np.vstack([x, y, z]))
            dens = self.kde3d(x, y, z,n = n)
        except Exception as e:
            print(e)
            print("Switch bandwidth to h = 1")
            dens = self.kde2d(np.vstack([x, y]), bw_method=1)
        
        ix = np.searchsorted(dens['x'], x, side='left')
        iy = np.searchsorted(dens['y'], y, side='left')
        iz = np.searchsorted(dens['z'], z, side='left')
        ii = np.vstack([ix, iy, iz]).T
        dens_d = [dens["d"][ii[i,0]-1,ii[i,1]-1,ii[i,2]-1] for i in range(len(ii))]
        return dens_d

    def bandwidth_nrd(self, v):
        """
        Calculate the normal reference distribution bandwidth.

        :param v: Data vector.
        :return: Bandwidth value.
        """
        r = np.quantile(v, [0.25, 0.75])
        h = (r[1] - r[0])/1.34
        h =  4 * 1.06 * min(np.sqrt(np.var(v, ddof=1)), h) * len(v)**(-1/5)  # R的MASS::bandwidth.nrd(X)
        return h

    def kde3d(self, x, y, z, h=None, n=20, lims=None):
        """
        Calculate 3D kernel density estimation.

        :param x: X coordinates.
        :param y: Y coordinates.
        :param z: Z coordinates.
        :param h: Bandwidth.
        :param n: Number of grid points.
        :param lims: Limits for the grid.
        :return: KDE results.
        """
        nx = len(x)
        if len(y) != nx or len(z) != nx:
            raise ValueError("data vectors must be the same length")
        
        if h is None:
            h = np.array([self.bandwidth_nrd(x), self.bandwidth_nrd(y), self.bandwidth_nrd(z)]) / 6
    
        elif len(h) != 3:
            h = np.repeat(h, 3)
        else:
            h = np.array(h)
            
        if len([n]) != 3:
            n = [n]*3
        
        if lims is None:
            lims = [[min(x), max(x)], [min(y), max(y)], [min(z), max(z)]]
        elif len(lims) == 2:
            lims = np.repeat(lims, 3, axis=0)
        else:
            lims = np.array(lims)

        gx = np.linspace(lims[0][0], lims[0][1], num=n[0])
        gy = np.linspace(lims[1][0], lims[1][1], num=n[1])
        gz = np.linspace(lims[2][0], lims[2][1], num=n[2])

        mx = np.transpose(
            [np.exp(-(x_i - gx)**2 / (2*h[0]**2)) for x_i in x]
        )
        my = np.transpose(
            [np.exp(-(y_i - gy)**2 / (2*h[1]**2)) for y_i in y]
        )
        mz = np.transpose(
            [np.exp(-(z_i - gz)**2 / (2*h[2]**2)) for z_i in z]
        )

        v = np.zeros(n)
        tmy_nx = my.T / nx
        for k in range(n[2]):
            tmy_nz_zk = tmy_nx.T * mz[ k,:]
            v[:, :, k] = mx.dot(tmy_nz_zk.T)

        return {"x": gx, "y": gy, "z": gz, "d": v}


    def kde2d(self,x, y, h=None, n=25, lims=None):
        nx = len(x)
        if len(y) != nx:
            raise ValueError("data vectors must be the same length")
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError("missing or infinite values in the data are not allowed")
        if lims is None:
            lims = [x.min(), x.max(), y.min(), y.max()]
        if not np.isfinite(lims).all():
            raise ValueError("only finite values are allowed in 'lims'")
        n_x, n_y = n, n
        if isinstance(n, int):
            n_x, n_y = n, n
        elif isinstance(n, tuple):
            n_x, n_y = n
        else:
            raise ValueError("'n' must be an integer or a tuple")
        gx = np.linspace(lims[0], lims[1], num=n_x)
        gy = np.linspace(lims[2], lims[3], num=n_y)
        if h is None:
            h_x, h_y = [self.bw_nrd(x), self.bw_nrd(y)]
        else:
            h_x, h_y = h if isinstance(h, tuple) else [h, h]
        if h_x <= 0 or h_y <= 0:
            raise ValueError("bandwidths must be strictly positive")
        h_x, h_y = h_x/4, h_y/4
        ax = np.outer(gx, x) - x.mean()
        ay = np.outer(gy, y) - y.mean()
        den = norm.pdf(ax/h_x).reshape(ax.shape[0], -1) @ norm.pdf(ay/h_y).reshape(ay.shape[0], -1)
        z = den / (nx*h_x*h_y)
        return {'x': gx, 'y': gy, 'z': z}

    def bw_nrd(self,x):
        """
        Calculate the normal reference distribution bandwidth.

        :param x: Data vector.
        :return: Bandwidth value.
        """
        IQR = np.percentile(x, q=75) - np.percentile(x, q=25)
        h = 0.9 * min(x.std(), IQR/1.34) / len(x)**0.2
        return h

    def loginfo(self,*args, printnow=True):
        """
        Log information with a timestamp.

        :param args: Arguments to log.
        :param printnow: Print the log immediately.
        :return: Log message.
        """
        msg = "".join(map(str, args))
        msg = "[{}]: {}\n".format(pd.Timestamp.now(), msg)
        if printnow:
            print(msg, end='')
        return msg

    #                    _ooOoo_
    #                   o8888888o
    #                   88" . "88
    #                   (| -_- |)
    #                   O\  =  /O
    #                ____/`---'\____
    #              .'  \\|     |//  `.
    #             /  \\|||  :  |||//  \
    #            /  _||||| -:- |||||-  \
    #            |   | \\\  -  /// |   |
    #            | \_|  ''\---/''  |   |
    #            \  .-\__  `-`  ___/-. /
    #          ___`. .'  /--.--\  `. . __
    #       ."" '<  `.___\_<|>_/___.'  >'"".
    #      | | :  `- \`.;`\ _ /`;.`/ - ` : | |
    #      \  \ `-.   \_ __\ /__ _/   .-` /  /
    # ======`-.____`-.___\_____/___.-`____.-'======
    #                    `=---='
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #  Blessed by the Buddha, may there be no bugs.

    def RedeTIL_plot(self,cellinfo_tbl,signif_results,outdir):
        """
        Plot results using RedeTIL.

        :param cellinfo_tbl: Cell information table.
        :param signif_results: Significance results.
        :param outdir: Output directory.
        """
        from matplotlib import pyplot as plt
        import seaborn as sns
        import plotly.express as px
        sns.set_theme(style="white")
        
        # 3D plot
        fig = px.scatter_3d(cellinfo_tbl, x='x', y='y', z='z',
                color='labels',hover_name ='cells' )
        fig.show()
        
        colors = dict((fig.data[i].name, fig.data[i].marker.color) for i in range(len(fig.data))) 
        # 2D section plot
        data_tmp = cellinfo_tbl[(cellinfo_tbl['z']<2)&(cellinfo_tbl['z']>-2)]
        plt.figure(figsize=[4,4])
        sns.scatterplot(data = data_tmp
                            ,hue='labels'
                            ,x='x',y='y',palette = colors
                            ,s=30 ,alpha=0.7)
        plt.xlim([-55,55])
        plt.ylim([-55,55])
        plt.legend(loc=[1.01,0.2],framealpha=0)
        plt.title('Section plot \nz = 0')
        plt.savefig(f'{outdir}/section.png',bbox_inches='tight',facecolor='w')
        plt.show()
        
        # 2D section density plot
        # data_tmp['-log10 density']=-np.log10(data_tmp['density'])
        # plt.figure(figsize=[5,4])
        # g= sns.scatterplot(data = data_tmp
        #                 ,hue='-log10 density'
        #                     ,x='x',y='y',palette = 'autumn_r'
        #                     ,s=30,alpha=0.9)
        # plt.legend(loc=[1.01,0.2],framealpha=0)
        # plt.title('Section plot \nz = 0')
        # norm = plt.Normalize( data_tmp['-log10 density'].max(),data_tmp['-log10 density'].min()) 
        # sm = plt.cm.ScalarMappable(cmap="autumn", norm=norm)
        # sm.set_array([])
        # # Remove the legend and add a colorbar
        # g.get_legend().remove()
        # g.figure.colorbar(sm,label='Density') 
        # plt.xlim([-55,55])
        # plt.ylim([-55,55])
        # plt.savefig(f'{outdir}/section_Density.png',bbox_inches='tight',facecolor='w')
        # plt.show()
        
        # Connection plot
        signif_results['connections'].columns.name,signif_results['qvalue'].columns.name = 'from','from'
        signif_results['connections'].index.name,signif_results['qvalue'].index.name = 'to','to'
        conn_mat = signif_results['connections'].stack()
        q_mat =  signif_results['qvalue'].stack()
        conn_mat = pd.concat([conn_mat,q_mat],axis=1).reset_index()
        conn_mat['type'] = 'other'
        conn_mat.loc[conn_mat[1]<0.05,'type'] = 'enriched'
        conn_mat.loc[conn_mat[1]>0.95,'type'] = 'depleted'
        conn_mat.rename(columns={0:'Connection',1:'qvalue'},inplace=True)
        g = sns.relplot(
                data=conn_mat,
                x="from", y="to", hue="type", hue_order=['enriched','other','depleted'],
                size='Connection',palette=['r','g','b'], height=5, 
                sizes=(100, 400), size_norm=(min(conn_mat['Connection']), max(conn_mat['Connection'])),
            )
        plt.xticks(rotation=90)
        plt.xlabel(None)
        plt.ylabel(None)
        plt.title('Connection map')
        plt.savefig(f'{outdir}/Connection.png',bbox_inches='tight',facecolor='w')
        plt.show()
        
        # Rank Distance
        signif_results['rank_distance'].columns.name,signif_results['rank_distance'].index.name = 'from','to'
        g = sns.heatmap(
                data=signif_results['rank_distance'],cmap='mako',square=True,annot=True
            )
        plt.xticks(rotation=90)
        plt.xlabel('From')
        plt.ylabel('To')
        plt.title('Rank distance map')
        plt.savefig(f'{outdir}/RankDistance_map.png',bbox_inches='tight',facecolor='w')
        plt.show()
        
        # Distance
        g = sns.heatmap(
                data=signif_results['euclidean_distance'],cmap='mako',square=True,annot=True
            )
        plt.xticks(rotation=90)
        plt.title('Distance map')
        plt.xlabel('')
        plt.ylabel('')
        plt.savefig(f'{outdir}/Distance_map.png',bbox_inches='tight',facecolor='w')
        plt.show()
        
    def get_LR(self):
        """
        Load ligand-receptor pairs from a file.

        :return: Ligand-receptor pairs.
        """
        current_dir = os.path.dirname(__file__)
        LR_path = os.path.join(current_dir, 'data/LR_pairs.txt')
        LR = pd.read_csv(LR_path, sep='\t', header=None)
        return LR

    def TPM_LR_dot(self, TPM,LR):
        """
        Calculate the affinity matrix for ligand-receptor pairs.

        :param TPM: TPM matrix.
        :param LR: Ligand-receptor pairs.
        :return: Affinity matrix.
        """
        genenames = TPM.index
        TPM = TPM.fillna(0)
        # find ligands and receptors TPM
        ligandsIndex = [np.where(genenames == x)[0][0] if x in genenames else np.nan for x in pd.concat([LR.loc[:, 0], LR.loc[:, 1]], axis=0).values]
        receptorIndex = [np.where(genenames == x)[0][0] if x in genenames else np.nan for x in pd.concat([LR.loc[:, 1], LR.loc[:, 0]], axis=0).values]
        validIndex = np.where(~np.isnan(ligandsIndex) & ~np.isnan(receptorIndex))[0]
        ligandsTPM = TPM.iloc[np.array(ligandsIndex)[validIndex], :].values
        receptorTPM = TPM.iloc[np.array(receptorIndex)[validIndex], :].values
        LRscores = pd.concat([LR.iloc[:, 2],LR.iloc[:, 2]]).reset_index()
        LRscores = LRscores.loc[validIndex,2]
        affinityMat = np.dot(np.dot(ligandsTPM.T, np.diag(LRscores)), receptorTPM)
        return affinityMat

    def get_affinity_matrix_naive(self,TPM,LR):
        """
        Calculate the affinity matrix for ligand-receptor pairs.

        :param TPM: TPM matrix.
        :param LR: Ligand-receptor pairs.
        :return: Affinity matrix.
        """
        affinityMat = self.TPM_LR_dot(TPM,LR)
        
        for i in range(affinityMat.shape[0]):
            affinityArray = affinityMat[i, :]
            affinityArray[i] = 0
            affinityArraySorted = np.sort(affinityArray)[::-1]
            affinityArray[affinityArray <= affinityArraySorted[50]] = 0 # Denoise
            affinityMat[i, :] = affinityArray
        return affinityMat

    def Abundance_features(self):
        """
        Calculate abundance features.

        :return: Abundance features.
        """
        cts = np.unique(self.label['labels'])
        cell_ratio = pd.DataFrame()
        for ct in cts:
            label_T = self.label[self.label['labels']==ct]
            data_T = self.TPM.loc[self.genesets,label_T['cells']].T>0
            tmp = 100*data_T.sum()/len(label_T)
            tmp = pd.DataFrame(tmp,
                            columns = [f'Frequency in {ct} (%)'])
            cell_ratio = pd.concat([cell_ratio,tmp],axis=1)
        cell_ratio.index = [cell_ratio.index[i]+'+' for i in range(len(cell_ratio))]
        os.makedirs(f'{self.outdir}/Abundance_features',exist_ok=True)
        # cell_ratio.to_csv(f'./results/Abundance_features/Abundance_features.txt',sep='\t')
        cell_ratio.to_csv(f'{self.outdir}/Abundance_features/Abundance_features.txt', sep='\t')
        return cell_ratio

    def Spatial_features(self, plot=False):
        """
        Calculate spatial features.

        :param plot: Plot results.
        :return: Spatial features.
        """
        LR = self.get_LR()
        os.makedirs(f'{self.outdir}/Spatial_features',exist_ok=True)
        affinityMat = self.get_affinity_matrix_naive(self.TPM,LR)
        # pd.DataFrame(affinityMat,index=TPM.columns,columns=TPM.columns).to_csv(f"{outdir}/Spatial_features/affinityMat.txt",  sep="\t")
        # optimization
        coords = self.optimization(affinityMat)
        coords = pd.DataFrame(coords,columns = ['x', 'y', 'z'])
        coords.index = self.TPM.columns
        cellinfo_tbl = coords.merge(self.label,left_on =coords.index,right_on='cells')
        density_obj = self.getDensity3D(cellinfo_tbl['x'], cellinfo_tbl['y'], cellinfo_tbl['z'])
        cellinfo_tbl = cellinfo_tbl.assign(density=density_obj)
        signif_results = self.getSignificance(coords.values, labels=cellinfo_tbl[['cells','labels']])
        contribution_list = self.getContribution(self.TPM, LR, signif_results['detailed_connections'])
        cellinfo_tbl.rename(columns={'cells':'cellName'}).to_csv(f'{self.outdir}/Spatial_features/cellinfo_tbl.txt', sep='\t', index=False)

        signif_results['connections'].to_csv(f'{self.outdir}/Spatial_features/Cell-Cell_connections.txt', sep='\t')
        signif_results['qvalue'].to_csv(f'{self.outdir}/Spatial_features/Cell-Cell_connections_qvalue.txt', sep='\t')
        signif_results['detailed_connections'].to_csv(f'{self.outdir}/Spatial_features/Detailed_connections.txt', sep='\t', index=False)
        signif_results['rank_distance'].to_csv(f'{self.outdir}/Spatial_features/Rank Cell-Cell_Distance.txt', sep='\t')
        signif_results['euclidean_distance'].to_csv(f'{self.outdir}/Spatial_features/Euclidean Cell-Cell_Distance.txt', sep='\t')
        contribution_list.to_csv(f'{self.outdir}/Spatial_features/Ligand-receptor_contribution.txt', sep='\t')
        if plot==True:
            self.RedeTIL_plot(cellinfo_tbl,signif_results,self.outdir)
            
        return signif_results
    
    def Dynamic_features(self, plot=False):
        """
        Calculate dynamic features.

        :param plot: Plot results.
        :return: Dynamic features.
        """
        LR = self.get_LR()
        if self.target is None:
            raise ValueError("Drug target was not provided")
        if self.target not in self.TPM.index:
            raise ValueError("Drug target was not detected in single-cell profiles")
        if self.combo_target not in self.TPM.index:
            raise ValueError("Combo target was not detected in single-cell profiles")
        if self.T_cells is None:
            raise ValueError("Identity of T cells was not provided")
        if self.Cancer_cells is None:
            raise ValueError("Identity of Cancer cells was not provided")
        if len(self.label[self.label['labels'] == self.T_cells])==0:
            raise ValueError(f"{self.T_cells} was not detected")
        if len(self.label[self.label['labels'] == self.Cancer_cells])==0:
            raise ValueError(f"{self.Cancer_cells} was not detected")
        
        if 'Spatial_features' in os.listdir(self.outdir):
            initial_config = pd.read_csv(f"{self.outdir}/Spatial_features/cellinfo_tbl.txt",sep='\t')
            if not os.path.exists(f'{self.outdir}/Spatial_features/Distance Tsub-Cancer.txt'):
                self.TILC_subset(initial_config, self.TPM, outdir = f'{self.outdir}/Spatial_features')
        else:
            self.Spatial_features()
            initial_config = pd.read_csv(f"{self.outdir}/Spatial_features/cellinfo_tbl.txt",sep='\t')
            self.TILC_subset(initial_config, self.TPM, outdir = f'{self.outdir}/Spatial_features')
        RDO_naive = pd.read_csv(f"{self.outdir}/Spatial_features/Distance Tsub-Cancer.txt",sep='\t',index_col=0)
            
        if self.perturbation == 'depletion':
            label = self.label
            TPM_target  = self.TPM.loc[self.target,label['cells']]
            target_Cells = TPM_target[TPM_target>0].index
            label_pert = label[~label['cells'].isin(target_Cells)]
            if self.combo_target is not None:
                TPM_target  = self.TPM.loc[self.combo_target,label['cells']]
                target_Cells = TPM_target[TPM_target>0].index
                label_pert = label_pert[~label_pert['cells'].isin(target_Cells)] 
            TPM_pert = self.TPM.loc[:,label_pert['cells']]
            affinityMat = self.get_affinity_matrix_naive(TPM_pert,LR)
        elif self.perturbation == 'block':
            label = self.label
            TPM_target  = self.TPM.loc[self.target,label['cells']]
            target_Cells = TPM_target[TPM_target>0].index
            if self.combo_target is not None:
                TPM_target  = self.TPM.loc[self.combo_target,label['cells']]
                target_Cells_combo = TPM_target[TPM_target>0].index
            label_pert = label
            TPM_pert = self.TPM.loc[:,label_pert['cells']]
            if self.combo_target is None:
                TPM_pert.loc[:,target_Cells] *=  0.7
            else:
                TPM_pert.loc[:,list(set(list(target_Cells) + list(target_Cells_combo)))] *= 0.7  
            affinityMat = self.get_affinity_matrix_naive(TPM_pert,LR)
            
        initial_config = initial_config[initial_config['cellName'].isin(label_pert['cells'])]     
        # optimization
        coords = self.optimization(affinityMat,initial_config=initial_config)
        coords = pd.DataFrame(coords,columns = ['x', 'y', 'z'])
        coords.index = TPM_pert.columns
        cellinfo_tbl = coords.merge(label_pert,left_on =coords.index,right_on='cells')
        density_obj = self.getDensity3D(cellinfo_tbl['x'], cellinfo_tbl['y'], cellinfo_tbl['z'])
        cellinfo_tbl = cellinfo_tbl.assign(density=density_obj)
        signif_results = self.getSignificance(coords.values, labels=cellinfo_tbl[['cells','labels']])
        contribution_list = self.getContribution(TPM_pert, LR, signif_results['detailed_connections'])
        cellinfo_tbl.rename(columns={'cells':'cellName'},inplace=True)
        
        if self.combo_target is None:
            out_path = f'{self.outdir}/Dynamic_features/{self.target}_{self.perturbation}'
        else:
            out_path = f'{self.outdir}/Dynamic_features/{self.target}+{self.combo_target}_{self.perturbation}'
        os.makedirs(f'{out_path}/Perturbated Spatial_features',exist_ok=True)
          
        cellinfo_tbl.to_csv(f'{out_path}/Perturbated Spatial_features/cellinfo_tbl.txt', sep='\t', index=False)
        signif_results['connections'].to_csv(f'{out_path}/Perturbated Spatial_features/Cell-Cell_connections.txt', sep='\t')
        signif_results['qvalue'].to_csv(f'{out_path}/Perturbated Spatial_features/Cell-Cell_connections_qvalue.txt', sep='\t')
        signif_results['detailed_connections'].to_csv(f'{out_path}/Perturbated Spatial_features/Detailed_connections.txt', sep='\t', index=False)
        signif_results['rank_distance'].to_csv(f'{out_path}/Perturbated Spatial_features/Rank Cell-Cell_Distance.txt', sep='\t')
        signif_results['euclidean_distance'].to_csv(f'{out_path}/Perturbated Spatial_features/Euclidean Cell-Cell_Distance.txt', sep='\t')
        contribution_list.to_csv(f'{out_path}/Perturbated Spatial_features/Ligand-receptor_contribution.txt', sep='\t')
        if plot==True:
            self.RedeTIL_plot(cellinfo_tbl,signif_results,self.outdir)
        self.TILC_subset(cellinfo_tbl, TPM_pert, outdir = f'{out_path}/')
        
        RDO_pert = pd.read_csv(f"{out_path}/Distance Tsub-Cancer.txt",sep='\t',index_col=0)
        RDO = RDO_naive-RDO_pert
        RDO.columns = ['Infiltration change']
        RDO.to_csv(f'{out_path}/Dynamic_features.txt',sep='\t')
        return signif_results
        
    def TILC_subset(self, cellinfo_tbl, TPM,  outdir):
        """
        Calculate rank distance between cancer cell and T cell subsets.

        :param cellinfo_tbl: Cell information (cell type and coordinate) table.
        :param TPM: TPM matrix.
        :param outdir: Output directory.
        """
        relabeled_cell = self.T_cells
        TPM_marker = TPM.loc[self.genesets,:].T
        label_naive = cellinfo_tbl
        label_naive = label_naive.merge(TPM_marker,left_on='cellName',right_on=TPM_marker.index)
        coordinates = label_naive[['x','y','z']].values
        dist_mt = distance.cdist(coordinates, coordinates, 'euclidean')
        #Mean rank distance
        out = pd.DataFrame()
        for gene in self.genesets:
            if gene in label_naive.columns:
                rank_distance = pd.DataFrame(dist_mt).rank()-1
                if gene == self.target:
                    label_naive.loc[(label_naive[gene]>0),'TMP'] = f"{gene}+"
                    label_naive.loc[(label_naive[gene]==0),'TMP'] = f"{gene}-"
                    label_naive.loc[:,gene] = label_naive.loc[:,'TMP'] 
                else:
                    label_naive['TMP'] = 'others'
                    label_naive.loc[(label_naive[gene]>0)&(label_naive[self.target]==0),'TMP'] = f"{gene}+"
                    label_naive.loc[(label_naive[gene]==0),'TMP'] = f"{gene}-"
                    label_naive.loc[~(label_naive['TMP'].isin([f"{gene}+",f"{gene}-"])),'TMP'] = "other"
                    label_naive.loc[:,gene] = label_naive.loc[:,'TMP'] 
                    
                label_naive['labels2'] = label_naive['labels'] 
                label_naive.loc[label_naive['labels'].isin([relabeled_cell]),'labels2'] = label_naive.loc[label_naive['labels'].isin([relabeled_cell]),gene]+label_naive.loc[label_naive['labels'].isin([relabeled_cell]),'labels2']
                
                rank_distance.index,rank_distance.columns = "To-"+label_naive['labels2'],label_naive['labels2']
                rank_distance = rank_distance.groupby(rank_distance.index).mean()
                rank_distance = rank_distance.T.groupby(rank_distance.T.index).mean().T #Mean rank distance
                
                signif_results_rm_PD1all_RDC = rank_distance/len(self.TPM.columns)
                if f'To-{gene}+{relabeled_cell}' in signif_results_rm_PD1all_RDC.index:
                    out.loc[f'Cancer to {gene}+ T cells','Rank Distance'] = signif_results_rm_PD1all_RDC.loc[f'To-{gene}+{relabeled_cell}',self.Cancer_cells]
        out.to_csv(f"{outdir}/Distance Tsub-Cancer.txt",sep='\t')

# demo
if __name__ == '__main__':
    adata = sc.read_h5ad('./demo/demo.h5ad')

    # One target perturbation
    redetil = RedeTIL_Features(adata, target = 'PDCD1', perturbation='block',
                    T_cells='T-cell', Cancer_cells='malignant')
    redetil.Abundance_features()
    redetil.Spatial_features()
    redetil.Dynamic_features()

    # Targets combo
    redetil = RedeTIL_Features(adata, target = 'CTLA4', combo_target = 'VEGFA', perturbation='block',
                    T_cells='T-cell', Cancer_cells='malignant')
    redetil.Dynamic_features()


