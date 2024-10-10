from scipy.spatial.distance import pdist, squareform
import pandas as pd
import numpy as np
from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
import os
def distance_table(df1,df2,metric="euclidean"):
    df = pd.concat([df1,df2])
    x = squareform(pdist(df.values,metric=metric))
    x = np.nan_to_num(x)
    m = x.max()
    x = x[:df1.shape[0],:]
    x = x[:,df1.shape[0]:]

    return pd.DataFrame(x,index=df1.index.values,columns=df2.index.values), m

def plot_loadings_components(pca, component1, component2, features, n_genes, ax):
    """
    Function to plot the loadings of two specified PCA components on the same 2D plot
    and select the top genes based on the square mean of the loadings over all components.

    Parameters:
        pca (PCA): Fitted PCA object.
        component1 (int): Index of the first component to plot.
        component2 (int): Index of the second component to plot.
        features (list): List of feature names.
        n_genes (int): Number of top genes to select.

    Returns:
        None
    """
    # Calculate the mean square of loadings for each gene across all components
    loadings_sq_mean = np.mean(pca.components_ ** 2, axis=0)
    # Get the indices of the top genes based on the square mean
    top_gene_indices = np.argsort(loadings_sq_mean)[-n_genes:]

    # Get the loadings for the specified components
    loadings_comp1 = pca.components_[component1, top_gene_indices]
    loadings_comp2 = pca.components_[component2, top_gene_indices]

    # Calculate the minimum and maximum values of the loadings for both axes
    min_loading_x = min(loadings_comp1)
    max_loading_x = max(loadings_comp1)
    min_loading_y = min(loadings_comp2)
    max_loading_y = max(loadings_comp2)

    # Add a small buffer to the limits
    buffer_x = 0.1 * (max_loading_x - min_loading_x)
    min_loading_x -= buffer_x
    max_loading_x += buffer_x

    buffer_y = 0.1 * (max_loading_y - min_loading_y)
    min_loading_y -= buffer_y
    max_loading_y += buffer_y

    ax.set_xlim(min_loading_x, max_loading_x)  # Set x-axis limits
    ax.set_ylim(min_loading_y, max_loading_y)  # Set y-axis limits
    
    ax.quiver(np.zeros(n_genes), np.zeros(n_genes), loadings_comp1, loadings_comp2,
              angles='xy', scale_units='xy', scale=1, color='b', width=0.005, label=f'Components {component1+1}-{component2+1}')

    # Annotate the top genes based on loadings
    texts = []
    for idx in range(n_genes):
        x_pos = loadings_comp1[idx]
        y_pos = loadings_comp2[idx]
        text = features[top_gene_indices[idx]]

        texts.append(ax.text(x_pos, y_pos, text, fontsize=10, ha='left', va='bottom', color='black', bbox=dict(facecolor='white', alpha=0.8)))

    # Adjust text positions to avoid overlap with quiver vectors
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))

    ax.set_xlabel(f'Component {component1+1} Loading')
    ax.set_ylabel(f'Component {component2+1} Loading')
    ax.set_title(f'Top {n_genes} Genes based on Loadings of PCA Components {component1+1} and {component2+1}')
    ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    ax.axvline(0, color='k', linestyle='--', linewidth=0.5)
    ax.grid(False)
    ax.legend()

    return ax

def plot_volcano(res_df: pd.DataFrame = None, 
                 lfc_threshold: float = 1, pvalue_threshold: float = 0.05, basemean_filter: float = 10,
                markers: list = None, gene_labels:bool=True, dot_size:int=20,
                palette1: list = ['red', 'green', 'blue', 'lightgrey'], palette2: list = ['purple'], title: str = '',
                get_degs:bool=False, save:bool=True,save_format: str = 'svg'):
    """
    Generates a volcano plot for visualizing differentially expressed genes (DEGs).

    Parameters:
    - res_df (pd.DataFrame): DataFrame containing the gene expression results with columns 'log2FoldChange', 'pvalue', and 'baseMean'.
    - lfc_threshold (float): Threshold for log2 fold change to consider a gene as differentially expressed.
    - pvalue_threshold (float): P-value threshold to consider a gene as differentially expressed.
    - basemean_filter (float): Minimum mean expression level to filter out lowly expressed genes.
    - markers (list): List of specific gene symbols to highlight on the plot.
    - gene_labels (bool): Whether to add gene labels to the plot.
    - size (float): Size for markers.
    - palette1 (list): List of colors for plotting categories: ['DEG', 'Log2FC', 'p-value', 'NS'].
    - palette2 (list): List of colors for highlighting specific markers.
    - title (str): Title of the plot.
    - get_degs (bool): Whether to save the list of differentially expressed genes (DEGs) to a CSV file.
    - save_format (str): File format for saving the plot (e.g., 'svg', 'png').

    Returns:
    - fig, ax: Matplotlib figure and axis objects of the generated plot.
    """
    
    # Check if res_df is provided
    if res_df is None:
        raise ValueError("A results DataFrame must be provided.")
        
    # Sorting the dataframe
    res = res_df.sort_values(by='log2FoldChange', ascending=False).copy()
    res['Symbol'] = res.index
    res['-log10P'] = -np.log10(res['pvalue'])
    
    # Filter lowly expressed genes
    res = res[res['baseMean'] > basemean_filter]
    
    # Define color mapping function
    def map_color(row, lfc_threshold, pvalue_threshold, markers):
        if markers and row['Symbol'] in markers:
            return 'Marker'
        if abs(row['log2FoldChange']) > lfc_threshold and row['pvalue'] < pvalue_threshold:
            return 'DEG'
        if abs(row['log2FoldChange']) > lfc_threshold:
            return 'Log2FC'
        if row['pvalue'] < pvalue_threshold:
            return 'p-value'
        return 'NS'

    # Apply color mapping
    res['Category'] = res.apply(map_color, axis=1, args=(lfc_threshold, pvalue_threshold, markers))
    # number of degs
    degs_n=res[res['Category']=='DEG'].shape[0]

    # Plotting
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    
    sns.scatterplot(data=res, x='log2FoldChange', y='-log10P', hue='Category', s=dot_size,
                    hue_order=['DEG', 'Log2FC', 'p-value', 'NS'],
                    palette=palette1, edgecolor='black', ax=ax)
    
    if markers:
        sns.scatterplot(data=res[res['Category'] == 'Marker'], x='log2FoldChange', y='-log10P', hue='Category', s=dot_size*1.5,
                        palette=palette2, edgecolor='black', ax=ax)
    
    # Axes titles and lines
    ax.set_title(f'LFC threshold: {lfc_threshold} // p-value threshold: {pvalue_threshold} // DEGs: {degs_n}/{res.shape[0]}', fontsize=10)
    ax.axhline(-np.log10(pvalue_threshold), zorder=0, color='black', linewidth=2, linestyle='--')
    ax.axvline(lfc_threshold, zorder=0, color='black', linewidth=2, linestyle='--')
    ax.axvline(-lfc_threshold, zorder=0, color='black', linewidth=2, linestyle='--')
    ax.legend().set_title('')
    
    # Add gene labels
    if gene_labels:
        if markers is not None:
            if len(markers)>0:
                m_df = res[res['Category'] == 'Marker']
                labels = []
                for l in range(m_df.shape[0]):
                    x_pos = m_df['log2FoldChange'].iloc[l]
                    y_pos = m_df['-log10P'].iloc[l]
                    label = m_df['Symbol'].iloc[l]
                    labels.append(ax.text(x_pos, y_pos, str(label), fontsize=10, ha='center', va='bottom', color='black', 
                                        bbox=dict(facecolor='white', alpha=0.5)))
                adjust_text(labels, arrowprops=dict(arrowstyle='-', color='gray'))
            else:
                raise ValueError('Empty marker list provided.')
        elif markers is None:
            m_df = res[res['Category'] == 'DEG']
            labels = []
            for l in range(m_df.shape[0]):
                x_pos = m_df['log2FoldChange'].iloc[l]
                y_pos = m_df['-log10P'].iloc[l]
                label = m_df['Symbol'].iloc[l]
                labels.append(ax.text(x_pos, y_pos, label, fontsize=10, ha='center', va='bottom', color='black', 
                                    bbox=dict(facecolor='white', alpha=0.5)))
            adjust_text(labels, arrowprops=dict(arrowstyle='-', color='gray'))
    
    # Overall title
    plt.suptitle(title, y=0.94)
    
    # Save plot
    if save:
        os.makedirs(f'results/volcano/', exist_ok=True)
        plt.savefig(f'results/volcano/volcano_{title}.{save_format}')

    if get_degs:
        degs= res[(abs(res['log2FoldChange']) > lfc_threshold) & (res['pvalue'] < pvalue_threshold)]
        degs.to_csv(f'results/volcano/volcano_degs_{title}.csv')

    return fig, ax


def create_mixed_populations(pop_sets, proportion, i1:int=None, i2:int=None):
    """
    Create a mixed population of cells from two different populations, given a proportion for each.

    Parameters:
    -----------
    pop_sets : dict
        A dictionary where the keys represent different cell populations and the values are lists of 
        AnnData objects (cell data) for each population.
    proportion : float
        The proportion of the first population (e.g., 0.7 for 70% of the first population and 30% of the second).
    i1 : int, optional
        Index for selecting a specific subset (e.g., time point or condition) within the first population in pop_sets.
    i2 : int, optional
        Index for selecting a specific subset within the second population in pop_sets.

    Returns:
    --------
    a_p_population : AnnData
        An AnnData object representing the concatenated mixed population of cells.
    """

    # Extract the keys (population names) from the dictionary pop_sets
    k = list(pop_sets.keys())
    
    # Determine the minimum number of cells between the two populations (pop1 and pop2) at the given indices i1, i2
    # Ensures that we don't sample more cells than are available from either population
    n_cells = min(pop_sets[k[0]][i1].n_obs, pop_sets[k[1]][i2].n_obs)

    # Determine how many cells will be taken from each population based on the given proportion
    n_s1 = int(proportion * n_cells)  # Number of cells from the first population
    n_s2 = n_cells - n_s1  # Remaining cells from the second population (complementary to n_s1)

    # Sample 'n_s1' cells randomly from the first population without replacement
    pop1 = pop_sets[k[0]][i1][np.random.choice(pop_sets[k[0]][i1].n_obs, n_s1, replace=False)].copy()

    # Sample 'n_s2' cells randomly from the second population without replacement
    pop2 = pop_sets[k[1]][i2][np.random.choice(pop_sets[k[1]][i2].n_obs, n_s2, replace=False)].copy()
    
    # Concatenate the sampled cells from the two populations to create the mixed population
    a_p_population = pop1.concatenate(pop2)
    
    # Return the mixed population
    return a_p_population