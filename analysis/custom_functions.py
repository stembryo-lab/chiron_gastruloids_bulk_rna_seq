from scipy.spatial.distance import pdist, squareform
import pandas as pd
import numpy as np
from adjustText import adjust_text
import matplotlib.pyplot as plt

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