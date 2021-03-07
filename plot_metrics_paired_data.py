import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 
from sklearn.preprocessing import MinMaxScaler
def plot_accuracy_per_cell_type(adata,
                               accuracy,
                               cell_type,
                               display_value=True,
                                show=True,
                               save=None):
    labels = sorted(accuracy.keys())
    colors= adata.uns[cell_type+'_colors']
    values = []
    for cell_type in labels:
        values.append(round(accuracy[cell_type], ndigits=2))
    x = np.arange(len(labels))

    f, ax = plt.subplots(figsize=(18,9)) # set the size that you'd like (width, height)

    value_plot = plt.bar(x, values, color=colors)
    # Create names on the x-axis
    #ax.legend(fontsize = 14)
    plt.xticks(x, labels, rotation=90)
    ax.set_ylabel('Percentage of matching barcodes clustered together')
    ax.set_title('Percentage of matching barcodes per cell type')
    if display_value == True:
        #"""Attach a text label above each bar in *rects*, displaying its height."""
        for rect in value_plot:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
    
    if save != None:
        plt.savefig(save, bbox_inches='tight')
    if show == True:
        plt.show()
    else:
        plt.close()
        
def metric_heatmap(dataframe, cmap='viridis', scale=False, display_values=True, show=True, save=None):
    """
    for raw we use 'viridis'
    for scaled we use 'plasma'
    """
    if scale== True:
        # scaling
        scaler = MinMaxScaler()
        index_save = dataframe.index.copy()
        dataframe = pd.DataFrame(scaler.fit_transform(dataframe), columns=dataframe.columns)
        dataframe.index = index_save.copy()
        
    # plot
    sns.set_style('ticks')
    fig, ax = plt.subplots()
    # the size of A4 paper
    fig.set_size_inches(11.7, 8.27)
    sns.heatmap(dataframe.transpose(copy=True),
                cmap=cmap,
                linewidths=0.2,
                annot=display_values,
                xticklabels=False,
                yticklabels=False,
                square=True)
    plt.yticks(np.arange(0.5, len(dataframe.columns), 1), dataframe.columns)
    plt.xticks(np.arange(0.5, len(dataframe.index), 1), dataframe.index, rotation=90)
    
    if save != None:
        fig.savefig(save, bbox_inches="tight")
        
    if show== True:
        plt.show()
    else:
        plt.close()
    
