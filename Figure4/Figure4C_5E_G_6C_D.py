

import plot_spatial

samples = []
clust_labels_combinations = [
  []
]



for sample in samples:
  
  file_path = f"./Cell2location/data/{sample}_filtered_visium_data.h5ad"
slide = sc.read_h5ad(file_path)

for clust_labels in clust_labels_combinations:
  
  clust_col = [str(label) for label in clust_labels]

with mpl.rc_context({'figure.figsize': (10, 10)}):
  fig = plot_spatial(
    adata=slide,
    color=clust_col, 
    labels=clust_labels, 
    show_img=True,
    style='fast',  
    max_color_quantile=0.95,  
    circle_diameter=3.2,  
    colorbar_position='right',  
  )

output_dir = "./Cell2location/res/"
os.makedirs(output_dir, exist_ok=True)  
os.chdir(output_dir)

output_file = f"{sample}_spatial_{'_'.join(clust_labels)}.pdf"
plt.savefig(output_file, dpi=500, format="pdf")
print(f"Saved plot for {sample}, {clust_labels} to {output_file}")
