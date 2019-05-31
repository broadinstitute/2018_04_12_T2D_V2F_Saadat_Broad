# Morpheus Figures

[Morpheus](https://software.broadinstitute.org/morpheus/) was developed at the Broad Institute and is a web based application.
The application can be used to visualize large datasets.
Here, we use Morpheus to visualize hierarchically clustered similarity heatmaps of our cell painting data.

## Datasets Used

We input the following `.gct` files into Morpheus:

* Plate 1 - `results/morpheus/BR00101075_normalized_variable_selected_2019_04_16_Batch1_morpheus.gct`
* Plate 2 - `results/morpheus/BR00101076_normalized_variable_selected_2019_04_16_Batch1_morpheus.gct`

These two plates were collected after 15 days.
One plate was treated with isoproterenol.

## Morpheus Parameters

We performed the following procedure to generate the morpheus heatmaps.

1. Compute similarity matrix with `Tools` > `Similarity Matrix` with options `"Metric"="Pearson correlation"` and `"Compute matrix for"="Columns"`
2. Perform `Hierarchical Clustering` with `Tools` > `Hierarchical Clustering` with options `"Metric": "Matrix values (for a precomputed similarity matrix)"`, `"Linkage method": "Average"`, `"Cluster": "Rows and Columns"`
3. Configure `Options` > `Annotations` > `Row annotations` to display `"Assay_Plate_Barcode"`, `"cell_line"`, `"diff_day"`, `"FFA"`, and `"patient"`
4. Configure `Options` > `Annotations` > `Column annotations` to display the same variables

## Output

We saved the heatmaps as `.png` and `.pdf` files (shown below)

### Plate 1 (15 Day)

![Plate 1](https://raw.githubusercontent.com/broadinstitute/2018_04_12_T2D_V2F_Saadat_Broad/master/figures/morpheus/BR00101075_normalized_variable_selected_2019_04_16_Batch1_morpheus.png)

### Plate 2 (15 Day + Iso)

![Plate 2](https://raw.githubusercontent.com/broadinstitute/2018_04_12_T2D_V2F_Saadat_Broad/master/figures/morpheus/BR00101076_normalized_variable_selected_2019_04_16_Batch1_morpheus.png)
