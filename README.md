# ROITumorDetect

## Code

Code part 1: spatial smoothing, collection of all variables needed for subsequent modeling.
Code part 2: modeling, selection of variables, prediction visualization, etc.

## Data
### Full samples (cell information)

Variable explanation:
"X" and "Y": x and y coordinates
"class": class of manually annotated region in which the cell nucleus lies
"area", "perimeter" and "orientation": area size, perimeter and orientation of cell nucleus
"_color_", "variance\__color_": mean and variance of _color_ intensity within cell nucleus
"_color_ 3", "_color_ 7var": mean and variance of _color_ intensity within cell nucleus and a number of (respectively three and seven) pixels of neighborhood around the cell nucleus

### GEOJSON

geojson files which outline per sample the region coordinates

### detailed grids

detailed 100x100 grid tessellations (.Rda) and which grid cells are empty ('which.na') (.Rda)

### cellcounts

cellcounts100: counts of different types of cells in the 100x100 grid (.Rda)
cellcounts_extrasplit_*: counts of different types of cells (split up further into alveolar and bronchial epithelium or healthy cells) in the * x * grid (.Rda)

### nx*

Following data for the * x * grid 
Data frame per sample of all relevant variables for modeling purposes ('complete_variables_nx_* ') (.Rda) (Further explained in the section "Generated variables resulting from aggregation ...")
Selected variables in the forward selection procedure ('sel.5var.list.tvntbin.* ') (.Rda)
Ranking of which variable fits best in a univariate model ('var.ranking.tvntbin* ') (.Rda)

### Error contribution

Contribution of each cell type to the error (and comparison to random distribution of non-tumor cells scenario) to assess whether certain cell types are more frequently associated with poorer/better prediction.

### Morisita-Horn

Each non-tumor cell types spatial association with tumor cells (and comparison to random distribution of non-tumor cells scenario) using Morisita-Horn index.

### Generated variables resulting from aggregation and Gaussian Process modelling

Variable explanation:
Variables as explained in "Full samples (cell information)" + 

"rm": real mean of variable
"sm": mean of smoothed estimate of variable
"sv": variance of smoothed estimate of variable
"gm": mean of gradient estimate of variable
"gv": variance of gradient estimate of variable
