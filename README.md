# ROITumorDetect

## Full samples (cell information)

Variable explanation:
"X" and "Y": x and y coordinates
"class": class of manually annotated region in which the cell nucleus lies
"area", "perimeter" and "orientation": area size, perimeter and orientation of cell nucleus
"_color_", "variance\__color_": mean and variance of _color_ intensity within cell nucleus
"_color_ 3", "_color_ 7var": mean and variance of _color_ intensity within cell nucleus and a number of (respectively three and seven) pixels of neighborhood around the cell nucleus

## Generated variables resulting from aggregation and Gaussian Process modelling

Variable explanation:
Variables as explained in "Full samples (cell information)" + 

"rm": real mean of variable
"sm": mean of smoothed estimate of variable
"sv": variance of smoothed estimate of variable
"gm": mean of gradient estimate of variable
"gv": variance of gradient estimate of variable
