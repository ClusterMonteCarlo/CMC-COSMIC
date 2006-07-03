# This requires 3d_density.rb, 2d_density.rb and quick_density_plot.sh. It just ties them all together in a single script. 
# Usage: make_density_plot.sh <file_prefix> <snapshot #> <half-mass radius>

echo "Extracting 3D density from snapshot..."
3d_density.rb $1 $2 $3 
echo "Extracting 2D luminosity/density from snapshot..."
2d_density.rb $1 $2 $3
echo "Generating plot..."
quick_density_plot.sh $1 $2
