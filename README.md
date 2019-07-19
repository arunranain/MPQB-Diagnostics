# MPQB-Diagnostics
Diagnostics for MPQB Evaluation on CDS
 
The following scripts are meant for MPQB Diagnostics and plotting.
The main script (after remapping all the datasets) has been divided into 2 parts i.e spatial comparison diagnostics (averaged over time length) and temporal comparison diagnostics (averagered over space).
The below space is used to define what/where you as ECV evaluator has to change in order to run it for your ECV.
The following should be modified according to needs of your ECV:

1. Install skill_metrics module - $ pip install SkillMetrics
2. Input Data Information - provide links to your datasets that are interpolated with "Remap_v0.bash" or similar. You need to define Reference, Dataset1, Dataset2 and so on.
3. In all spatial plots please use your ECV specific/desired projection for mapping. In the following script we have used Mercator. Other Projection are listed and linked in documentation.
4. I will leave other lat/lon information same unless compelling to keep the plots homogenised and same goes for filling continents and coastlines.
5. In all spatial plots you need to adjust colormap according to ECV in cmap='yourcolor' in e.g. m.pcolormesh(lon, lat, PlotVariable,
             latlon=True, cmap='Blues_r')
6. Similarly you need to define the range for the displayed variable accoding to diagnostic presented by changing plt.clim range values in e.g. plt.clim(0, 100)
7. Please define the plot tiles as per need and convinience in all the plots.
8. Means and variability - Annual Mean - Position of the colorbar (cax1 = fig.add_axes([0.9, 0.2, 0.02, 0.6])) in the plot needs to be adjusted to linking and based on number of datasets in the plot.
9. Means and variability - grid-point correlation (and in all futher spatial plots) - Both position of the colorbars need to be fixed based on the plot/datasets in cax1 and cax2.
10. Time series of global/hemisphere annual/monthly mean (line plot) - Start and end year needs to be defined and changed. If daily dataset then we need to account for days as well. All the labels/titles/linecolor(ranges and position) needs to be defined.
