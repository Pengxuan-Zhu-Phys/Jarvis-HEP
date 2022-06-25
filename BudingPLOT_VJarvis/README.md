# ![](src/image-src/BPicon.png)

					By: 	Buding

## Introduction 
An easy to use plot tool. 

`BudingPlot` using `.csv` format data, embedding 2D statistical method to show the `MultiNest` Data. 

Current Version Support:

​	1)	`2D_Stat_Profile` Method; 

​	2)	`2D_Scatter` Method; 

​	3)	`2DC_Scatter` Method;

​	4)	`1D_Stat` Method;

## Version Info

### V Pulsar

The biggest update since the publication of Version 2. 

1.)	Fixed bugs from user's report;

2.)	Support Chinese User Manual;

3.)	Redesigned the Logo.

### V2.4.4

1.)	Support `&pwd` method the option`[PLOT_CONFI]:path` . In this `&pwd` method, the option will automatically set the `[PLOT_CONFI]:path` as the directory of the plotting configuration file `*.ini`.

2.)	Support multi-format in the plotting. For example, if you want get `.pdf`, `.png`, `.jpg` and `.eps` format files in the your work, you could set the option `[PLOT_CONFI]:save_format`  by 
`save_format	= pdf, png, jpg, eps`.
Current supported version of `BudingPLOT` is listed as:
`['ps', 'eps', 'pdf', 'pgf', 'png', 'raw', 'rgba', 'svg', 'svgz', 'jpg', 'jpeg', 'tif', 'tiff']`.
If `save_format` is not declared in your `.ini` plotting configuration file, the default `.pdf` will be saved.

### V2.4.1

1.)	Fix bugs when the data is empty!

2.)	Change the path dependency of `COLORMAP` section, the options `colorsetting` and `StyleSetting` do not dependent on the path you declare in the option `[PLOT_CONFI] path`

### V2.4.0
1.)	Add `1D_Stat` Method, one dimensional PL and PDF distribution plots。

2.)	Fix bugs

### V2.3.0
1.)	Add `2DC_Scatter` Method, 2D Scatter with color bar;

2.)	`Manu` method for coordinates label is support.

### V2.2.0
1.)	Add Line Method and Text Method;

2.)	Add Function method `FUNCTION1D` . 	

### V2.1.1
1.)	fix bugs for the $1\sigma$ or $2\sigma$ curves when using `2D_Stat_Profile` method.

2.)	Support `Log` scale for `2D_Stat_Profile` method. 

### V2.1.0
1.)	fix bugs

2.)	Support `2D_Scstter` Mode. 

3.)	Standard `*.csv` format Data file now is supported.  format `*.dat` is not supported. 

4.)	Standard Scatter figure format `Default_Scatter` now is online. Color Setting of this format can be seen in file `src/preference.info` and the corresponding color can be seen in the `src/image-src/Default_Scatter_Format.png` . Any question please contact me first time!!!!

## Copyright Declaration 
Final interpretation right and copyright belongs to Author, if this code used in your work, please contact with me first!

2019-11-13

