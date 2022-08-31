# Pystata
- 1. Save python data to csv
- 2. Use python to write stata do file following regression specifications
- 3. Report regression results and re-read it into python 

# Important!: __You need to have Stata installed and Stata license__

# Install Stata library if you are estimating fixed effects (enter these command in Stata terminal)
- ssc install reghdfe, ftools, esout

# Specified the Stata path so that Python can find it
Please specified the Stata path in the config.ini. Note that do not put GUI Stata path here

----
## Example (See example.ipynb for more details)

```python
from src.pystata import summary_col

# some random combinations of fixed effects
fx_1 = {'Stock fixed effects': 'fx1', 'Year fixed effects': 'fx2'}
fx_2 = {'Stock fixed effects': 'fx1', 'Industry Fixed effects': 'fx3'}
fx_3 = {'Stock fixed effects': 'fx1', 'Year fixed effects': 'fx2', 'Industry Fixed effects': 'fx3'}
# Syntax: [data, regression specification, covariance type (enter cluster list),fixed effects]
reg_inputs = [[data, 'Y  ~ 1  + x1+ x2', 'covariance type', {Fixed Effects}],  # This is an example (column 1)
[data1, 'y  ~ 1  + x1+ x2 ', 'robust', fx_1],  # (column 2) 
[data2, 'y  ~ 1  + x1+ x2 + x3 + x4', 'robust', fx_2],  # (column 3)
[data2, 'y  ~ 1  + x1+ x2 + x3 ', 'robust', fx_2],  # (column 4)
[data1, 'y  ~ 1  + x1+ x2 + x4', ['fx1', 'fx2'], fx_2],  # (column 5)
[data2, 'y  ~ 1  + x1+ x2 + x3 + x4', ['fx1', 'fx2'], fx_3]  # (column 6)
]
outputDir = '/home/user/pystata'  # set the directory to save Stata output (log and results)
table = summary_col(reg_inputs)  # read regression specification
table.set_dir(outputDir)  # set the directory to save Stata output (log and results)
table.name = 'table_pystata'  # set the name of the table
table.modelname = ["Y1", "Y1", "Var", "Variable", "Model name", "Y", ]  # set the name for columns
table.order = ['x1', 'x2', 'x3', 'x4']  # Determine independent variables order
table._main_()  # transit data from python to Stata and write Stata do file accordingly
table.run_do()  # run Stata do file
```

