class OLS():
    ''' Write Stata do file for a given OLS regression specification '''
    def __init__(self, reg: str, cov_type=None, **kwargs):
        self.reg = reg # dependent and independent variables, e.g. 'y ~ x1 + x2 + x3'
        self.cov_type = cov_type # robust or clustered
        self.cluster_list = [] # list of clusters if cluster standard error is choosen, e.g. ['cluster1', 'cluster2']
        self.read_script = ''
        self.winsor_script = ''
        self.reg_script = ''
        self.do_script = ''
        self.kwargs = kwargs
        # dict input of fixed effects, e.g. {'stock fixed effects': 'x1', 'Date fixed effects': 'x2'}
        self.fixed_effects = self.kwargs.get("fx")
        self.fixed_effects_space = self.kwargs.get("fx_space")  # dict input: this is for multi-regression reporting
        self.fx_dict = {}

        # Get the name and group of fixed effects.
        if self.fixed_effects:
            self.fx_name = list(self.fixed_effects.keys())
            self.fx_var = list(self.fixed_effects.values())
        else:
            self.fx_name, self.fx = None, None

        # Text for fixed effects reported in regression table
        for i, _fx in enumerate(self.fixed_effects_space):
            self.fx_dict[f'fx_{i}'] = "Yes" if _fx in self.fx_name else "No"

        # check panel data
        self.boolpanel = True if self.fixed_effects else False

        # replace column names
        self.column = self.kwargs.get("column")

        # aggreate all parts and run
        self._agg_()

    def _parse(self, formula):
        ''' Parse the formula and return the dependent variable and independent variables '''
        parts = formula.split('~')
        dep_var = parts[0].strip()
        indep_vars = [var.strip() for var in parts[1].split('+')]
        try:
            indep_vars.remove('1')
        except Exception:
            pass
        return [dep_var] + indep_vars

    def _olsreg_(self, quietly=True):
        ''' Write the do file for the OLS regression '''
        self.dep_var, *indep_vars = self._parse(self.reg)
        self.indep_vars = indep_vars

        reg_temp = f'regress {self.dep_var} ' + ' '.join(indep_vars)

        line_header = 'quietly ' if quietly else ''
        comment_head = f'* OLS regression for {self.reg} \n'
        reg_line = line_header + reg_temp + '\n'
        comment_end = f'* End of regression \n'  # for {self.reg}

        reg_do = comment_head + reg_line + comment_end
        self.reg_script = reg_do

    def _panelreg_(self, quietly=True):
        ''' Write the do file for the panel regression '''
        self.dep_var, *indep_vars = self._parse(self.reg)
        self.indep_vars = indep_vars

        reg_temp = f'reghdfe {self.dep_var} ' + ' '.join(indep_vars)
        ff_list = ' '.join(list(self.fixed_effects.values()))
        ff_temp = f'absorb({ff_list})'

        line_header = 'quietly ' if quietly else ''
        comment_head = f'* Panel regression for {self.reg} with fixed effect: {ff_list} \n'
        reg_line = f'{line_header} {reg_temp}, {ff_temp}' + '\n'
        comment_end = f'* End of regression \n'  # for {self.reg}

        reg_do = comment_head + reg_line + comment_end
        self.reg_script = reg_do

    def _se_type(self):
        ''' Specifiy the standard error type '''
        if self.cov_type == 'robust':
            se_line = 'vce(robust)'
            egen_cluster = ''
        elif self.cov_type == 'clustered':
            cluster = self.kwargs.get("cluster")
            if len(cluster) == 1:
                egen_cluster = f'egen cluster_level1 = group({cluster})\n'
                self.cluster_list.append(cluster)
                se_line = 'vce(cluster  cluster_level1)'

            if len(cluster) == 2:
                double_cluster = ' '.join(cluster)
                egen_cluster = f'egen cluster_level2 = group({double_cluster})\n'
                self.cluster_list.extend(cluster)
                se_line = 'vce(cluster  cluster_level2)'

        ind = self.reg_script.find('\n* End of regression')

        if ',' in self.reg_script:
            self.reg_script = self.reg_script[:ind] + f' {se_line}' + self.reg_script[ind::]
        else:
            self.reg_script = self.reg_script[:ind] + f', {se_line}' + self.reg_script[ind::]

        self.reg_script = egen_cluster + self.reg_script


    def winsorize(self, var, percentile=1):
        ''' Winsorize the data: default by 1% '''
        comment_head = f'* Winsor variable {var} by {percentile}% each tail \n'
        winsor_temp = f'drop if {var} == . \n' \
                      f'egen pLow = pctile({var}), p({percentile}) \n' \
                      f'egen pHigh = pctile({var}), p({100 - percentile}) \n' \
                      f'replace {var} = pLow if {var} <= pLow \n' \
                      f'replace {var} = pHigh if {var} >= pHigh \n' \
                      f'drop pLow pHigh \n'
        comment_end = f'* End of winsorization of variable\n \n'

        winsor_line = comment_head + winsor_temp + comment_end
        self.winsor_script = winsor_line

    def _add_read(self):
        ''' Read data from python's output '''
        str_filename = f'"{self.data}"'
        comment_head = f'* Read the csv file from python\n'
        read = f'use {str_filename}, clear \n \n'
        self.read_script = comment_head + read

    def _agg_(self):
        ''' Aggregate all parts of the do file '''
        # check the type of OLS regression
        if self.boolpanel:
            self._panelreg_()
        else:
            self._olsreg_()

        # get covariance type
        if self.cov_type:
            self._se_type()

        # winsorization
        if self.kwargs.get('winsor'):
            var_win = self.kwargs.get('winsor')
        else:
            var_win = self.dep_var
        self.winsorize(var_win)

        # write the do file
        self.do_script = (self.winsor_script + self.reg_script).lower()
        self._store_()

    def _store_(self):
        ''' Store the regression result in do file '''
        store_line = f'estimates store reg{self.column} \n'
        for key, value in self.fx_dict.items():
            store_line += f'estadd local {key} "{value}", replace \n'
        store_line += f'estadd local tcov "{self.cov_type}", replace \n'
        self.do_script = self.do_script + store_line

    def all_lower(reglist):
        ''' Make all the regressor names lower case '''
        return [x.lower() for x in reglist]
