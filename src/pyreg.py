import numpy as np

class OLS(object):
    def __init__(self, reg:str, cov_type=None, **kwargs):
        # regression specification in string
        self.reg = reg

        # covariance type, and cluster level if clusted
        self.cov_type = cov_type
        self.cluster_list = []

        # keywords arg including cluster level and fixed effects
        self.kwargs = kwargs
        self.fixed_effects = self.kwargs.get("fx")  # dict input
        self.fixed_effects_space = self.kwargs.get("fx_space")  # dict input
        self.fx_dict = {}

        # winsorzation: default 1%
        self.winsor_per = 1

        # Estimation method
        self.est_method = 'OLS'

        # initiate blank scripts
        self.read_script = ''
        self.winsor_script = ''
        self.reg_script = ''
        self.do_script = ''

        # fixed effect indicator
        if self.fixed_effects:
            self.fx_name = list(self.fixed_effects.keys())
            self.fx_var = list(self.fixed_effects.values())
        else:
            self.fx_name, self.fx = None, None
        # fixed effect reporting

        if self.fx_name:
            for i, _fx in enumerate(self.fixed_effects_space):
                self.fx_dict[f'fx_{i}'] = "Yes" if _fx in self.fx_name else "No"
        else:
            for i, _fx in enumerate(self.fixed_effects_space):
                self.fx_dict[f'fx_{i}'] = "No"

        # parse str regression
        self._parse(self.reg)

        # check regression type
        self._regtype_()
        self.column = self.kwargs.get("column")
        self._regdict_ = {'PanelIV': self._ivpanelreg_, 'Panel': self._panelreg_, 'IV': self._ivreg_,
                          'OLS': self._olsreg_}

        # run
        self._agg_()

    def __str__(self):
        return f'Regression specification: {self.reg}'

    def _parse(self, formula: str, ret = False) -> list:
        '''
        parse the str format regression to dependent variables and independent variables
        '''
        parts = formula.split('~')
        self.dep_var = parts[0].strip()
        self.indep_vars = [var.strip() for var in parts[1].split('+')]

        # check if instrument variables regression
        self.bool_iv = np.any(['(' in var for var in self.indep_vars])

        # add instrumented variable and instrumented variable if specified
        # remove instrument and instrumented variables from self.independent variables
        if self.bool_iv:
            ivar_term = [v for v in self.indep_vars if '(' in v][0]
            # remove from self.indep_vars
            self.indep_vars.remove(ivar_term)

            # set instrument and instrumented variables
            self.hatv = ivar_term.split('(')[0].strip()
            self.iv = ivar_term.split('(')[1].replace(')','').strip()

        # remove constant
        try:
            self.indep_vars.remove('1')
        except ValueError:
            pass

        # if return is a must, return based on regression type
        if ret:
            if self.bool_iv:
                return [self.dep_var] + self.indep_vars + [self.hatv, self.iv]
            else:
                return [self.dep_var] + self.indep_vars

    def _regtype_(self):
        '''
        identify regression type
        '''
        # check if panel regression
        self.boolpanel = True if self.fixed_effects else False

        # identify regression type
        if self.boolpanel & self.bool_iv:
            self.regtype = 'PanelIV'
        elif self.boolpanel & (~self.bool_iv):
            self.regtype = 'Panel'
        elif (~self.boolpanel) & self.bool_iv:
            self.regtype = 'IV'
        else:
            self.regtype = 'OLS'

    def _olsreg_(self, quietly=True):
        '''
        estimate OLS as default
        '''
        reg_temp = f'regress {self.dep_var} ' + ' '.join(self.indep_vars)

        line_header = 'quietly ' if quietly else ''
        comment_head = f'* OLS regression for {self.reg} \n'
        reg_line = line_header + reg_temp + '\n'
        comment_end = f'* End of regression \n'  # for {self.reg}

        reg_do = comment_head + reg_line + comment_end
        self.reg_script = reg_do

    def _panelreg_(self, quietly=True):
        '''
        estimate panel regression if fixed effects are detected
        '''
        reg_temp = f'reghdfe {self.dep_var} ' + ' '.join(self.indep_vars)
        ff_list = ' '.join(list(self.fixed_effects.values()))
        ff_temp = f'absorb({ff_list})'

        line_header = 'quietly ' if quietly else ''
        comment_head = f'* Panel regression for {self.reg} with fixed effect: {ff_list} \n'
        reg_line = f'{line_header} {reg_temp}, {ff_temp}' + '\n'
        comment_end = f'* End of regression \n'  # for {self.reg}

        reg_do = comment_head + reg_line + comment_end
        self.reg_script = reg_do

    def _ivpanelreg_(self, quietly=True):
        '''
        estimate IV regression if instrumented variable is detected
        '''
        reg_temp = f'ivreghdfe {self.dep_var} ' + f'({self.hatv}={self.iv})' + ' '.join(self.indep_vars)

        line_header = 'quietly ' if quietly else ''
        comment_head = f'* IV regression for {self.reg} \n'
        reg_line = line_header + reg_temp + '\n'
        comment_end = f'* End of regression \n'  # for {self.reg}

        reg_do = comment_head + reg_line + comment_end
        self.reg_script = reg_do
        self.est_method = '2SLS'


    def _ivreg_(self, quietly=True):
        '''
        estimate Panel IV regression if instrumented variable is detected
        '''
        reg_temp = f'ivreghdfe {self.dep_var} ' + f'({self.hatv}={self.iv})' + ' '.join(self.indep_vars)
        ff_list = ' '.join(list(self.fixed_effects.values()))
        ff_temp = f'absorb({ff_list})'

        line_header = 'quietly ' if quietly else ''
        comment_head = f'* IVPanel regression for {self.reg} with fixed effect: {ff_list} \n'
        reg_line = f'{line_header} {reg_temp}, {ff_temp}' + '\n'
        comment_end = f'* End of regression \n'  # for {self.reg}

        reg_do = comment_head + reg_line + comment_end
        self.reg_script = reg_do
        self.est_method = '2SLS'

    def _se_type(self):
        '''
        parse standard error type and generate corresponding do script
        '''
        if self.cov_type == 'robust':
            se_line = 'vce(robust)'
            egen_cluster = ''
        elif self.cov_type == 'clustered':
            # check cluster level
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

        # find insertion place for covariance type
        ind = self.reg_script.find('\n* End of regression')

        # check regression line if comma is included
        if ',' in self.reg_script:
            self.reg_script = self.reg_script[:ind] + f' {se_line}' + self.reg_script[ind::]
        else:
            self.reg_script = self.reg_script[:ind] + f', {se_line}' + self.reg_script[ind::]

        self.reg_script = egen_cluster + self.reg_script

    def winsorize(self, var:str):
        '''
        :param var: winsorized variable
        :param percentile: winsorized percentage
        '''
        percentile = self.winsor_per
        comment_head = f'* Winsor variable {var} by {percentile}% each tail \n'
        winsor_temp = f'drop if {var} == . \n' \
                      f'egen pLow = pctile({var}), p({percentile}) \n' \
                      f'egen pHigh = pctile({var}), p({100 - percentile}) \n' \
                      f'replace {var} = pLow if {var} <= pLow \n' \
                      f'replace {var} = pHigh if {var} >= pHigh \n' \
                      f'drop pLow pHigh \n'
        comment_end = f'* End of winsorization of variable\n \n'

        # do script for winsoriztion
        winsor_line = comment_head + winsor_temp + comment_end
        self.winsor_script = winsor_line

    def _add_read(self):
        '''
        this is redundent
        '''
        str_filename = f'"{self.data}"'
        comment_head = f'* Read the csv file from python\n'
        read = f'use {str_filename}, clear \n \n'
        self.read_script = comment_head + read

    def _agg_(self):
        '''
        aggregate all the function
        '''
        # identify the corresponding regression type and estimate accordingly
        regression = self._regdict_[self.regtype]
        regression()

        # identify the covariance type and add accordingly to do script
        if self.cov_type:
            self._se_type()

        # add winsoriztion to do script
        if self.kwargs.get('winsor'):
            var_win = self.kwargs.get('winsor')
        else:
            var_win = self.dep_var

        self.winsorize(var_win)

        # clean do script and save
        self.do_script = (self.winsor_script + self.reg_script).lower()
        self._store_()

    def _store_(self):
        '''
        Ask stata to store regression results in Stata before printout
        '''
        # store regression results
        store_line = f'estimates store reg{self.column} \n'
        # store fixed effects
        for key, value in self.fx_dict.items():
            store_line += f'estadd local {key} "{value}", replace \n'
        # store covariance type
        store_line += f'estadd local tcov "{self.cov_type}", replace \n'
        # store estimation type
        store_line += f'estadd local estmethod "{self.est_method}", replace \n'
        self.do_script = self.do_script + store_line

    def all_lower(reglist):
        return [x.lower() for x in reglist]
