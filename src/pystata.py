import os
import random
import re
import string
import numpy as np
from IPython.core.display import display, HTML
import pandas as pd
from src.pyreg import OLS

# import configure
import configparser

config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
config.read('config.ini')
stata_path =  config.get('Directory', 'StataDir')

class pystata(OLS):
    '''Tunnel between python and Stata '''
    def __init__(self, dat:pd.DataFrame, reg, cov_type, **kwargs):
        super().__init__(reg, cov_type, **kwargs)
        self.data = dat #python pd.DataFrame data
        self.inputDir = '' # Input Directory
        self.outputDir = '' # output Directory
        self.filename = 'temp.csv' # temporary filename for saved python pd.DataFrame data
        self.pydir = self.kwargs.get('dir')

        # initialize the pystata
        self.set_pydir(self.pydir[0], self.pydir[1], self.pydir[2])
        self._main_()

    def _main_(self):
        ''' Save data from python and read it on Stata'''
        all_vars = self._get_vars()
        self.filename = self._savedf_(self.data, all_vars)

        self.read_script = self._load_()
        self.do_script = self.read_script + self.do_script

    def set_pydir(self, workDir, _inputDir_, _outputDir_):
        ''' Set the working directory and input and output directory'''
        self.inputDir = workDir
        self.inputDir = _inputDir_
        self.outputDir = _outputDir_

    def _load_(self):
        ''' Load the data from python to Stata'''
        filename = f'"$input/temp_{self.filename}.csv"'
        comment = '* Load CSV saved from python \n'
        load_file = f'import delimited {filename}, clear \n'
        return comment + load_file

    def _get_vars(self):
        ''' Get the variables from the regression specification'''
        variables = self._parse(self.reg)
        unique_vars = list(set(variables))

        # get cluster and fixed effects variables
        try:
            cluster_vars = list(set(self.cluster_list))
        except TypeError:
            cluster_vars = []

        try:
            fe_vars = list(set(self.fixed_effects.values()))
        except AttributeError or TypeError:
            fe_vars = []
        all_var = unique_vars + cluster_vars + fe_vars

        return all_var

    def _save_(self):
        ''' Save the data in Stata'''
        filename = f'"$output/temp/temp.dta"'
        comment = '* Save the dta file for regression \n'
        save_line = f'save {filename}, replace \n'
        return comment + save_line

    def _savedf_(self, df, cols):
        ''' Save the python dataframe in a csv file'''
        filename = _randstr_()
        df[cols].to_csv(self.inputDir + f'temp_{filename}.csv', index=False)
        # print (f'Save pd.DataFrame to {self.inputDir}temp_{filename}.csv')
        return filename

    def save_dofile(self, _workDir_, _filename_):
        ''' Save the do file'''
        filename = _filename_ + '.do'
        with open(_workDir_ + filename, 'w') as f:
            for row in self.do:
                f.write(row)


class summary_col():
    '''
    Summary of the columns in the dataframe
    '''
    def __init__(self, reg_inputs: list):
        global stata_path
        self.reg_len = len(reg_inputs)
        self.reg_dict = dict(zip(range(0, self.reg_len), reg_inputs))
        self.data = reg_inputs
        self.workDir = ''
        self.inputDir = ''
        self.outputDir = ''
        self.do_script = ''
        self.macro = ''
        self.name = ''
        self.logDir = ''
        self.modelname = []
        self.order = []
        self.depvars = []
        self.indep_vars = []
        self.allvars = []
        self.fx = []
        self.outtype = 'tex'
        self.title = ''
        self.desc = ''
        self.fx_space = {}
        self.stata_path = stata_path

    def _main_(self):
        do_file = []

        # collect fixed effects info before running regression
        # for i, reg_input in enumerate(self.data):
        for reg_input in self.data:
            kwargs = self.get_fxcluster(reg_input)
            fx_temp = kwargs['fx']
            if fx_temp:
                self.fx_space.update(kwargs['fx'])


        # loop over each regression
        for i, reg_input in enumerate(self.data):
            kwargs = self.get_fxcluster(reg_input)
            kwargs['dir'] = [self.workDir, self.inputDir, self.outputDir]
            kwargs['column'] = i
            comments = self.comments(reg_input, **kwargs)
            kwargs['fx_space'] = list(self.fx_space.keys())
            if isinstance(reg_input[2], list):
                temp= pystata(reg_input[0], reg_input[1], 'clustered', **kwargs)
            else:
                temp= pystata(reg_input[0], reg_input[1], 'robust', **kwargs)

            temp.do_script = comments + temp.do_script
            do_file.append(temp.do_script)
            self.depvars.append(temp.dep_var)
            self.indep_vars.extend(temp.indep_vars)

        self.allvars = list(set(self.depvars + self.indep_vars))
        self.do_script = self.macro + '\n'.join(do_file)
        self.add_header
        self._summary_()
        self.save_dofile(self.outputDir, self.name)
        self.results = ''

    # generate comments
    def comments(self, reg, **kwargs):
        n = 300
        star_line = '*' * n + '\n'
        data = get_df_name(reg[0])
        reg = reg[1].strip()
        try:
            cov_type = reg[2]
        except IndexError:
            cov_type = 'Standard'
        fx = kwargs.get('fx')
        cluster = kwargs.get('cluster')

        firstline_temp = f' Data: {data} ' + '*' * 4 + f' Regression: {reg} '
        secondline_temp = f' Cov_type: {cov_type} ' + '*' * 4 + f' Cluster: {cluster} ' \
                          + '*' * 4 + f' Fixed Effects: {fx}'
        firstline = cfill(firstline_temp, n) + '\n'
        secondline = cfill(secondline_temp, n) + '\n'

        comment = star_line * 2 + firstline + secondline + star_line * 2
        return '\n' + comment

    # get cluster and fixed-effects
    def get_fxcluster(self, _input_):
        kwargs = {}
        if len(_input_) <= 2:
            return kwargs

        if len(_input_) == 3:
            if type(_input_[2]) is list:
                kwargs['cluster'] = _input_[2]
            elif type(_input_[2]) is dict:
                kwargs['fx'] = _input_[2]

        if len(_input_) == 4:
            kwargs['cluster'] = _input_[2]
            kwargs['fx'] = _input_[3]

        return kwargs

    @property
    def add_header(self):
        ''' Add header to the do file'''
        header = '*Please clean all the file in the output directory before running the code\n' \
                 'clear all\n' \
                 'macro drop _all\n' \
                 '* End of the header \n'
        self.do_script = header + self.do_script

    def set_dir(self, workDir, _inputsubDir_='input', _outputsubDir_='output', log='log'):
        ''' Set the directory for Stata'''
        globaldir = f'"{workDir}"'
        inputdir = f'"$dir/{_inputsubDir_}"'
        outputdir = f'"$dir/{_outputsubDir_}"'
        comment = '* Set global directory for code \n'
        setdir = f'global dir {globaldir} \n' \
                 f'global input {inputdir} \n' \
                 f'global output {outputdir} \n'

        self.workDir = workDir
        self.inputDir = folder_space(workDir, _inputsubDir_)
        self.outputDir = folder_space(workDir, _outputsubDir_)
        self.logDir = folder_space(workDir, log)
        self.macro = comment + setdir

    def _savedf_(self, df, cols):
        filename = _randstr_()
        df[cols].to_csv(self.inputDir + f'temp_{filename}.csv', index=False)
        return filename

    def save_dofile(self, _workDir_, _filename_):
        filename = _filename_ + '.do'
        with open(_workDir_ + filename, 'w') as f:
            for row in self.do_script:
                f.write(row)

    def _get_vars(self):
        variables = self._parse(self.reg)

        unique_vars = list(set(variables))
        try:
            cluster_vars = list(set(self.cluster_list))
        except TypeError:
            cluster_vars = []

        try:
            fe_vars = list(set(self.fixed_effects.values()))
        except AttributeError or TypeError:
            fe_vars = []
        all_var = unique_vars + cluster_vars + fe_vars


        return all_var


    def run_do(self):
        ''' Run the do file in Stata'''
        filename = self.outputDir + f'{self.name}'
        os.chdir(self.logDir)
        os.system(f'{self.stata_path} -e do {filename}')
        self.cleancache
        self.results = open(self.outputDir + f'{self.name}.{self.outtype}').read()

    def _summary_(self):
        ''' Stata summarize regression results'''
        models = [f'"{model}"' for model in self.modelname] if self.modelname else []
        modelname = ' mtitles (' + ' '.join(models) + ')' if models else ''

        order = [f'"{model}"' for model in self.order] if self.order else []
        orderlist = ' order (' + ' '.join(order) + ')' if models else ''

        reg_list = ' '.join(['reg' + str(i) for i in range(0, self.reg_len)])
        comment = '\n* Summary all regression in one table \n'

        _fx_stata = ' '.join([f'fx_{i}' for i in range(0, len(self.fx_space.keys()))])
        _fx_label = ' '.join([f'"{k}"' for k in self.fx_space.keys()])
        _format_stata = ' '.join(['0' for i in range(0, len(self.fx_space.keys()))])
        self.summary = f'esttab {reg_list} using "$output/{self.name}.{self.outtype}", replace ' \
                       'star(* 0.10 ** 0.05 *** 0.01) stat(' + _fx_stata + ' N r2_a tcov, ' \
                       'fmt('+ _format_stata + ' 0 3 0) label(' + _fx_label +  \
                       ' "Observations" "Adjusted R2" "SE Type")) noconstant' + modelname + orderlist.lower()
        self.do_script = self.do_script + comment + self.summary + '\n exit, STATA clear'

    def _readhtml_(self, columns=[]):
        ''' Read the html file'''
        self.results = open(self.outputDir + f'{self.name}.html').read()
        ind = self.results.find('td colspan=')
        ncol = findint(self.results[ind + 11:ind + 13])
        self.results = self.results.replace(f'<td colspan={ncol}><hr></td></tr>',
                                            f'<td colspan="{ncol}" style="border-bottom: 1px solid black"</td></tr>')
        # h = h.replace('colspan=10', 'colspan="10" style="border-bottom: 1px solid black"')
        if columns:
            self.rename(columns)
        display(HTML(self.results))

    def _readtex_(self, columns=[]):
        ''' Read the tex file'''
        self.results = open(self.outputDir + f'{self.name}.tex').read()
        if columns:
            self.rename(columns)
        # display(HTML(self.results))
        return self.results

    def rename(self, dep_dict):
        ''' Rename the independent variables'''
        lower_dict = dict((k.replace("_", "\_").lower(), v) for k, v in dep_dict.items())
        for oldname in lower_dict.keys():
            self.results = self.results.replace(oldname, lower_dict[oldname][0])


    def format(self, ratio=1):
        ''' Format the tex file, this is customized for my project (will be removed)'''
        table_start = self.results.find('\\begin{tabular}')
        self.results = self.results[:table_start] + '\n \\begin{adjustbox}{width=' + str(
            ratio) + '\\textwidth} \n' + self.results[table_start:]
        table_end = self.results.find('\\end{tabular}')
        self.results = self.results[:table_end + 13] + '\n \\end{adjustbox}' + self.results[table_end + 13:]

        title = self.title if self.title else 'Title'
        desc = self.desc if self.desc else 'Description'
        header = '\\begin{table} \n \\begin{center} \n'
        caption = '\\caption{{\\bf ' + title + ' }' + desc + '}' + '\\label{' + self.name + '} \n'
        footer = ' \\end{center} \n \\end{table}'

        self.results = header + self.results + caption + footer


    def add_avg(self, avg_list, name='Average'):
        ''' Add mean in first row, this is customized for my project (will be removed)'''
        ind = self.results.find('\\hline\\hline') + 14
        ind_add = self.results[ind::].find('\\hline') + 8
        avg = '&' + '&'.join([str(x) for x in avg_list])
        panela_head = '\\textit{Panel A: ' + name + '} \\\ \n'
        avg_line = avg + '\\\ \\\ \n'
        panela_foot = ' \\textit{Panel B: Estimate effect} \\\ \n'
        self.results = self.results[:ind + ind_add - 1] + panela_head + avg_line + panela_foot + self.results[
                                                                                                 ind + ind_add - 1::]

    def add_group(self, group):
        ''' add subgroup for columns, this is customized for my project (will be removed)'''
        tab_ind = self.results.find('\\begin{tabular}{') + 16
        tab_mod = '@{\extracolsep{4pt}}'
        self.results = self.results[:tab_ind] + tab_mod + self.results[tab_ind::]

        ind = self.results.find('\\hline\\hline') + 14
        tex = []
        cline = []
        start = 2
        for key in list(group.keys()):
            temp = group[key]
            length = temp[1] - temp[0] + 1
            col = '\\multicolumn{' + str(length) + '}{c}{' + key + '}'
            tex.append(col)
            cline_temp = '\\cline{' + f'{str(start)}-{str(start + length - 1)}' + '}'
            start += length
            cline.append(cline_temp)

        group_line = '&' + '&'.join(tex) + '\\\ ' + ' '.join(cline)
        self.results = self.results[:ind] + group_line + self.results[ind::]

    @property
    def cleancache(self):
        ''' Clean the cache'''
        temp_files = [temp for temp in os.listdir(self.inputDir) if 'temp_' in temp]
        os.chdir(self.inputDir)
        for item in temp_files:
            try:
                os.remove(item)
            except OSError as detail:
                print(f'{detail} : {item}')

    def save(self, **kwargs):
        ''' Save the tex file'''
        tempdir = kwargs.get('outputDir')
        tempname = kwargs.get('name')
        outputDir = tempdir if tempdir else self.outputDir
        filename = tempname if tempname else self.name

        filename = f'{filename}.{self.outtype}'
        with open(outputDir + filename, 'w') as f:
            for row in self.results:
                f.write(row)

    def __repr__(self):
        return self.results

    def __str__(self):
        self.results = open(self.outputDir + f'{self.name}.tex').read()
        printout = easy_print(self.results)
        return printout

    def print(self):
        '''
        Print tex table before complie in latex
        It is useful in ipython but not in jupyter
        '''
        self.results = open(self.outputDir + f'{self.name}.tex').read()
        printout = getmodel(self.results)
        return printout


def folder_space(_workDir_: str, subfolder_name: str, local_folder: bool = True):
    '''
    :param name: folder name
    :return: folder path
    '''

    if _workDir_ == '':
        CurrDir = os.getcwd().replace('\\', '/') + '/'
        workDir = CurrDir if local_folder else CurrDir + 'New folder/'
    else:
        workDir = _workDir_ if _workDir_[-1] == '/' else _workDir_ + '/'


    filename = workDir + subfolder_name + '/'

    try:
        os.makedirs(filename)
        # print("Directory " , dirName ,  " Created ")
    except FileExistsError:
        # print("Directory " , name ,  " already exists")
        pass
    return filename


def get_df_name(df):
    name = [x for x in locals() if locals()[x] is df][0]
    return name


def cfill(line, N):
    fill = int(np.round((N - len(line)) / 2))
    line = '*' * fill + line + '*' * fill
    return line


def _randstr_(N=10):
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=N))


def _readhtml_(h):
    h = h.replace('<hr>', '')
    h = h.replace('colspan=10', 'colspan="10" style="border-bottom: 1px solid black"')
    display(HTML(h))
    return h


from itertools import islice


def sread(fname, nlines):
    try:
        with open(fname) as f:
            for line in islice(f, nlines):
                print(line)
    except:
        with open(fname, encoding="utf8") as f:
            for line in islice(f, nlines):
                print(line)


def _cut_text_(tex, keyword, head = True):
    if len(keyword) == 2:
        s_ind = tex.find(keyword[0])
        e_ind = tex[s_ind+ len(keyword[0])::].find(keyword[1]) + s_ind + len(keyword[0])
        if s_ind == -1:
            s_ind == 0
        if e_ind != -1:
            return tex[s_ind + len(keyword[0]):e_ind]
        else:
            return tex[s_ind + len(keyword[0])::]
    elif len(keyword) == 1:
        if head:
            # keyword is the head
            ind = tex.find(keyword[0])
            _len = len(keyword[0])
            if ind == -1:
                ind = 0
            return tex[ind + _len::]
        else:
            # keyword is the end
            ind = tex.find(keyword[0])
            if ind != -1:
                return tex[:ind]
            else:
                return tex


def _clean_coefficients(coefficients, val_star = True):
    temp = [c.replace(" ","") for c in coefficients]
    coeff_temp = [_cut_text_(c,['\\sym'],head = False) for c in temp]
    if val_star:
        stars = [_cut_text_(c,['\\sym{','}']) for c in temp]
        coeff = [f'{coeff_temp[i]}{stars[i]}' for i in range(len(coeff_temp))]
    else:
        coeff = coeff_temp
    return coeff


def easy_print(tex):
    # split the tex by lines
    lines = tex.split('\n')

    # from table start and end
    ind = [i for i, x in enumerate(lines) if x == '\\hline\\hline']
    tables = lines[ind[0]+1:ind[1]]
    return tables

def getmodel(tex):
    # split the tex by lines
    lines = tex.split('\n')

    # from table start and end
    ind = [i for i, x in enumerate(lines) if x == '\\hline\\hline']
    tables = lines[ind[0]+1:ind[1]]

    # find model title
    models_ind = [i for i, x in enumerate(tables) if x == '\\hline']
    temp = tables[:models_ind[0]]

    models_temp = [x.strip(' \\').split('&') for x in temp]
    model_number = [_cut_text_(x,['{(',')}'] )
                    for x in models_temp[0] if x != '']
    model_name = [_cut_text_(x,['c}{','}'] )
                    for x in models_temp[1] if x != '']
    models = [f'{model_name[i]} ({model_number[i]})' for i in range(len(model_name))]

    # get table info
    rows = [x.strip('\\\\') for x in tables[models_ind[0]:models_ind[1]] if x not in ['\\hline','[1em]']]
    coefficients_list_temp = [rows[i].split("&")[1::] for i in range(len(rows))]
    coefficients_list = [_clean_coefficients(c) for c in coefficients_list_temp]

    vars =  [rows[i].split("&")[0].strip(' \\') for i in range(len(rows))]


    ## Statistics
    rows = [x.strip('\\\\') for x in tables[models_ind[1]::] if x not in ['\\hline','[1em]']]
    stat_list_temp = [rows[i].split("&")[1::] for i in range(len(rows))]
    stat_list = [_clean_coefficients(c) for c in stat_list_temp]

    stats_name =  [rows[i].split("&")[0].strip(' \\') for i in range(len(rows))]

    # put together
    print = pd.DataFrame(coefficients_list, columns = models)
    print['Var'] = vars

    # add line between table and statistics
    s_line = {f'{c}':'---' for c in print.columns.tolist()}
    print = print.append(s_line, ignore_index = True)

    stat = pd.DataFrame(stat_list, columns = models)
    stat['Var'] = stats_name

    print = print.append(stat)
    print[['Var'] + models]

    return print[['Var'] + models]

def findint(s):
    return int(re.search(r'\d+', s).group())
