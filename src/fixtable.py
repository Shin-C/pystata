import numpy as np


class fixtable:
    def __init__(self, object, fullsample, **kwargs):
        self.table = object.results
        self.outputDir = object.outputDir
        self.depvars = object.depvars
        self.indep_vars = object.indep_vars
        self.outtype = object.outtype
        self.name = object.name
        self.fx = object.fx
        self.model = self.depvars
        self.fullsample = fullsample
        self.scale = 0

        # Panel A title
        self.paname = 'Mean'
        self.groupname = 'Default Name'
        # self.group_ind = kwargs.get('group')
        self.modelname_map = kwargs.get('modelmap')
        self.indep_map = self.modelname_map
        # self.columns = kwargs.get('columns')
        self.title = ''
        self.desc = ''
        self.key_var = []
        self.subgroup = False

    def read(self, file):
        self.table = open(file).read()

    def rename(self):
        lower_dict = dict((k.replace("_", "\_").lower(), v) for k, v in self.modelname_map.items())
        for oldname in lower_dict.keys():
            old1 = '\n' + oldname + ' '
            old2 = '\n' + oldname + '&'
            self.table = self.table.replace(old1, '\n' + lower_dict[oldname][0] + ' ')
            self.table = self.table.replace(old2, '\n' + lower_dict[oldname][0] + '&')

    def format(self, ratio=1):
        self.table = self.table[:2] + '\\tablebodyfont \n' + self.table[2:]

        title = self.title if self.title else 'Title'
        desc = self.desc if self.desc else 'Description'
        header = '\\begin{table} \n \\begin{center} \n'
        # caption = '\\caption{{\\tablecaptionfont{\\bf ' + title + ' }}' + '\\label{' + self.name + '}} \n'
        caption = '\\caption{ ' + title + ' }' + '\\label{' + self.name + '} \n'
        description = '\\caption*{{\\scriptsize ' + desc + '}} \n'
        footer = ' \\end{center} \n \\end{table}'

        self.table = header + caption + description + self.table + footer

    def cal_avg(self, pre=False):
        pre_avg = []

        multi_sample = len(self.fullsample) > 1
        for i, dep in enumerate(self.depvars):
            if pre:
                if multi_sample:
                    temp = self.fullsample[i]
                else:
                    temp = self.fullsample[0]
                temp = temp.query('after == 0')[dep].mean()
                pre_avg.append(temp)
            else:
                if multi_sample:
                    temp = self.fullsample[i]
                else:
                    temp = self.fullsample[0]
                temp = temp[dep].mean()
                pre_avg.append(temp)
        preavg = np.round(pre_avg, 2)
        return preavg

    def add_avg(self, avg_list, name='Average'):
        ind = self.table.find('\\hline\\hline') + 14
        ind_add = self.table[ind::].find('\\hline') + 8
        avg = '&' + '&'.join([str(x) for x in avg_list])
        panela_head = '\\textit{Panel A: ' + name + '} \\\ \n'
        avg_line = avg + '\\\ \\\ \n'
        panela_foot = ' \\textit{Panel B: Estimate effect} \\\ \n'
        self.table = self.table[:ind + ind_add - 1] + panela_head + avg_line + panela_foot + self.table[
                                                                                             ind + ind_add - 1::]

    def add_group(self, group):
        tab_ind = self.table.find('\\begin{tabular}{') + 16
        tab_mod = '@{\extracolsep{4pt}}'
        self.table = self.table[:tab_ind] + tab_mod + self.table[tab_ind::]

        ind = self.table.find('\\hline\\hline') + 14
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
        self.table = self.table[:ind] + group_line + self.table[ind::]

    def save(self, **kwargs):
        tempdir = kwargs.get('outputDir')
        tempname = kwargs.get('name')
        outputDir = tempdir if tempdir else self.outputDir
        filename = tempname if tempname else self.name

        filename = f'{filename}.{self.outtype}'
        with open(outputDir + filename, 'w') as f:
            for row in self.table:
                f.write(row)

    def _control_desc(self):
        key_var = [self.indep_map[c][1] for c in self.key_var]

        if len(key_var) == 1:
            inde_desc = 'the variable of interest is ' + key_var[0]
        else:
            inde_desc = 'the variable of interest is ' + ', '.join(key_var[0:-1]) + ' and ' + key_var[-1] + '.'

        control = unilist(self.indep_vars)
        try:
            control.remove(self.key_var)
        except ValueError:
            pass

        control_desc = [self.indep_map[c][1] for c in control]
        self.control_desc = 'For independent variable(s), ' + inde_desc + 'Control variables include ' + ', '.join(
            control_desc) + '. '

    def model_name(self):
        try:
            self.model = [self.modelname_map.get(c)[0] for c in self.depvars]
        except TypeError:
            pass

    def _dep_desc(self):
        uni_dep = unilist(self.depvars)
        dep_desc = [self.modelname_map.get(c) for c in uni_dep]
        depvar_desc = [f'{desc[0]} is {desc[1]}' for desc in dep_desc]
        self.dep_desc = 'For dependent variable(s), ' + '. '.join(depvar_desc) + '. '

    def titile(self):
        self.table

    @property
    def add_footer(self):
        desc_footer = 'Table A.1 in the online appendix defines the variables. ' \
                      't-statistics from robust standard errors clustered by date and stock are reported in the parentheses, ' \
                      'and *,** and *** indicate statistical significance at 10\%, 5\% and 1\% level, respectively.'
        self.desc = self.desc + desc_footer

    @property
    def add_header(self):
        title = self.title.split(' and')[0]
        # desc_header = 'This table presents coefficient estimates from fixed effects OLS regression of dependent variable ' \
        #               'on the independent variables. '
        desc_header = f'This table presents coefficient estimates from fixed effects OLS regression of {title}. '
        self.desc = desc_header + self.desc

    @property
    def remove_b_pvalue(self):
        self.table = self.table.replace(
            "\\multicolumn{6}{l}{\\footnotesize \\sym{*} \\(p<0.10\\), \\sym{**} \\(p<0.05\\), \\sym{***} \\(p<0.01\\)}\\\\\n",
            "")

    def main(self):
        # rename indepent
        self.rename()
        # set caption
        self._control_desc()
        self._dep_desc()

        self.desc = self.dep_desc + self.control_desc
        self.add_footer
        self.add_header

        # set the table scale
        if self.scale == 0:
            if len(self.depvars) >= 8:
                self.format(1.1)
            elif len(self.depvars) >= 6:
                self.format(0.85)
            elif len(self.depvars) >= 4:
                self.format(0.75)
            else:
                self.format(0.7)
        else:
            self.format(self.scale)

        # preavg = self.cal_avg()
        # self.add_avg(preavg, name=self.paname)

        if self.subgroup:
            group_ind = self.group_ind if self.group_ind else int(len(self.depvars) / 3 * 2)
            # self.add_group(
            #     group={'SMS Increase': [1, group_ind[0]], 'Liq $\\xrightarrow$ Illiq': [group_ind[0] + 1, group_ind[1]],
            #                 'Illiq $\\xrightarrow$ Liq': [group_ind[1] + 1, len(self.depvars)]})
            self.add_group(self._add_group_tool())
        else:
            group_ind = int(len(self.depvars))
            self.add_group(
                group={self.groupname: [1, group_ind]})

        self.remove_b_pvalue
        self.save()

    def _add_group_tool(self):
        _len = len(self.groupname)
        _len_ = len(self.group_ind) + 1
        assert _len == _len_

        inner_ind = [[self.group_ind[i] + 1, self.group_ind[i + 1]] for i in range(_len_ - 2)]
        full_ind = [[1, self.group_ind[0]]] + inner_ind + [[self.group_ind[-1] + 1, len(self.depvars)]]

        return {self.groupname[i]: full_ind[i] for i in range(_len_)}


def unilist(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
