import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import re


class DataFrame_parser(object):
    def __init__(self, df, col):
        self.df = df
        self.col = col


    ### FLAGGING FOR REPLACE

    def duplicates(self):
        ret = self.df[self.df.duplicated(subset=self.col, keep=False)].groupby(by=self.col)
        return ret

    def implement_duplicates(self):
        r = {}
        for name, group in self.duplicates():
            r[name] = ", ".join(group.index.values.tolist())
        s = pd.Series(data=r)
        return s

###############################

    def punctuation(self):
        ret = self.df[self.df[self.col].str.contains("_|\+|:")][self.col]
        comment = pd.Series(index=ret.index, data="P")
        return ret, comment

    # def plus(self):
    #     ret = self.df[self.df[self.col].str.contains("\+")][self.col]
    #     comment = pd.Series(index=ret.index, data="P")
    #     return ret, comment

    def bad_words(self):
        bad_words = ["[Pp]roteins",
        "[Ss]ubunits",
        "[Ww]ith",
        "[Ww]ithin",
        "[Oo]n",
        "[Aa]n",
        "[Tt]o",
        "[Ii]n",
        "[Ii]nvolved",
        "[Ff]or",
        "[Aa]nd",
        "[Nn]or",
        "[Bb]ut",
        "[Oo]r",
        "[Ss]o"]
        bad_words = ['\W{0}\W'.format(elem) for elem in bad_words]
        ret = self.df[self.df[self.col].str.contains("|".join(bad_words))]
        comment = pd.Series(index=ret.index, data="F")
        return ret, comment

    def bad_start(self):
        ret = self.df[self.df[self.col].str.contains("^[P|p]redicted|[P|p]robable|[P|]utative")]
        comment = pd.Series(index=ret.index, data="S")
        return ret, comment

    def pref_words(self):
        words = 'precursor', 'homolog', 'paralog', 'ortholog', 'gene'
        words = [' {0} '.format(elem) for elem in words]
        ret = self.df[self.df[self.col].str.contains("|".join(words), regex=True)]
        comment = pd.Series(index=ret.index, data='W')
        return ret, comment

    def compile_flags(self):
        ret_df = self.df[[self.col,'COMMENT']]
        for r, c in [self.punctuation(), self.bad_words(), self.bad_start(), self.pref_words()]:
            ret_df['COMMENT'] = c.combine(ret_df['COMMENT'], lambda c, r:str(c)+str(r))
        ret_df['COMMENT'] = ret_df["COMMENT"].str.replace("nan", '')
        return ret_df.replace('', np.nan, regex=True)


def plot_pie(func, savedname=False, title=False, legend=False):
    val = func.dropna()['COMMENT'].value_counts()
    threshold = val.sum()/50
    val['Else'] = val[val < threshold].sum()
    labels = val[val>threshold].index
    values = val[val>threshold].values
    matplotlib.rcParams['font.size'] = 12.0
    fig, ax = plt.subplots(figsize=(10,10))
    ax.pie(values,
          labels=labels,
          startangle=90,
          autopct='%1.1f%%')
    handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles, [legend[a] for a in labels], bbox_to_anchor=(1, 1),
    #        bbox_transform=plt.gcf().transFigure)
    if title:
        plt.title(title)
    if savedname:
        plt.savefig(savedname, bbox_inches='tight')
    plt.show()

### S - replace semicolon, L - lowercase start, T - trailing stop, C - other stop, R - excessive capitals

### AUTOMATED REPLACE

def semicolon(df): #replace semicolons with commas
    ret = df[df['NAME'].str.contains(";")]['NAME'].str.replace(";", ',')
    comment = pd.Series(index=ret.index, name='COMMENT', data="S")
    return ret, comment

def lowercase_start(df): #replace lowercase start with capital
    st_lower = df[df['NAME'].str[0].str.islower()]['NAME']
    st_lower = st_lower.mask(st_lower.str.contains(r'^[m|t|r|ss|ds][R|D]NA|^cAMP', regex=True)).dropna()
    ret = st_lower.str[0].str.upper() + st_lower.str[1:]
    comment = pd.Series(index=ret.index,name='COMMENT', data="L")
    return ret, comment

def trailing_stop(df): #remove trailing dots
    ret = df[df["NAME"].str.contains('\.$|,$|;$')]['NAME'].str[:-1]
    comment = pd.Series(index=ret.index, name='COMMENT', data="T")
    return ret, comment

def other_stop(df): #replace other dots with commas
    s = df[df['NAME'].str.contains("\.")]['NAME']
    s = s.mask(s.str.contains(r'\d\.\d|\.$')).dropna()
    ret = s.str.replace('.', ',')
    comment = pd.Series(index=ret.index, name='COMMENT', data="C")
    return ret, comment

def ii_to_II(df):
    ret = df[df.NAME.str.contains('ii[\W|b]')].NAME.str.replace('ii', 'II')
    comment = pd.Series(index=ret.index, name='COMMENT', data="I")
    return ret, comment

def run_rename(df):
    cocRe = re.compile
    skip_list = [r'\w*[A-Z]\w*[A-Z]\w*',
                r'C-[T|t]erminal',
                r'N-[T|t]erminal',
                r'^[A-Z]-\w+',
                r'Hippel\-Lindau',
                r'Willebrand',
                r'Kunitz',
                r'cocc',
                r'^[A-Z]\W?$',
                r'^\d*[A-Z]\d*\W?$',
                r'\d+',
                r'[Kk][Dd][Aa]',
                r'[Zz][Nn]\W',
                r'[Mm][Uu]\W',
                r'Fe\W',
                r'Hint',
                r'gpD',
                r'Rab'
                ]
    acronRegex = re.compile('|'.join(skip_list))
    ret = pd.Series()
    for sfam in df.itertuples():
        l = sfam.NAME.split()
        new_name = [l[0]]
        for word in l[1:]:
            if acronRegex.search(word):
                new_name.append(word)
            else:
                new_name.append(word.lower())
        ret[sfam.Index] = " ".join(new_name)
    ret = ret[ret != df.NAME]
    comment = pd.Series(index=ret.index,name='COMMENT', data="R")
    return ret, comment

def org_names(df):
    ret = pd.Series()
    cocRE = re.compile('cocc')
    for sfam in df[df.NAME.str.contains('cocc')].itertuples():
        l = sfam.NAME.split()
        new_name = [l[0]]
        for word in l[1:]:
            if cocRE.search(word):
                new_name.append(word[0].upper() + word[1:])
            else:
                new_name.append(word)
        ret[sfam.Index] = " ".join(new_name)
    ret = ret[ret != df[df.NAME.str.contains('cocc')].NAME]
    comment = pd.Series(index=ret.index,name='COMMENT', data="O")
    return ret, comment

def terminal(df):
    reg = re.compile('(\W)([a-z])(\W)')
    def repl(mat):
        return mat.group(1) + mat.group(2).upper() + mat.group(3)
    ret = df[df.NAME.str.contains(reg)].NAME.str.replace(reg, repl)
    comment = pd.Series(index=ret.index,name='COMMENT', data="N")
    return ret, comment

def implement_replacements(df):
    ret_df = df[['NAME','COMMENT']]
    ret_df['OLD_NAME'] = ret_df['NAME']
    for f in [run_rename, semicolon, lowercase_start, trailing_stop, other_stop, org_names, ii_to_II, terminal]:
        r, c = f(ret_df)
        ret_df['NAME'] = r.combine_first(ret_df['NAME'])
        ret_df['COMMENT'] = c.combine(ret_df['COMMENT'], lambda c, r:str(c)+str(r))
    ret_df['COMMENT'] = ret_df["COMMENT"].str.replace("nan", '')
    ret_df = ret_df.rename(columns={'NAME':'NEW_NAME'})[['OLD_NAME', 'COMMENT', 'NEW_NAME']]
    return ret_df.replace('', np.nan, regex=True)
