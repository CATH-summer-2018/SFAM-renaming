import pandas as pd
import numpy as np

class DataFrame_parser(object):
    def __init__(self, df):
        self.df = df

    ### AUTOMATED REPLACE

    def semicolon(self): #replace semicolons with commas
        ret = self.df[self.df['NAME'].str.contains(";")]['NAME'].str.replace(";", ',')
        comment = pd.Series(index=ret.index, name='COMMENT', data="S")
        return ret, comment

    def lowercase_start(self): #replace lowercase start with capital
        st_lower = self.df[self.df['NAME'].str[0].str.islower()]['NAME']
        st_lower = st_lower.mask(st_lower.str.contains(r'^[m|t|r|ss|ds][R|D]NA|^cAMP', regex=True)).dropna()
        ret = st_lower.str[0].str.upper() + st_lower.str[1:]
        comment = pd.Series(index=ret.index,name='COMMENT', data="L")
        return ret, comment

    def trailing_stop(self): #remove trailing dots
        ret = self.df[self.df["NAME"].str.contains('\.$|,$|;$')]['NAME'].str[:-1]
        comment = pd.Series(index=ret.index, name='COMMENT', data="T")
        return ret, comment

    def other_stop(self): #replace other dots with commas
        s = self.df[self.df['NAME'].str.contains("\.")]['NAME']
        s = s.mask(s.str.contains(r'\d\.\d|\.$')).dropna()
        ret = s.str.replace(".", ',')
        comment = pd.Series(index=ret.index, name='COMMENT', data="C")
        return ret, comment

    def implement_replacements(self): #combine replacements with
        ret_df = self.df[['NAME','COMMENT']]
        ret_df['NEW_NAME'] = ret_df['NAME']
        for r, c in [self.semicolon(), self.lowercase_start(), self.other_stop(), self.trailing_stop()]:
            ret_df['NEW_NAME'] = r.combine_first(ret_df['NEW_NAME'])
            ret_df['COMMENT'] = c.combine(ret_df['COMMENT'], lambda c, r:str(c)+str(r))
        ret_df['COMMENT'] = ret_df["COMMENT"].str.replace("nan", '')
        return ret_df.replace('', np.nan, regex=True)


    ### S - replace semicolon, L - lowercase start, T - trailing stop, C - other stop

    ### FLAGGING FOR REPLACE

    def duplicates(self):
        ret = self.df[self.df.duplicated(subset="NAME", keep=False)].groupby(by='NAME')
        return ret

    def implement_duplicates(self):
        r = {}
        for name, group in self.duplicates():
            r[name] = ", ".join(group.index.values.tolist())
        s = pd.Series(data=r)
        return s

    def underscore(self):
        ret = self.df[self.df['NAME'].str.contains("_")]['NAME']
        comment = pd.Series(index=ret.index, data="U")
        return ret, comment

    def plus(self):
        ret = self.df[self.df['NAME'].str.contains("\+")]['NAME']
        comment = pd.Series(index=ret.index, data="P")
        return ret, comment

    def bad_words(self):
        bad_words = ["proteins", "subunits", "with", "within", "on", "an", "to", "in", "involved", "for", "and", "nor", "but", "or", "so"]
        bad_words = [' {0} '.format(elem) for elem in bad_words]
        ret = self.df[self.df['NAME'].str.contains("|".join(bad_words), regex=True)]
        comment = pd.Series(index=ret.index, data="F")
        return ret, comment

    def bad_start(self):
        ret = self.df[self.df['NAME'].str.contains("^[P|p]redicted|[P|p]robable|[P|]utative")]
        comment = pd.Series(index=ret.index, data="S")
        return ret, comment

    def pref_words(self):
        words = 'precursor', 'homolog', 'paralog', 'ortholog', 'gene'
        words = [' {0} '.format(elem) for elem in words]
        ret = self.df[self.df['NAME'].str.contains("|".join(words), regex=True)]
        comment = pd.Series(index=ret.index, data='W')
        return ret, comment

    def compile_flags(self):
        ret_df = self.df[['NAME','COMMENT']]
        for r, c in [self.underscore(), self.plus(), self.bad_words(), self.bad_start(), self.pref_words()]:
            ret_df['COMMENT'] = c.combine(ret_df['COMMENT'], lambda c, r:str(c)+str(r))
        ret_df['COMMENT'] = ret_df["COMMENT"].str.replace("nan", '')
        return ret_df.replace('', np.nan, regex=True)
