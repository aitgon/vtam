from datetime import datetime


class SummaryFileMerge:

    def __init__(self, params_dic, stats_df):

        self.stats_df = stats_df
        self.params_dic = params_dic

    def write(self, fname):

        title_str = "# Command Merge\n\n"
        date_str = datetime.now().strftime("%m/%d/%Y, %H:%M:%S\n\n")
        params_str = ""
        for k in self.params_dic:
            params_str = params_str + "{}: {}\n".format(k, self.params_dic[k])
        stats_str = self.stats_df.to_string(index=False, justify='left')

        out_str = title_str + date_str + params_str + stats_str + "\n\n\n"

        with open(fname, 'a') as fout:
            fout.write(out_str)
