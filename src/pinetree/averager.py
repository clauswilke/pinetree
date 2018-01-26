import pandas
import math

my_df1 = pd.read_csv("run0_counts.tsv", sep="\t", names=['time', 'protein', 'count'])
my_df2 = pd.read_csv("run1_counts.tsv", sep="\t", names=['time', 'protein', 'count'])

all_df = pd.concat(['my_df1', 'my_df2'])

all_df['time'] = all_df['time'].apply(math.floor)

all_df = all_df.groupby(['time', 'protein'])

all_df.agg({'count':['mean', 'std']}).reset_index()

