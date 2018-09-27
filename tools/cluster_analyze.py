import argparse
import json
#import plotly.plotly as py 
from plotly.offline import download_plotlyjs, plot, iplot
import plotly.figure_factory as ff 
import plotly.graph_objs as go 
import numpy as np 

def main():
    parser = argparse.ArgumentParser(description="Cluster Analyze")
    parser.add_argument('-g', '--graph', default=False, action="store_true", help='produce graphs of cluster distributions, etc. ')
    parser.add_argument('clusters', type=str, nargs='+', help='list of dirs containing cluster json')

    args = parser.parse_args()

    if (args.graph):
        for i, c in enumerate(args.clusters):
            with open(c) as f:
                cluster_json = json.load(f)
                # clusters.append(cluster_json)
            cluster_sizes = []
            for s in cluster_json:
                cluster_sizes.append(len(s))
            
            arr = np.asarray(cluster_sizes)

            # group_labels = ['cluster size']
            # hist_data = [arr]

            # fig = ff.create_distplot(hist_data, group_labels)
            data = [go.Histogram(x=arr)]
            plot(data, filename="cluster_sizes_{}.html".format(i))


    cluster_sets = []
    rep_sets = []

    sorted_lists = []

    for c in args.clusters:
        with open(c) as f:
            cluster_json = json.load(f)
        cs = set()
        reps = set()
        print("there are {} clusters in {}".format(len(cluster_json), c))
        for cluster in cluster_json:
            cs.add(frozenset(cluster))
            reps.add(cluster[0])  # rep is first seq
        print("there are {} cluster sets in {}".format(len(cs), c))
        cluster_sets.append(cs)
        rep_sets.append(reps)
        sorted_lists.append(sorted(cluster_json, key=lambda x: len(x)))

    final_set = cluster_sets[0]
    print("set 0 size is {}".format(len(final_set)))

    for i, c in enumerate(cluster_sets[1:]):
        print("set {} size is {}".format(i+1, len(c)))
        final_set = final_set.difference(c)

    print("final diff set size is {}".format(len(final_set)))

    rep_intersect = rep_sets[0]
    for i, r in enumerate(rep_sets[1:]):
        rep_intersect = rep_intersect.intersection(r)
    
    print("rep intersect set size is {}".format(len(rep_intersect)))

    for l in sorted_lists:
        print("biggest:\n")
        for i in range(5):
            print(l[-(i+1)])
            print("\n")





if __name__ == "__main__":
    main()