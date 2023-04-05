from pathlib import Path
import logging
import itertools

logger = logging.getLogger('prepare_data_for_clustering')
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    "%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

# auxiliary functions for probability calculation

xs = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
ys = [0., 0.005, 0.01, 0.02, 0.05, 0.14, 0.25, 0.5, 1., 1., 1.]

def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def prob(x):
    """ Given a similarity measure, it computes the probability 
        that this interaction exists given that it has that
        similarity measure (data extrapolated from Yu et al.,
        2004).
    """
    if x >= 1.0:
        return 1.0
    if x <= 0.0:
        return 0.0
    for i, (xmin, xmax) in list(enumerate(pairwise(xs))):
        if xmin <= x <= xmax:
            break
    ymin, ymax = ys[i], ys[i+1]
    return ymin + (ymax - ymin) * (x - xmin) / (xmax - xmin)

def run(interaction_file, interolog_file, outfile):
    outfile = open(outfile, "w")
    logger.info(f"Loading interolog file {interolog_file}")
    with open(interolog_file) as interof:
        for line in interof:
            p1, p2, _, _, quality, *_ = line.strip().split()
            outfile.write(f"{p1}\t{p2}\t{prob(float(quality) / 100.0)}\n")
    logger.info(f"Loading interactions file {interaction_file}")
    with open(interaction_file) as interof:
        for line in interof:
            p1, p2, *_ = line.strip().split()
            outfile.write(f"{p1}\t{p2}\t1.0\n")
    
    outfile.close()
    logger.info("Done")
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="formats the raw output of blast into the homolog information for ProteoBOOSTER")
    parser.add_argument("interactions_file", help="file with experimental interactions")
    parser.add_argument("interologs_file", help="file with predicted interactions")
    parser.add_argument("output", help="path to write the graph")
    args = parser.parse_args()
    run(args.interactions_file, args.interologs_file, args.output)

