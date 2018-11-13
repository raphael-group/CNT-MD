#!/usr/bin/python2.7

import os, sys
import math
from collections import deque
import argparse

dec2 = (lambda x : int(round((x * 100.0))) / 100.0)
dec4 = (lambda x : int(round((x * 10000.0))) / 10000.0)

order = (lambda y : int(''.join([l for l in y if l.isdigit()])))


def parse_args():
    """
    Parse command line arguments
    Returns:
    """
    description = "The following command allows to explicitly compute the fractional copy numbers from standard log2Ratio values of multiple samples and to output (in standard output) the corresponding input to execute CNT-MD. The procedure is similar to the one used by current methods, as ABSOLUTE and Battenberg; as such, this method aims to find the best values of tumor purity and tumor ploidy by assuming that solutions with more clonal CNAs are more likely and by assuming to know the most common total copy number in the genome of tumor cells. The required input is a tab-separated file in the format 'SAMPLE\tCHR\tSTART\tEND\t...\tLOG2RATIO'. The command also requires to provide the most common total copy number in the genome of tumor cells; this value is typically equal to 2 for diploid tumors that did not undergo whole-genome duplications, and is equal to 4 when a single WGD occurs. The commands does not automatically predicts the presence of a WGD."
    parser = argparse.ArgumentParser(description=description)#, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("INPUT", type=str, help="Input TAB-separated file in the format 'SAMPLE\tCHR\tSTART\tEND\t...\tLOG2RATIO'")
    parser.add_argument("-c", "--copynumber", required=False, type=int, default=2, help="Most-prominent copy number for tumor cells (integer, 2 by default for normal diploid segments, e.g. 2 is for diploid tumors and 4 for tetraploid tumors with a WGD occurred)")
    parser.add_argument("-b", "--bins", required=False, type=int, default=30, help="Number of bins user to discretize the distribution of logRatios and to fine the highest peak (default: 30, value can be increased to improve accuracy when the distribution is dense)")
    parser.add_argument("-m", "--maximumcn", required=False, type=int, default=12, help="Maximum copy number used to fit the logRatios (default: 12)")
    parser.add_argument("-u", "--minimumu", required=False, type=float, default=0.3, help="Minimum tumor purity used to fit the logRatios (default: 0.3)")
    args = parser.parse_args()

    if not os.path.isfile(args.INPUT):
        sys.stderr.write("ERROR: input file does not exist!\n")
        sys.exit(-1)
    if not args.copynumber > 0:
        sys.stderr.write("ERROR: most prominent copy number must be strictly greater than zero!\n")
        sys.exit(-1)
    if not args.bins > 0:
        sys.stderr.write("ERROR: number of bins must be strictly greater than zero!\n")
        sys.exit(-1)
    if not args.maximumcn >= args.copynumber:
        sys.stderr.write("ERROR: maximum copy number must be at least equal to the base copy number!\n")
        sys.exit(-1)
    if not 0.0 <= args.minimumu <= 1.00:
        sys.stderr.write("ERROR: minimum tumor purity must be within [0, 1]!\n")
        sys.exit(-1)


    return {
        "input" : args.INPUT,
        "base" : args.copynumber,
        "bins" : args.bins,
        "maxcn" : args.maximumcn,
        "minu" : args.minimumu
    }


def read_input(pt):
    read = {}
    with open(pt, 'r') as i:
        i.readline()
        for l in i:
            p = l.strip().split()
            if len(l) > 1:
                sample = p[0]
                chro = p[1]
                start = int(p[2])
                end = int(p[3])
                ratio = math.pow(2.0, float(p[-1]))

                if sample not in read:
                    read[sample] = {}
                if chro not in read[sample]:
                    read[sample][chro] = {}
                assert (start, end) not in read[sample][chro]
                read[sample][chro][(start, end)] = ratio

    return read


def joint_seg(d):
    samples = set(d.keys())
    chros = set(c for s in d for c in d[s])
    bks = {c : sorted(set(b for s in samples for e in d[s][c] for b in e)) for c in chros}

    counts = {c : {b : 0 for b in zip(bks[c][:-1], bks[c][1:])} for c in bks}
    mapb = {s : {c : {b : None for b in counts[c]} for c in counts} for s in samples}

    for s in samples:
        for c in counts:
            bk = deque(bks[c])
            left = -1
            right = bk.popleft()
            for (l, r) in sorted(d[s][c], key=(lambda x : x[0])):
                while right != r:
                    left = right
                    right = bk.popleft()
                    if l <= left and right <= r:
                        counts[c][left, right] += 1
                        assert counts[c][left, right] <= len(samples)
                        mapb[s][c][left, right] = (l, r)

    taken = {c : set(b for b in counts[c] if counts[c][b] == len(samples)) for c in counts}
    before = sum(float(s[1] - s[0]) for c in counts for s in counts[c])
    now = sum(float(s[1] - s[0]) for c in taken for s in taken[c])
    sys.stderr.write('## Proportion of joint covered genome from joint segmentation = {}%\n'.format(now / before * 100))

    res = {s : {c : {b : d[s][c][mapb[s][c][b]] for b in taken[c]} for c in taken} for s in samples}

    return res


def get_fractions(d, bins, base, maxcn, minu):
    wins = {s : sorted([d[s][c][e] for c in d[s] for e in d[s][c]]) for s in d}
    perc = {s : int(round(float(len(wins[s])) * 0.1)) for s in d}
    wins = {s : wins[s][perc[s]:-perc[s]] for s in d}
    minr = {s : min(wins[s]) for s in d}
    maxr = {s : max(wins[s]) for s in d}

    sbin = {s : (maxr[s] - minr[s]) / bins for s in d}
    lb = (lambda s, i : minr[s] + i * sbin[s])
    ub = (lambda s, i : minr[s] + (i + 1) * sbin[s])
    p = (lambda s, i : (lb(s, i), ub(s, i)))

    check = (lambda s, i, x : max(minr[s], lb(s, i)) <= x <= min(ub(s, i), maxr[s]))
    ppeak = (lambda s, i : sum(e[1] - e[0] for c in d[s] for e in d[s][c] if check(s, i, d[s][c][e])))
    peak = {s : {p(s, i) : ppeak(s, i) for i in range(bins)} for s in d}
    top = {s : max(peak[s].keys(), key=(lambda x : peak[s])) for s in d}
    top = {s : (top[s][1] + top[s][0]) / 2.0 for s in d}

#         totsize = {s : sum(e[1] - e[0] for c in d[s] for e in d[s][c]) for s in d}
#         top = {s : sum((e[1] - e[0]) * d[s][c][e] for c in d[s] for e in d[s][c]) / totsize[s] for s in d}

    scale = None
    obj_scale = None
    for i in range(100, int(round(minu*100)), -1):
        pur = float(i) / 100.0
        center = (lambda y : 2.0 * (1 - pur) + y * pur)
        cbase = center(base)
        gamma = {s : cbase / top[s] for s in d}
        dist = (lambda x : min(abs(x - center(i)) / x for i in range(0, maxcn + 1)))
        totdist = (lambda s : sum(dist(d[s][c][e] * gamma[s]) * (e[1] - e[0]) for c in d[s] for e in d[s][c]))
        obj = {s : totdist(s) for s in d}

        if scale is None:
            scale = gamma
            obj_scale = obj
        else:
            for s in d:
                if obj[s] < (obj_scale[s] - 0.0001):
                    scale[s] = gamma[s]
                    obj_scale[s] = obj[s]

    return {s : {c : {e : dec2(d[s][c][e] * scale[s]) for e in d[s][c]} for c in d[s]} for s in d}


def write_input(frac):
    print '#PARAMS'
    print '{} #number of chromosomes'.format(len(set(c for s in frac for c in frac[s])))
    print '{} #number of samples'.format(len(frac.keys()))
    tg = frac.keys()[0]
    print '{} #number of segments for each chromosome'.format(' '.join([str(len(frac[tg][c])) for c in frac[tg]]))
    print '#SAMPLES',
    lbe = (lambda c : ' '.join(['{},{}'.format(e[0], e[1]) for e in sorted(frac[tg][c], key=(lambda x : x[0]))]))
    print ' | '.join(['{} : {}'.format(c, lbe(c)) for c in sorted(frac[tg], key=order)])

    fmc = (lambda s, c : ' '.join([str(frac[s][c][e]) for e in sorted(frac[s][c], key=(lambda x : x[0]))]))
    fms = (lambda s : ' | '.join([fmc(s, c) for c in sorted(frac[s], key=order)]))
    print '\n'.join(['{}-{} : {}'.format(i, s.replace(' ', '_'), fms(s)) for i, s in enumerate(sorted(frac.keys()))])


def main():
    sys.stderr.write("# Reading and parsing arguments\n")
    args = parse_args()

    sys.stderr.write("# Reading input log ratios\n")
    rdr = read_input(args['input'])
    sys.stderr.write('## {} samples found\n'.format(len(rdr)))

    sys.stderr.write("# Performing joint segmentation of all samples\n")
    join = joint_seg(rdr)

    sys.stderr.write("# Computing fractional copy numbers\n")
    frac = get_fractions(d=join, bins=args['bins'], base=args['base'], maxcn=args['maxcn'], minu=args['minu'])

    sys.stderr.write("# Writing prepared input for CNT-MD\n")
    write_input(frac)


if __name__ == "__main__":
    main()
