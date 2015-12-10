import

def merge_files(fns_in, fn_out):
    """
    writing all files fns_ins into fn_out
    """
    with open(fn_out, 'w') as fout:
        for fn_in in fns_in:
            with open(fn_in) as fin:
                for line in fin:
                    fout.write(line)


def merge_assoc0_scan(assoc0file):
    """
    merging associations files
    """
    fn_beta = glob.glob(assoc0file + '*.beta.matrix')
    fn_pv = glob.glob(assoc0file + '*.pv.matrix')
    merge_files(fn_beta, assoc0file + '.beta.matrix')
    merge_files(fn_pv, assoc0file + '.pv.matrix')


