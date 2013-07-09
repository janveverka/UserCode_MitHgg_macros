'''
TTH Tag Performance Calculator

Reads output log file of a Hgg bambu run macro with TTH tag debug info
and calculates the efficiencies of the tth tag based on the
ttbar decay modes determined from the gen-level info:
  * fully hadronic
  * semi-leptonic
  * semi-leptonic (no taus)
  * semi-leptonic (1 tau)
  * leptonic
  * leptonic (no taus)
  * leptonic (1 tau)
  * leptonic (2 taus)
  
USAGE: python tth_tag_perf.py
'''

import commands
import math
import ROOT
ROOT.gROOT.SetBatch()

log_file = 'test100k.log'

alias_name_formula_map = dict(
    isQuark1 = 'abs(Wdau1)<10',
    isLepton1 = '(abs(Wdau1)==11 | abs(Wdau1)==13 | abs(Wdau1)==15)',
    isTau1 = 'abs(Wdau1)==15',
    isElemuon1 = '(abs(Wdau1)==11 | abs(Wdau1)==13)',
    )

## Add 1->2 symmetric aliases
for name1, formula1 in alias_name_formula_map.items():
    name2 = name1.replace('1', '2')
    formula2 = formula1.replace('Wdau1', 'Wdau2')
    alias_name_formula_map[name2] = formula2

category_name_selection_list = [
    ('All', ''),
    ('Hadronic'      , 'isQuark1 & isQuark2'),
    ('Semileptonic'  , '(isQuark1 & isLepton2) | (isLepton1 & isQuark2)'),
    ('Semileptonic (no tau)' , '(isElemuon1 & isQuark2) | (isQuark1 & isElemuon2)'),
    ('Semileptonic (w/ tau)'   , '(isTau1 & isQuark2) | (isQuark1 & isTau2)'),
    ('Leptonic'      , 'isLepton1 & isLepton2'),
    ('Leptonic (no taus)'     , 'isElemuon1 & isElemuon2'),
    ('Leptonic (1 tau)'     , '(isTau1 & isElemuon2) | (isElemuon1 & isTau2)'),
    ('Leptonic (2 taus)'       , 'isTau1 & isTau2'),
    ]

#______________________________________________________________________________
def main():
    global tree
    data_file = log_file + '.dat'
    create_data_file(source=log_file, destination=data_file)
    tree = get_tree(filename=data_file)
    print_yield_table(tree)
    print
    print_efficiency_table(tree)
## End of main


#______________________________________________________________________________
def create_data_file(source, destination):
    command = "grep WplusDau %(source)s" % locals()
    command += " | awk -F: '{print $3}' > %(destination)s" % locals()
    status, output = commands.getstatusoutput(command)
## End of create_data_file


#______________________________________________________________________________
def get_tree(filename):
    tree = ROOT.TTree('tree', 'tth tag info')
    tree.ReadFile(filename, 'tag/I:Wdau1:Wdau2')
    for name, formula in alias_name_formula_map.items():
        tree.SetAlias(name, formula)
    return tree
## End of get_tree(..)


#______________________________________________________________________________
def get_tag_counts(tree, selection):
    '''
    Calculates the counts of events with tthTag = 0, 1, and 2 
    '''
    ## Apply selection
    tree = tree.CopyTree(selection)
    if tree.Draw('', 'tag < 0 | tag > 2', 'goff') > 0:
        raise RuntimeError, 'Found untagged events!'
    return (tree.Draw('', 'tag==1', 'goff'),
            tree.Draw('', 'tag==2', 'goff'),
            tree.Draw('', 'tag==0', 'goff'))
## End of get_tag_counts(..)

#______________________________________________________________________________
def print_yield_table(tree):
    ## Print the header line
    headers = ('Type of ttbar decay', 'Lep. Tag', 'Had. Tag', 'No Tag', 'Total')
    print '%-21s   %-7s   %-7s   %-7s   %-7s' % headers
    ## Loop over various categories
    for name, selection in category_name_selection_list:
        counts = get_tag_counts(tree, selection)
        print format_line(name, counts, percent=False)


#______________________________________________________________________________
def print_efficiency_table(tree):
    headers = ('Type of ttbar decay', 'Lep. Tag (%)', 'Had. Tag (%)', 'No Tag (%)')
    print '%-21s   %-13s   %-13s   %-13s' % headers
    ## Loop over various categories
    for name, selection in category_name_selection_list:
        counts = get_tag_counts(tree, selection)
        print format_line(name, counts, percent=True)


#______________________________________________________________________________
def format_line(name, counts, percent=False):
    total = sum(counts)
    if percent:
        values = tuple(
            [name,] +
            [format_quantity(get_efficiency(x, total)) for x in counts]
            )
        line = '%-21s   %7s   %7s   %7s' % values
    else:
        values = (name,) + counts + (total,)
        line = '%-21s   %7d   %7d   %7d   %7d' % values
    return line
## End of format_line(name, counts)


#______________________________________________________________________________
def get_efficiency(selected, total, percent=True):
    value = float(selected) / total
    ## Normal approximation of binomial confidence interval
    error = math.sqrt(value * (1. - value) / total)
    if percent == True:
        value *= 100.
        error *= 100.
    return dict(value=value, error=error)
## End of get_efficiency(selected, total)


#______________________________________________________________________________
def format_quantity(x, precision=1):
    return '%*.*f +/- %*.*f' % (precision + 4, precision, x['value'], 
                                precision + 2, precision, x['error'])
## End of format_quantity(value, error, precision)


#______________________________________________________________________________
if __name__ == '__main__':
    main()
    import user

