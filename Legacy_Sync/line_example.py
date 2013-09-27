# -*- coding: utf-8 -*-
'''
Prints an example line with sync info describing a single event.
USAGE: python line_example.py
Christopher Palmer, Rishi Patel, Jan Veverka, 24 September 2013.
'''
import math

## These are without a prefix.
per_event_variables = '''
    run
    lumi
    event
    rho
    '''.split()

## These are without a prefix.
per_diphoton_evariables = '''
    mass
    rVtxSigmaMoM
    wVtxSigmaMoM
    vtxIndex
    vtxProb
    ptBal
    ptAsym
    logSPt2
    p2Conv
    nConv
    cosDPhi
    rho
    '''.split()

## Will be prepended with pho{1,2}_ for the leading and subleading photon.
per_photon_variables = '''
    ind
    scInd
    pt
    eta
    phi
    e
    eErr
    isConv
    HoE
    hcalIso03
    trkIso03
    pfChargedIsoGood02
    pfChargedIsoGood03
    pfChargedIsoBad03
    pfPhotonIso03
    pfNeutralIso03
    sieie
    sieip
    etaWidth
    phiWidth
    r9
    lambdaRatio
    s4Ratio
    scEta
    ESEffSigmaRR
    ptOverM
    '''.split()
leading_photon_variables  = ['pho1_' + var for var in per_photon_variables]
trailing_photon_variables = ['pho2_' + var for var in per_photon_variables]

## Will be prepended with "jet{1,2}_ for the {leading,subleading} jet.
per_jet_variables = '''
    ind
    pt
    eta
    '''.split()
leading_jet_variables  = ['jet1_' + var for var in per_jet_variables]
trailing_jet_variables = ['jet2_' + var for var in per_jet_variables]

## Will be prepended with "dijet_".
per_dijet_variables = '''
    dEta
    Zep
    dPhi
    Mjj
    '''.split()
per_dijet_variables = ['dijet_' + var for var in per_dijet_variables]

all_variables = (per_event_variables +
                 leading_photon_variables +
                 trailing_photon_variables +
                 per_diphoton_variables +
                 leading_jet_variables +
                 trailing_jet_variables +
                 per_dijet_variables)

## Maximum length of a variable name
max_var_len = max([len(var) for var in all_variables])

## Maximum length of a variable number (index)
max_num_len = int(math.floor(math.log10(len(all_variables)))) + 1


#_______________________________________________________________________________
def main():
    #print_example_line()
    #print ''
    #print_all_variables()
    print_all_cplusplus_snippets()
## End of main()


#_______________________________________________________________________________
def print_header():
    print '%*s   Variable' % (max_num_len, '#')
    line_len = max_num_len + 3 + max_var_len
    print ''.join(['-'] * line_len)


#_______________________________________________________________________________
def print_all_variables():
    print_header()
    print "Event Variables"
    print_variable_list(per_event_variables      )
    print "\nLeading Photon Variables"
    print_variable_list(leading_photon_variables )
    print "\nTrailing Photon Variables"
    print_variable_list(trailing_photon_variables)
    print "\nDiphoton Variables"
    print_variable_list(per_diphoton_variables   )
    print "\nLeading Jet Variables"
    print_variable_list(leading_jet_variables    )
    print "\nTrailing Jet Variables"
    print_variable_list(trailing_jet_variables   )
    print "\nDijet Variables"
    print_variable_list(per_dijet_variables      )
## End of print_all_variables()


#_______________________________________________________________________________
def print_all_cplusplus_snippets():
    print "// Event Variables"
    print_cplusplus_snippet(per_event_variables      )
    
    print "\n// Leading Photon Variables"
    print_cplusplus_snippet(leading_photon_variables )
    
    print "\n// Trailing Photon Variables"
    print_cplusplus_snippet(trailing_photon_variables)

    print "\n// Diphoton Variables"
    print_cplusplus_snippet(per_diphoton_variables)

    print "\n// Leading Jet Variables"
    print_cplusplus_snippet(leading_jet_variables    )
    
    print "\n// Trailing Jet Variables"
    print_cplusplus_snippet(trailing_jet_variables   )
    
    print "\n// Dijet Variables"
    print_cplusplus_snippet(per_dijet_variables      )
## End of print_all_variables()


#_______________________________________________________________________________
def print_variable_list(variables):
    for var in variables:
        var_num = all_variables.index(var) + 1
        print '%*d   %s' % (max_num_len, var_num, var)
## End of print_variable_list()


#_______________________________________________________________________________
def print_cplusplus_snippet(variables):
    for var in variables:
        var_num = all_variables.index(var) + 1
        mask = 'dumpVar("%s"%s, %s%s); // %*d'
        var_padding = ' ' * (max_var_len - len(var))
        args = (var, var_padding, var, var_padding, max_num_len, var_num)
        print mask % args
## End of print_cplusplus_snippet()


#_______________________________________________________________________________
def print_example_line():
    print '\t'.join([var + ':-999' for var in all_variables])
## End of print_example_line()

#_______________________________________________________________________________
if __name__ == '__main__':
    main()
    import user
    


